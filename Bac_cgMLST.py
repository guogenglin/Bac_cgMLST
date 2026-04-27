# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 08:46:19 2023

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

import re
import sys
import uuid
import gzip
import shutil
import pathlib
import tempfile
import argparse
import subprocess
from tqdm import tqdm
from Bio import SeqIO
import multiprocessing
from Bio.Seq import Seq
from multiprocessing import Pool

__version__ = '1.2'

def get_argument() -> argparse.ArgumentParser:
    # Parsers
    parser = argparse.ArgumentParser(description = 'Bac_cgMLST', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-r', '--reference', required = True, type = str, 
                                help = 'reference core genome file or path toprepared database')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'cgMLST_output.txt',
                              help = 'Output file')

    # Parameters
    parser_group_2.add_argument('-m', '--mode', required = False, type = str, default = 'core_genome', 
                                choices = ['core_genome', 'database'], 
                                help = 'Both core genome file and prepared database can be used as reference')
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'cgMLST v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def check_programs_shutil(progs: list[str]):
    '''
    Check if programs are installed and executable using shutil.which.
    '''
    missing = []
    for prog in progs:
        path = shutil.which(prog)
        if not path:
            missing.append(prog)
    if missing:
        print(f'Error: could not find {", ".join(missing)}', file = sys.stderr)
        sys.exit(1)

class cgMLSTError(RuntimeError):
    pass

def parse_database(inputfile: str) -> dict[str, str]:
    '''
    Parse the reference core genome file into a dictionary of {gene_name: sequence}.
    '''
    input_seq = {}
    for contig in SeqIO.parse(inputfile, 'fasta'):
        input_seq[contig.name] = contig.seq
    if not input_seq:
        print('invalid FASTA file: %s', inputfile)
        sys.exit(1)
    return input_seq

def separate_core(reference_cgMLST: str, workpath: pathlib.Path) -> tuple[pathlib.Path, list[str]]:
    '''
    Separate the core genome sequences into individual gene files and return the path and gene list.
    '''
    if not pathlib.Path('cgMLST').is_dir():
        pathlib.Path('cgMLST').mkdir()
    cgpath = workpath / 'cgMLST'
    cglist = []
    reference = parse_database(reference_cgMLST)
    for key, value in reference.items():
        gene_name = key
        cglist.append(gene_name)
        filepath = cgpath / f'{gene_name}.fasta'
        sequence = str(value)
        sequence_for_write = ''
        while len(sequence) > 60:
            sequence_for_write += sequence[:60] + '\n'
            sequence = sequence[60:]
        if sequence:
            sequence_for_write += sequence
        with open(filepath, 'wt') as file:
            file.write('>1')
            file.write('\n')
            file.write(sequence_for_write)
            file.write('\n')
    
    return cgpath, cglist

def get_compression_type(filename: pathlib.Path) -> str:
    '''
    Determine the compression type of a file based on its magic bytes. Supports gzip, bzip2, and zip.
    '''
    filename = pathlib.Path(filename)
    if not filename.is_file():
        raise cgMLSTError(f"{filename} does not exist or is not a file")
    magic_dict = {
        'gz': b'\x1f\x8b\x08',
        'bz2': b'BZh',
        'zip': b'PK\x03\x04',
    }

    max_len = max(len(m) for m in magic_dict.values())
    with filename.open('rb') as f:
        file_start = f.read(max_len)

    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            if file_type in ['bz2', 'zip']:
                raise cgMLSTError(f"{file_type} format is not supported; use gzip instead")
            return file_type

    return 'plain'

def gunzip_assembly(assembly: pathlib.Path, temp_dir: str) -> pathlib.Path:
    '''
    If the assembly is gzipped, gunzip it to a temporary file and yield the path. 
    Otherwise, yield the original path.
    '''
    temp_path = pathlib.Path(temp_dir)
    assembly = pathlib.Path(assembly)

    assembly_name = assembly.stem
    compression_type = get_compression_type(assembly)

    if compression_type == 'gz':
        unzipped_path = temp_path / f"{uuid.uuid4().hex}.fasta"

        with gzip.open(assembly, 'rb') as f_in, unzipped_path.open('wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        return unzipped_path
    
    return assembly

class Orf:
    # Parse the prodigal results
    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.sequence = str()

class BlastResult:
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0].split('_')[-1]
        self.sqseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4]) - 1
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])

class Assembly:
    '''
    Represents an assembly with multiple contigs. Provides methods to predict ORFs, create a BLAST database, 
    and map query sequences to the assembly.
    '''
    def __init__(self, path: pathlib.Path, name: str | None = None, contigs: dict[str, Seq] | None = None):
        self.path = pathlib.Path(path).resolve()
        self.name = name or self.path.stem
        self.contigs = contigs or {}
        self.blastdb_name = ''
        self.blast_hits = []
        self.fna = dict()
        self.orfs = []

    def __repr__(self) -> str:
        return f"<Assembly {self.name}, {len(self.contigs)} contigs>"

    def __len__(self) -> int:
        return sum(len(c) for c in self.contigs.values())

    def prodigal(self) -> list[Orf]:
        '''
        Predict ORFs in the assembly and return a list of Orf objects.
        '''
        command = 'prodigal -f sco -i %s -m' % self.path
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
        raw_orfs = result.stdout.decode().replace('\r', '')
        orfs = []
        contig = ''
        RESULT_RE = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-+])$')
        CONTIG_RE = re.compile(r'^# Sequence.+?seqhdr="(.+?)"(?:;|$)')
        for line in raw_orfs.rstrip().split('\n'):
            if line.startswith('# Sequence Data'):
                name = CONTIG_RE.match(line).group(1)
                contig = name.split(' ')[0]
            elif line.startswith('# Model Data'):
                continue
            else:
                result = RESULT_RE.match(line).groups()
                orfs.append(Orf(contig, *result))
        for orf in orfs:
            orf.sequence = self.contigs[orf.contig][orf.start-1:orf.end]

        self.orfs = orfs

        return self.orfs
    
    def make_input_blastdb(self, temp_dir: str) -> tuple[pathlib.Path, dict[str, str]]:
        '''
        Create a BLAST database from the predicted ORFs and return the database name and a dictionary of 
        {gene_code: sequence}.
        '''
        blastdb_name = pathlib.Path(temp_dir) / 'input_fna.fna'
        input_fna_dict = dict()
        input_fna = open(blastdb_name, 'wt')
        for orf in self.orfs:
            gene_code = str(orf.contig) + '_' + str(orf.start) + '_' + str(orf.end) + '_' + str(orf.strand)
            sequence = str(orf.sequence)
            input_fna_dict[gene_code] = sequence
            sequence_for_write = ''
            while len(sequence) > 60:
                sequence_for_write += sequence[:60] + '\n'
                sequence = sequence[60:]
            if sequence:
                sequence_for_write += sequence
            input_fna.write('>')
            input_fna.write(gene_code)
            input_fna.write('\n')
            input_fna.write(sequence_for_write)
            input_fna.write('\n')
        input_fna.close()
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', blastdb_name]
        makeblastdb_process = subprocess.run(makeblastdb_command, stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)
        self.fna = input_fna_dict
        self.blastdb_name = blastdb_name
        return self.blastdb_name, self.fna

    def map(self, query: pathlib.Path | str | None) -> list:
        '''
        Map a query sequence to the assembly using BLAST and return a list of BlastResult objects
        representing the hits.
        '''
        command = ['blastn', '-query', query, '-db', self.blastdb_name, '-num_threads', '1', '-outfmt', 
                '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
        result = subprocess.run(command, capture_output = True, text = True)

        self.blast_hits = [BlastResult(line) for line in result.stdout.splitlines() if line]
        return self.blast_hits
    
def parse_input_assembly(assembly: pathlib.Path) -> Assembly:
    '''
    Parse the input assembly FASTA file into an Assembly object.
    '''
    basename = assembly.stem

    # Parse fasta file
    input_seq = {}

    try:
        for contig in SeqIO.parse(assembly, 'fasta'):
            input_seq[contig.name] = contig.seq

    except Exception as e:
        raise cgMLSTError(f'Error parsing {assembly.name}: {e}') from e

    return Assembly(path = assembly, name = basename, contigs = input_seq)

def add_new_allele(best_match: str, gene_file: pathlib.Path, input_file_obj: Assembly) -> str:
    '''
    Add a new allele to the gene file based on the best hit and return the new allele number.
    '''
    last_id = None
    for record in SeqIO.parse(gene_file, 'fasta'):
        last_id = int(record.id.strip().split('_')[-1]) + 1

    sequence = str(input_file_obj.fna[best_match])
    sequence_for_write = '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))

    with open(gene_file, 'a') as f:
        f.write(f'\n>{last_id}\n')
        f.write(sequence_for_write + '\n')

    return str(last_id)

def pending_result(gene_file: pathlib.Path, input_file_obj: Assembly) -> str:
    '''
    Determine the allele number for a gene based on hits.
    '''
    best_match = ''
    best_coverage = 0.0
    best_identity = 0.0
    best_ST = ''
    for i in input_file_obj.blast_hits:
        coverage = i.query_cov
        identity = i.pident
        if coverage >= best_coverage and identity >= best_identity:
            best_coverage = coverage
            best_identity = identity
            best_match = str(i.sseqid)
            best_ST = str(i.qseqid)
    if int(best_identity) < 100:
        best_ST = ''
    if best_ST:
        return best_ST
    else:
        if best_match:
            best_ST = add_new_allele(best_match, gene_file, input_file_obj)
        else:
            best_ST = '0'
        return best_ST

def generate_output(output: pathlib.Path, cgMLST_list: list[str]):
    '''
    Generate the output file with header if it does not exist.
    '''
    if pathlib.Path(output).is_file():
        return
    headers = ['Isolate']
    headers += cgMLST_list
    with open(output, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(output: pathlib.Path, input_name: str, cgMLST_list: list[str], cgMLST_type: dict[str, str]):
    '''
    Write the cgMLST type for the input genome to the output file.
    '''
    line = [input_name]
    for i in cgMLST_list:
        line.append(cgMLST_type[i])
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')  

def process_gene(gene: str, cgMLST_path: pathlib.Path, input_file_obj: Assembly) -> tuple[str, str]:
    '''
    Process a single gene by mapping it to the input assembly and determining the allele number.
    '''
    gene_fasta = cgMLST_path / f'{gene}.fasta'
    input_file_obj.map(gene_fasta)
    best_ST = pending_result(gene_fasta, input_file_obj)
    return gene, best_ST

def main(argv: list[str] | None = None):
    argv = sys.argv[1:] if argv is None else argv
    args = get_argument().parse_args(argv)
    check_programs_shutil(['makeblastdb', 'blastn', 'blastx', 'prodigal'])
    if args.mode == 'core_genome':
        workpath = pathlib.Path.cwd()
        # Separate the reference core genome to every single file
        cgMLST_path, cgMLST_list = separate_core(args.reference, workpath)
    else:
        cgMLST_path = pathlib.Path(args.reference)
        cgMLST_list = [file.stem for file in pathlib.Path(args.reference).iterdir()]
    # Run this pipeline for each single input genome
    for inputfile in tqdm(args.input, desc = 'Processing files', unit='genome'):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = gunzip_assembly(inputfile, temp_dir)
            input_file_obj = parse_input_assembly(input_file)
            if input_file_obj.name == 'core_genome':   # skip the reference core genome file itself.
                continue
            else:
                input_file_obj.prodigal()
                input_file_obj.make_input_blastdb(temp_dir)

                allele_files = list(cgMLST_path.glob('*.fasta'))
                gene_names = [f.stem for f in allele_files]
                cgMLST_type = dict()
                # Process each gene in parallel using multiprocessing Pool
                with Pool(processes = args.threads) as pool:
                    results = pool.starmap(
                        process_gene,
                        [(gene, cgMLST_path, input_file_obj) for gene in gene_names]
                    )

                cgMLST_type = dict(results)

                generate_output(args.output, cgMLST_list)
                output(args.output, input_file_obj.name, cgMLST_list, cgMLST_type)
   
if __name__ == '__main__':
    raise SystemExit(main())

    
