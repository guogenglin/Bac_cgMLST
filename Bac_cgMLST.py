# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 08:46:19 2023

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

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
from typing import Generator
from subprocess import Popen, PIPE
from Bio.SeqIO.FastaIO import SimpleFastaParser


__version__ = '1.1'

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
        gene_name = ''
        if str(key).strip().split('_')[-1].isdigit():
            gene_name = str(key).strip().split('_')[-2] + '_' + str(key).strip().split('_')[-1]
        else:
            gene_name = str(key).strip().split('_')[-1]
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

class Contig:
    """
    Represents a single contig in an assembly, with a name, description, and sequence.
    """
    def __init__(self, name: str, description: str, seq: Seq):
        self.name = name
        self.description = description
        self.seq = seq

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:
        return f"<Contig {self.name}, length={len(self)}>"

class Alignment:
    '''
    Represents a single alignment between a query and a contig, with coordinates, strand, and optional tags.
    '''
    def __init__(
            self, q: str | None = None, q_len: int | None = 0, q_st: int | None = 0,
            q_en: int | None = 0, strand: str | None = None, ctg: str | None = None,
            ctg_len: int | None = 0, r_st: int | None = 0, r_en: int | None = 0,
            mlen: int | None = 0, blen: int | None = 0, mapq: int | None = 0,
            tags: dict | None = None):
        self.q = q or 'unknown'  # Query sequence name
        self.q_len = q_len or 0  # Query sequence length, None -> 0
        self.q_st = q_st or 0  # Query start coordinate (0-based)
        self.q_en = q_en or 0  # Query end coordinate (0-based)
        self.strand = strand or 'unknown'  # ‘+’ if query/target on the same strand; ‘-’ if opposite
        self.ctg = ctg or 'unknown'  # Target sequence name
        self.ctg_len = ctg_len or 0  # Target sequence length, None -> 0
        self.r_st = r_st or 0  # Target start coordinate on the original strand (0-based)
        self.r_en = r_en or 0  # Target end coordinate on the original strand (0-based)
        self.mlen = mlen or 0  # Number of matching bases in the alignment
        self.blen = blen or 0  # Number bases, including gaps, in the alignment
        self.mapq = mapq or 0  # Mapping quality (0-255 with 255 for missing)
        self.tags = tags or {}  # {tag: value} pairs

    @classmethod
    def from_paf_line(cls, line: str):
        '''
        Parse a PAF line into an Alignment object. Expects at least 12 columns, with optional tags after.
        '''
        fields = line.rstrip().split('\t')

        if len(fields) < 12:
            raise cgMLSTError(f'Invalid PAF line (<12 columns): {line}')
        
        try:
            def parse_tag(tag: str):
                key, type_, value = tag.split(':', 2)
                if type_ == 'i':
                    return key, int(value)
                if type_ == 'f':
                    return key, float(value)
                return key, value
            
            tags = dict(parse_tag(t) for t in fields[12:])
            
            return Alignment(q = fields[0], q_len = int(fields[1]), q_st = int(fields[2]), q_en = int(fields[3]), 
                             strand = fields[4], ctg = fields[5], ctg_len = int(fields[6]), r_st = int(fields[7]), 
                             r_en = int(fields[8]), mlen = int(fields[9]), blen = int(fields[10]), 
                             mapq = int(fields[11]), tags = tags)
        
        except Exception as e:
            raise cgMLSTError(f'Error parsing PAF line: {line}') from e

    @property
    def identity(self) -> float:
        ''' Returns percent identity (0-100). Requires --eqx in minimap2. '''
        return 100.0 * self.mlen / self.blen if self.blen > 0 else 0.0

    @property
    def coverage(self) -> float:
        ''' Returns query coverage percentage (0-100). '''
        return 100.0 * (self.q_en - self.q_st) / self.q_len if self.q_len > 0 else 0.0
    
    def is_exact(self) -> bool:
        ''' Returns True if the alignment is a 100% match across the full query length. '''
        return self.mlen == self.blen and self.q_en - self.q_st == self.q_len

    @property
    def score(self) -> int:
        ''' Returns the alignment score (AS tag) if available. '''
        return self.tags.get('AS', 0)
    
    def __repr__(self) -> str:
        return f'{self.q}:{self.q_st}-{self.q_en} {self.ctg}:{self.r_st}-{self.r_en} {self.strand}'

    def __len__(self) -> int:
        return self.q_en - self.q_st

    def __getattr__(self, item):
        # First look in tags if not a normal attribute
        if item in self.tags:
            return self.tags[item]
        raise AttributeError(f'{self.__class__.__name__} object has no attribute {item}')
    
class Assembly:
    '''
    Represents an assembly with multiple contigs. Provides methods to build minimap2 index and map queries.
    '''
    def __init__(self, path: pathlib.Path, name: str | None = None, contigs: dict[str, Contig] | None = None):
        self.path = pathlib.Path(path).resolve()
        self.name = name or self.path.stem
        self.contigs = contigs or {}
        self.mmi_index: pathlib.Path | None = None

    def __repr__(self) -> str:
        return f"<Assembly {self.name}, {len(self.contigs)} contigs>"

    def __len__(self) -> int:
        return sum(len(c) for c in self.contigs.values())

    def seq(self, ctg: str, start: int, end: int, strand: str = '+') -> Seq:
        '''
        Retrieve the sequence for a given contig and coordinates, accounting for strand.
        '''
        if ctg not in self.contigs:
            raise cgMLSTError(f"Contig {ctg} not found in assembly {self.name}")
        
        seq = self.contigs[ctg].seq[start:end]
        return seq if strand == '+' else seq.reverse_complement()

    def build_minimap2_index(self, temp_dir: str) -> pathlib.Path:
        '''
        Build a minimap2 index for the assembly and return the path to the index file.
        '''
        if self.mmi_index and self.mmi_index.exists():
            return self.mmi_index

        temp_path = pathlib.Path(temp_dir)
        self.mmi_index = (temp_path / (uuid.uuid4().hex + '.mmi')).resolve()

        cmd = ['minimap2', '-d', str(self.mmi_index), str(self.path)]


        p = subprocess.run(cmd, capture_output = True, text = True)

        if p.returncode != 0:
            raise cgMLSTError(f'minimap2 failed indexing {self.name}:\n{p.stderr}')

        return self.mmi_index
    
    def map(self, query: pathlib.Path | str | None, threads: int = 1, extra_args: str = '') -> Generator[Alignment, None, None]:
        '''
        Map a query sequence to the assembly using minimap2 and yield Alignment objects.
        '''
        if not self.mmi_index:
            raise cgMLSTError(f'Index not built for assembly {self.name}')
        
        cmd = f"minimap2 -c {extra_args} -t {threads} '{self.mmi_index}' '{query}'"


        proc = Popen(cmd, shell = True, stdin = PIPE, stdout = PIPE, stderr = PIPE, text = True)
        
        stdout, stderr = proc.communicate()

        for line in stdout.splitlines():
            yield Alignment.from_paf_line(line)
    
def parse_input_assembly(assembly: pathlib.Path) -> Assembly:
    '''
    Parse the input assembly FASTA file into an Assembly object.
    '''
    basename = assembly.stem

    # Parse fasta file
    contigs = {}

    try:
        with open(assembly, mode = 'rt') as f:
            for header, seq in SimpleFastaParser(f):
                parts = header.split(maxsplit = 1)
                name = parts[0]
                description = parts[1] if len(parts) > 1 else ''
                contigs[name] = Contig(name, description, Seq(seq))

    except Exception as e:
        raise cgMLSTError(f'Error parsing {assembly.name}: {e}') from e

    return Assembly(path = assembly, name = basename, contigs = contigs)

def add_new_allele(hit: Alignment, gene_file: pathlib.Path, input_file_obj: Assembly) -> str:
    '''
    Add a new allele to the gene file based on the best hit and return the new allele number.
    '''
    last_id = None
    for record in SeqIO.parse(gene_file, 'fasta'):
        last_id = int(record.id.strip().split('_')[-1]) + 1

    sequence = str(input_file_obj.seq(hit.ctg, hit.r_st, hit.r_en, hit.strand))
    sequence_for_write = '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))

    with open(gene_file, 'a') as f:
        f.write(f'\n>{last_id}\n')
        f.write(sequence_for_write + '\n')

    return str(last_id)

def pending_result(hits: list[Alignment], gene_file: pathlib.Path, input_file_obj: Assembly) -> str:
    '''
    Determine the allele number for a gene based on hits.
    '''

    if not hits:
        return '0'  # No hits, assign allele number 0
    
    max_id = max(h.identity for h in hits)
    top_hits = [h for h in hits if h.identity == max_id]   # Multiple hits could have the same max identity
    max_score = max(h.score for h in top_hits)   # Filter further by score
    filtered_hits = [h for h in top_hits if h.score == max_score]

    for hit in filtered_hits:
        if hit.is_exact():
            return hit.q.split('_')[-1]  # Exact match, return allele number
    return add_new_allele(filtered_hits[0], gene_file, input_file_obj)  # No exact match, add new allele based on best hit

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
        
def main(argv: list[str] | None = None):
    argv = sys.argv[1:] if argv is None else argv
    args = get_argument().parse_args(argv)
    check_programs_shutil(['minimap2'])
    if args.mode == 'core_genome':
        workpath = pathlib.Path.cwd()
        # Separate the reference core genome to every single file
        cgMLST_path, cgMLST_list = separate_core(args.reference, workpath)
    else:
        cgMLST_path = pathlib.Path(args.reference)
        cgMLST_list = [file.stem for file in pathlib.Path(args.reference).iterdir()]
    # Run this pipeline for each single input genome
    for inputfile in tqdm(args.input, desc = 'Processing files'):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = gunzip_assembly(inputfile, temp_dir)
            input_file_obj = parse_input_assembly(input_file)
            if input_file_obj.name == 'core_genome':   # skip the reference core genome file itself.
                continue
            else:
                input_file_obj.build_minimap2_index(temp_dir)

                allele_files = list(cgMLST_path.glob('*.fasta'))
                gene_names = [f.stem for f in allele_files]
                cgMLST_type = dict()

                for gene in gene_names:
                    total_hits = []
                    gene_fasta = cgMLST_path / f'{gene}.fasta'
                    for aln in input_file_obj.map(gene_fasta, extra_args = '--end-bonus=10 --eqx -x asm5'):
                        # Filters
                        if aln.identity >= 90.0 and aln.coverage >= 80.0:
                            total_hits.append(aln)
                    best_ST = pending_result(total_hits, gene_fasta, input_file_obj)
                    cgMLST_type[gene] = best_ST

                generate_output(args.output, cgMLST_list)
                output(args.output, input_file_obj.name, cgMLST_list, cgMLST_type)
   
if __name__ == '__main__':
    raise SystemExit(main())

    