# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 08:46:19 2023

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

import argparse
import sys
import pathlib
import multiprocessing
import subprocess
import time
import re
from Bio import SeqIO

__version__ = '1.0'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'Bac_cgMLST', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-r', '--reference', required = True, type = str, 
                                help = 'reference core genome file')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'cgMLST_output.txt',
                              help = 'Output file')

    # Parameters
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'cgMLST v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['makeblastdb', 'blastn', 'blastx', 'prodigal']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Error: could not find %s tool' % i, file=sys.stderr)
            sys.exit(1)

def separate_core(reference_cgMLST, workpath):
    if not pathlib.Path('cgMLST').is_dir():
        pathlib.Path('cgMLST').mkdir()
    cgpath = workpath / 'cgMLST'
    cglist = []
    reference = parse_inputfile(reference_cgMLST)
    for key, value in reference.items():
        gene_name = ''
        if str(key).strip().split('_')[-1].isdigit():
            gene_name = str(key).strip().split('_')[-2] + '_' + str(key).strip().split('_')[-1]
        else:
            gene_name = str(key).strip().split('_')[-1]
        cglist.append(gene_name)
        filepath = cgpath / gene_name
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
    
def parse_inputfile(inputfile):
    # Generate a dict, with contig_name as key and contig_sequence as value
    input_seq = {}
    for contig in SeqIO.parse(inputfile, 'fasta'):
        input_seq[contig.name] = contig.seq
    if not input_seq:
        print('invalid FASTA file: %s', inputfile)
        sys.exit(1)
    return input_seq
        
def prodigal(inputfile, input_seq):
    # Generate a list with the informations and sequence of each orf in the inputfile
    command = 'prodigal -f sco -i %s -m' % inputfile
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
    raw_orfs = result.stdout.decode().replace('\r', '')
    orfs = group_prodigal(raw_orfs)
    for orf in orfs:
        orf.sequence = input_seq[orf.contig][orf.start-1:orf.end]
    return orfs

def group_prodigal(raw_orfs):
    # group every orf to a list
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
    return orfs

class Orf:
    # Parse the prodigal results
    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.sequence = str()
        
def make_input_blastdb(input_orfs):
    file_name = 'input_faa.faa'
    input_faa_dict = dict()
    input_faa = open(file_name, 'wt')
    for orf in input_orfs:
        gene_code = str(orf.contig) + '_' + str(orf.start) + '_' + str(orf.end) + '_' + str(orf.strand)
        sequence = str(orf.sequence)
        input_faa_dict[gene_code] = sequence
        sequence_for_write = ''
        while len(sequence) > 60:
            sequence_for_write += sequence[:60] + '\n'
            sequence = sequence[60:]
        if sequence:
            sequence_for_write += sequence
        input_faa.write('>')
        input_faa.write(gene_code)
        input_faa.write('\n')
        input_faa.write(sequence_for_write)
        input_faa.write('\n')
    input_faa.close()
    makeblastdb(file_name)
    pathlib.Path(file_name).unlink()
    return file_name, input_faa_dict
    
def makeblastdb(faa):
    makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', faa]
    makeblastdb_process = subprocess.run(makeblastdb_command, stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)

def pending_result(blast_hits, input_faa_dict, file):
    # Find the best serotype of inputfile by sequence alignment
    best_match = ''
    best_coverage = 0.0
    best_identity = 0.0
    best_ST = ''
    for i in blast_hits:
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
            best_ST = add_new_allele(best_match, input_faa_dict, file)
        else:
            best_ST = '0'
        return best_ST

def add_new_allele(best_match, input_faa_dict, file):
    single_gene = open(file, 'rt')
    count = ''
    for line in single_gene:
        if line.startswith('>'):
            count = line.strip()[1:]
    single_gene.close()
    number = eval(count) + 1
    sequence = str(input_faa_dict[best_match])
    sequence_for_write = ''
    while len(sequence) > 60:
        sequence_for_write += sequence[:60] + '\n'
        sequence = sequence[60:]
    if sequence:
        sequence_for_write += sequence
    single_gene = open(file, 'at')
    single_gene.write('>')
    single_gene.write(str(number))
    single_gene.write('\n')
    single_gene.write(sequence_for_write)
    single_gene.write('\n')
    return str(number)
    
def run_blast(query, subject, threads):
    # Do blast, iterator the result to a list
    blast_hits = []
    command = ['blastn', '-query', query, '-db', subject, '-num_threads', str(threads), '-outfmt', 
               '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    for line in line_iterator(out):
        blast_hits.append(BlastResult(line))
    return blast_hits

def line_iterator(line_breaks):
    # Handle the BLAST output and remove the line breaks 
    line = -1
    while True:
        nextline = line_breaks.find('\n', line + 1)
        if nextline < 0:
            break
        yield line_breaks[line + 1:nextline]
        line = nextline

class BlastResult(object):
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

def strip_suffix(inputfile):
    # Strip the suffix of inputfile
    if inputfile.lower().endswith('.fa'):
        inputfile = inputfile[:-3]
    elif inputfile.lower().endswith('.fna'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fas'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fasta'):
        inputfile = inputfile[:-6]
    return inputfile

def generate_output(output, cgMLST_list):
    # Generate a blank output table file
    if pathlib.Path(output).is_file():
        return
    headers = ['Isolate']
    headers += cgMLST_list
    with open(output, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(output, input_name, cgMLST_list, cgMLST_type):
    # Generate output
    line = [input_name]
    for i in cgMLST_list:
        line.append(cgMLST_type[i])
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')  
        
def main():
    print('If you have any questions or suggestions for Bac_cgMLST, please contact Genglin Guo, e-mail: 2019207025@njau.edu.cn')
    starttime = time.perf_counter()
    # Initialize
    args = get_argument().parse_args()
    check_dependencies()
    workpath = pathlib.Path.cwd()
    # Separate the reference core genome to every single file
    cgMLST_path, cgMLST_list = separate_core(args.reference, workpath)
    # Run this pipeline for each single input genome
    for inputfile in args.input:
        if inputfile == 'core_genome.fasta':
            continue
        else:
            input_name = strip_suffix(inputfile)
            print('start to typing {}'.format(input_name))
            input_seq = parse_inputfile(inputfile)
            input_orfs = prodigal(inputfile, input_seq)
            inputdb, input_faa_dict = make_input_blastdb(input_orfs)
            cgMLST_type = dict()
            for file in pathlib.Path(cgMLST_path).iterdir():
                blast_hits = run_blast(str(file), inputdb, args.threads)
                best_ST = pending_result(blast_hits, input_faa_dict, file)
                gene_name = str(file).strip().split('/')[-1]
                cgMLST_type[gene_name] = best_ST
            for file in pathlib.Path(workpath).iterdir():
                if str(file).endswith(('ndb', 'nhr', 'nin', 'not', 'nsq', 'ntf', 'nto')):
                    pathlib.Path(file).unlink()
            generate_output(args.output, cgMLST_list)
            output(args.output, input_name, cgMLST_list, cgMLST_type)
            print('Typing for {} is complete'.format(input_name))
    endtime = time.perf_counter() - starttime
    per_genome_time = endtime / len(args.input)
    print('{:.1f}h{:.1f}m{:.1f}s for one genome'.format(per_genome_time // 3600, per_genome_time % 3600 // 60, per_genome_time % 60))
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
main()


