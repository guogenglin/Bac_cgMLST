# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:56:01 2023

@author: Genglin Guo
@e-mail: 2019207025.njau.edu.cn
"""

import pathlib
from Bio import SeqIO
from Bio.Seq import Seq

core_genome = open('core_genome.fasta', 'wt')
accessory_genome = open('accessory_genome.fasta', 'wt')
reference_genes = open('pan_genome_reference.fa', 'rt')
core_name = open('core_name.txt', 'rt')
reference_genes_revise = open('pan_genome_reference.fasta', 'wt')
for line in reference_genes:
    if line.startswith('>'):
        new_line = line.replace(' ', '+')
        reference_genes_revise.write(new_line)
    else:
        reference_genes_revise.write(line)
reference_genes.close()
reference_genes_revise.close()
reference_genes_revise = open('pan_genome_reference.fasta', 'rt')
all_genes = dict()
for gene in SeqIO.parse(reference_genes_revise, 'fasta'):
    all_genes[gene.name] = gene.seq
reference_genes_revise.close()
pathlib.Path('pan_genome_reference.fasta').unlink()
core_genes = list()
for line in core_name:
    core_genes.append(line.strip())
core_name.close()
core_gene_num = 0
accessory_gene_num = 0
for name, seq in all_genes.items():
    separate_name = name.split('+')
    name = name.replace('+', '_')
    sequence_for_write = ''
    sequence = str(seq)
    while len(sequence) > 60:
        sequence_for_write += sequence[:60] + '\n'
        sequence = sequence[60:]
    if sequence:
        sequence_for_write += sequence
        sequence_for_write += '\n'
    if separate_name[-1] in core_genes:
        core_gene_num += 1
        core_genome.write('>')
        core_genome.write(name)
        core_genome.write('\n')
        core_genome.write(sequence_for_write)
        core_genome.write('\n')
    else:
        accessory_gene_num += 1
        accessory_genome.write('>')
        accessory_genome.write(name)
        accessory_genome.write('\n')
        accessory_genome.write(sequence_for_write)
        accessory_genome.write('\n')
core_genome.close()
accessory_genome.close()
print('{} core genes and {} accessory genes were identified'.format(str(core_gene_num), str(accessory_gene_num))) 
core_proteome = open('core_proteome.faa', 'wt')
accessory_proteome = open('accessory_proteome.faa', 'wt')
core_genes = dict()
accessory_genes = dict()
for gene in SeqIO.parse('core_genome.fasta', 'fasta'):
    core_genes[gene.name] = gene.seq
for gene in SeqIO.parse('accessory_genome.fasta', 'fasta'):
    accessory_genes[gene.name] = gene.seq
core_prot_num = 0
accessory_prot_num = 0
for name, seq in core_genes.items():
    core_prot_num += 1
    core_proteome.write('>')
    core_proteome.write(name)
    core_proteome.write('\n')
    aa_seq = str(seq.translate())
    aa_seq_for_write = ''
    while len(aa_seq) > 60:
        aa_seq_for_write += aa_seq[:60] + '\n'
        aa_seq = aa_seq[60:]
    if aa_seq:
        aa_seq_for_write += aa_seq
        aa_seq_for_write += '\n'
    core_proteome.write(aa_seq_for_write)
    core_proteome.write('\n')
for name, seq in accessory_genes.items():
    accessory_prot_num += 1
    accessory_proteome.write('>')
    accessory_proteome.write(name)
    accessory_proteome.write('\n')
    aa_seq = str(seq.translate())
    aa_seq_for_write = ''
    while len(aa_seq) > 60:
        aa_seq_for_write += aa_seq[:60] + '\n'
        aa_seq = aa_seq[60:]
    if aa_seq:
        aa_seq_for_write += aa_seq
        aa_seq_for_write += '\n'
    accessory_proteome.write(aa_seq_for_write)
    accessory_proteome.write('\n')
core_proteome.close()
accessory_proteome.close()
print('{} core genes and {} accessory genes were translated to proteins'.format(str(core_prot_num), str(accessory_prot_num))) 
