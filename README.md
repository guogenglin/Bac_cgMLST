# Bac_cgMLST
A tool for generate cgMLST from pan-genome analysis data

# External Dependencies
BLAST+

Prodigal

# Guide for users
This tool could generated a personalized cgMLST result for your sequence data, at first, you should find the core-genome of your data, the author was used Roary to do the Pan-genome analyse normally.

After you finished the analyse, the core genome and the accessory genome of your sequence could be separated based on the "gene_presence_and_absence.csv"

First:

Copy the core genome name from first coloum from "gene_presence_and_absence.csv" to a text file (Strongly suggest you only select the gene presenced in every genome sequence) and renamed it as "core_name.txt", put the "pan_genome_reference.fa" file generated by Roary in the same filefolder, run 'core_accessory_separate.py'
Four files will be generated, 'core_genome.fasta', 'accessory_genome.fasta', 'core_proteome.faa' and 'accessory_proteome.faa', only 'core_genome.fasta' will be used in cgMLST, others could be used in your further analysis of Pan-genome.

Second:

Run Bac_cgMLST.py

```
Bac_cgMLST.py [-i] [-r] [-o] [-t] [-v]
Input and Output:
  -i, --input             Input FASTA file
  -r, --reference         Reference core genome file
  -o, --output            Output file
Parameters:
  -t, --threads           Threads to use for BLAST searches
  -v, --version           Show version number and exit
```

There will be a folder, which contain the sequence of every alleles of every genes, and a summary table named 'cgMLST_output.txt', you can use this table to generate a minimum spanning tree using GrapeTree or other tools, and combining analyse with some extra data, such as virulence, serotype, isolation source, etc. Here is an example output:

![GPMS](https://github.com/guogenglin/Bac_cgMLST/assets/108860907/c6f70bcf-097e-4b5b-81d6-84d28382a6a4)

# Quick usage
``` Python
python Bac_cgMLST.py -i *.fasta -r core_genome.fasta
```

