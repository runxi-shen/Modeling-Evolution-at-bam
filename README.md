# Modeling-Evolution-at-bam

This repository stores all the analysis scripts used in our publication:

> Runxi Shen*, Miwa Wenzel*, Philipp W. Messer, Charles F. Aquadro (2022) A model of functionally buffered deleterious mutations can lead to signatures of positive selection distinguishable from an evolutionary conflict model.

__Description of analysis scripts:__

python_wrt_\*.py: scripts used to run the SLiM simulations by writing SLiM scripts first and then running the script in SLiM program;

vcf2fasta.py: convert vcf files output from SLiM to fasta files for downstream iMKT analysis;

sfsFromFasta_2fold.py: output the input files for iMKT;

iMKT.sh: bash script for running iMKT analysis
