#!/bin/bash -l
#
#$ -o blast+-2.2.29.out
#$ -e blast+2.2.29.err
#$ -N blast_test
#$ -cwd

source /etc/profile.d/modules.sh
module load ncbi-blast+/2.2.29 

blastn -query test.fasta  -db /blast/db/pdbnt -num_descriptions 10 -num_alignments 20 -out test_vs_nt.blastn.29

module unload ncbi-blast+/2.2.29 
