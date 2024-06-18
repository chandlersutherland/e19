#!/bin/bash
module load minimap2
module load samtools/1.14 

#define variables 
genome=/global/scratch/users/chandlersutherland/e19/re_digest/input_genome/Athaliana.fa
gene_bed=/global/scratch/users/chandlersutherland/Athaliana/Atha_genes.bed
nlr_bed=/global/home/users/chandlersutherland/e14/data/all_NLR.bed
base_dir=/global/scratch/users/chandlersutherland/e19/re_digest

#loop through short list 
for enzyme in 'AluI' 'BfaI' 'CviAII' 'CviQI' 'DpnII' 'FatI' 'HaeIII' 'HpyCH4V' 'MluCl' 'MseI' 'MspI' 'TaqI'
do 
    echo $enzyme
    cd ${base_dir}/${enzyme}

    cat *_seq_1.fasta > ${enzyme}_seq1.fasta 
    cat *_seq_2.fasta > ${enzyme}_seq2.fasta

    #map and sort 
    minimap2 -a $genome ${enzyme}_seq1.fasta ${enzyme}_seq2.fasta |
    samtools view -h |
    samtools sort > ${enzyme}_seq.sorted.bam 

    #index for IGV 
    samtools index ${enzyme}_seq.sorted.bam

    #calculate coverage of all genes and of NLRs
    samtools bedcov $gene_bed ${enzyme}_seq.sorted.bam > gene_seq_coverage.tsv
    samtools bedcov $nlr_bed ${enzyme}_seq.sorted.bam > nlr_seq_coverage.tsv 
    #all genome coverage
    samtools coverage ${enzyme}_seq.sorted.bam > ${enzyme}_seq_coverage.tsv
done 