#!/bin/bash

#0# juntar los fastq files generados por nextseq#
mkdir fastq_original ;
mv *.gz fastq_original/ ; 
cd fastq_original/ ; 
ls *fastq.gz | sed "s/_L.*//g" | sort | uniq > file.txt ; 
lista=( `ls *gz | sed "s/_L.*//g" | sort -d | uniq`) ;
for i in ${lista[@]} ;
do cat ${i}_L001_R2_001.fastq.gz ${i}_L002_R2_001.fastq.gz ${i}_L003_R2_001.fastq.gz ${i}_L004_R2_001.fastq.gz > ../${i}_L001_R2_001.fastq.gz ; 
cat ${i}_L001_R1_001.fastq.gz ${i}_L002_R1_001.fastq.gz ${i}_L003_R1_001.fastq.gz ${i}_L004_R1_001.fastq.gz > ../${i}_L001_R1_001.fastq.gz ;
done ;
cd .. ;
exit
