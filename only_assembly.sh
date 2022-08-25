#!/bin/bash

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 15 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 25, frecuencia de nucleotido predominante: 60%#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix} -q 25 -t 0.6 -m 20 ;
samtools mpileup -aa -A -d 600000 -B -Q 0 ${prefix}.bam | ivar variants -p ${prefix}.tsv -q 25 -t 0.6 -r reference.fasta ;
done ;

#6# remover subproductos y generar un multifasta#
rm *.fastq.gz.fa *.qual.txt *.fastq.gz.tsv ;
cat *.fa > secuencias.fasta ;

#7# identificar linajes con PANGOLIN#
conda activate pangolin ;
pangolin secuencias.fasta -t 15 --outfile lineage_report.csv --max-ambig 0.05 --min-length 28000 ;
conda deactivate ;
mkdir fasta fastq ;
mv *.gz fastq/ ;
mv *.fasta fasta/ ; 

#8# estimar la profunidad de cobertura y el porcentaje de Ns#
for r1 in *bam
do
prefix=$(basename $r1 .bam)
DEP=( `samtools depth $r1 | awk '{sum+=$3}END{print sum/29903}' `)
NPE=( `seqtk comp ${prefix}.fa | awk '{x=$3+$4+$5+$6;y=29903;print 1-(y-x)/y}'`)
echo "${prefix} ${DEP}x $NPE" >> profundidad_ns.txt ;
done ;

#9# tabla final R#
cp secuencias.fasta lineage_report.csv reference.fasta profundidad_ns.txt tabla_final.R fasta/ ; 
cd fasta/ ;
Rscript tabla_final.R ;
cat reference.fasta secuencias.fasta > input.fasta ;
mafft --6merpair --thread 30 --addfragments secuencias.fasta reference.fasta > alineado.fasta ; 
ls -lh ;
exit
