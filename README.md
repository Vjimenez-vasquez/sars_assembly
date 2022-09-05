## sars_assembly (INS)
A collection of commands for sars-cov-2 genome assembly derived from ILLUMINA NEXTSEQ platform.

## STEP -1: explanations
```r
merge_fastq.sh : merge 4 forward fastq files and 4 forward fastq files (only for NextSeq)
only_assembly.sh : assembly pipeline
tabla_final.R : prepare coverage, completeness and lineage identification table
```

## STEP 0: just in case generate a naked script.sh file
```r
cd $HOME && touch script.sh && chmod +x script.sh ;
cd $HOME && echo '#!/bin/bash' > script.sh && echo '# -*- ENCODING: UTF-8 -*-' >> script.sh ;
mv script.sh paste/your/working/directory/ ;

or 

chmod +x script.sh 
```

## STEP 1: merge_fastq.sh
```r
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
```

## STEP 2: only_assembly.sh
```r
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

#10# alinear con referencia#
mafft --6merpair --thread 15 --addfragments secuencias.fasta reference.fasta > alineado.fasta ; 

#11# header solo NETLAB#
sed 's/Consensus_//g' alineado.fasta > uno.fasta ; 
sed 's/_S.*//g' uno.fasta > alineado_clean.fasta ; 
rm uno.fasta ; 
aliview alineado_clean.fasta ;
ls -lh ;
```

## STEP 3: tabla_final.R
```r
run <- 1055
a <- read.csv("lineage_report.csv", header=TRUE)
b <- gsub("Consensus_","",a$taxon)
a$strain <- gsub("_threshold_0.6_quality_25","",b)
d <- read.csv("profundidad_ns.txt", header=FALSE, sep=" ")
names(d) <- c("strain","Mean_Coverage","ch")
d$ch <- d$ch*100
e <- merge(a,d,by="strain", all.x=TRUE)
f <- e[,c(1,3,6,18,19)]
names(f) <- c("taxon","lineage","scorpio_call","MeanCoverage","NonN")
g <- data.frame(f$taxon,corrida=rep(run,length(f$taxon)),f$lineage,f$scorpio_call,f$MeanCoverage,f$NonN)
names(g) <- c("taxon","corrida","lineage","scorpio_call","MeanCoverage","NonN")
write.csv(g,"nexstrain2gisaid_vero.csv",row.names=FALSE)
g$taxon <- gsub("_.*","",g$taxon)
write.csv(g,"nexstrain2gisaid.csv", row.names=FALSE)
q("no") ;
```
