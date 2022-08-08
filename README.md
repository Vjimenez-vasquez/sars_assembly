# sars_assembly
a collection of commands for sars-cov-2 genome assembly 

#0# juntar los fastq files generados por nextseq#
mkdir fastq ;
mv */*fastq.gz fastq ;
cd fastq ;
ls -lh ;
ls *fastq.gz | sed "s/_L.*//g" | sort | uniq > file.txt ; 
lista=( `ls *gz | sed "s/_L.*//g" | sort -d | uniq`) ;
for i in ${lista[@]};
do cat ${i}_L001_R2_001.fastq.gz ${i}_L002_R2_001.fastq.gz ${i}_L003_R2_001.fastq.gz ${i}_L004_R2_001.fastq.gz > ../${i}_R2_001.fastq.gz; 
cat ${i}_L001_R1_001.fastq.gz ${i}_L002_R1_001.fastq.gz ${i}_L003_R1_001.fastq.gz ${i}_L004_R1_001.fastq.gz > ../${i}_R1_001.fastq.gz ;
done & pwd ;
cd .. ;
ls -lh ;

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz 
do
prefix=$(basename $r1 _R1_001.fastq.gz)
r2=${prefix}_R2_001.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 2 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 2 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 25, frecuencia de nucleotido predominante: 55%#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 20 ;
done ;

#6# remover subproductos#
rm *.fastq.gz.fa *.qual.txt ; 
cat *.fa > secuencias.fasta ;

#7# estimate mean depth of coverage#
for r1 in *bam
do
prefix=$(basename $r1 .bam)
DEP=( `samtools depth $r1 | awk '{sum+=$3}END{print sum/29903}' `)
NPE=( `seqtk comp ${prefix}.fa | awk '{x=$3+$4+$5+$6;y=29903;print 1-(y-x)/y}'`)
echo "${prefix} ${DEP}x $NPE" >> echo.txt ;
done ;

