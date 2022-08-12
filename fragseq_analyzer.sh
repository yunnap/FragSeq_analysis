#!/bin/bash

#we need trimmomatic, fastp, bwa, samtools, special java script
#regarding the last one link here ...
#we give 4 arguments: r1, r2 for PE; file with adapter seqs; directory where our ref files are (there should be only ref files)
#script creates dir Results_fragseq

R1=$1
R2=$2
adapter_file=$3

#NOW WE WORK WITH ONE REF FILE 
#ref_dir=$4
f=$4
basename_ref=${f##*/}
ref_pre=${basename_ref%%.*}

basename_r1=${R1##*/}
noext_r1=${basename_r1%%.*}
basename_r2=${R2##*/}
noext_r2=${basename_r2%%.*}

path=$1
basename=${path##*/}
noext=${basename%%.*}_R2

mkdir Results_fragseq


trimmomatic PE $R1 $R2 out_fw_paired.fq.gz out_fw_unpaired.fq.gz out_rv_paired.fq.gz out_rv_unpaired.fq.gz CROP:100
rm out_fw_unpaired.fq.gz
rm out_rv_unpaired.fq.gz

fastp --in1 out_fw_paired.fq.gz --in2 out_rv_paired.fq.gz --out1 f_paired_trimq.fastq --out2 r_paired_trimq.fastq --qualified_quality_phred 20 --unqualified_percent_limit 10 -A -L 
rm out_fw_paired.fq.gz
rm out_rv_paired.fq.gz
#| xargs -L2 bash -c 'fastp --in1 ${0} --in2 ${1} --out1 f_paired_trimq_a.fastq --out2 r_paired_trimq_a.fastq -Q -L $adapter_file'
fastp --in1 f_paired_trimq.fastq --in2 r_paired_trimq.fastq --out1 Results_fragseq/${noext_r1}_paired_trimq_a.fastq --out2 Results_fragseq/${noext_r2}_paired_trimq_a.fastq -Q -L $adapter_file
#--adapter_fasta adapters.fasta'
rm f_paired_trimq.fastq
rm r_paired_trimq.fastq

#we need to index a reference genome
bwa index ${f} #kd403-real.fasta

bwa mem -t 4 -k 20 ${f} Results_fragseq/${noext_r1}_paired_trimq_a.fastq Results_fragseq/${noext_r2}_paired_trimq_a.fastq | samtools view -S -b | java -jar jvarkit/dist/biostar84452.jar | samtools view -S -b | samtools view --threads 4 -b -f 2 > Results_fragseq/${noext}_${ref_pre}_mapped.bam
samtools view -bq 1 Results_fragseq/${noext}_${ref_pre}_mapped.bam > Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam


#Rscript --vanilla bam_to_csv.R Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam Results_fragseq/${noext}_${ref_pre}_uni_mapped.csv

#Rscript --vanilla bam_to_csv.R Results_fragseq/Oct14_S10_R1_001_R2_kd403-real_uni_mapped.bam Results_fragseq/Oct14_S10_R1_001_R2_kd403-real_uni_mapped.csv









