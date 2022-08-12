#!/bin/bash

#we need trimmomatic, fastp, bwa, samtools, special java script (will be installed automatically)
#if not here is a link http://lindenb.github.io/jvarkit/Biostar84452.html

#we give 4 arguments: r1, r2 for PE; file with adapter seqs; directory with ref files (there should be only ref files in fasta format)
#script creates dir Results_fragseq with bam files and statistics file

R1=$1
R2=$2
adapter_file=$3

#NOW WE WORK WITH ONE REF FILE 
#several references from the pointed dir will be combined into one
ref_dir=$4
awk 1 ${ref_dir}/*.fasta > ${ref_dir}/combi_ref.fasta
f=${ref_dir}/combi_ref.fasta
#f=$4


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


#piece of the code responsible for installation of the script that removes soft clipped bases from the sam (sequences)
CDIR=$(pwd)

if ! [ -f jvarkit/dist/biostar84452.jar ]; then
echo 'No file biostar84452.jar'
git clone "https://github.com/lindenb/jvarkit.git"
cd jvarkit
./gradlew biostar84452
cd $CDIR
fi



bwa mem -t 4 -k 20 ${f} Results_fragseq/${noext_r1}_paired_trimq_a.fastq Results_fragseq/${noext_r2}_paired_trimq_a.fastq | samtools view -S -b | java -jar jvarkit/dist/biostar84452.jar | samtools view -S -b | samtools view --threads 4 -b -f 2 > Results_fragseq/${noext}_${ref_pre}_mapped.bam
samtools view -bq 1 Results_fragseq/${noext}_${ref_pre}_mapped.bam > Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam








