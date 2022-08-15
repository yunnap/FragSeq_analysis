#!/bin/bash

#we need trimmomatic, fastp, bwa, samtools, special java script (will be installed automatically)
#if not here is a link http://lindenb.github.io/jvarkit/Biostar84452.html

#also we need three additional python scripts for statistics

#we give 4 arguments: r1, r2 for PE; file with adapter seqs; directory with ref files (there should be only ref files in fasta format)
#script creates dir Results_fragseq

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


echo "The report for ${noext}" >> Results_fragseq/report
echo "" >> Results_fragseq/report
echo "-filter reads-" >> Results_fragseq/report
echo "" >> Results_fragseq/report
echo "" >> Results_fragseq/report


trimmomatic PE $R1 $R2 out_fw_paired.fq.gz out_fw_unpaired.fq.gz out_rv_paired.fq.gz out_rv_unpaired.fq.gz CROP:100 >> Results_fragseq/report 2>&1
rm out_fw_unpaired.fq.gz
rm out_rv_unpaired.fq.gz

echo "" >> Results_fragseq/report

fastp --in1 out_fw_paired.fq.gz --in2 out_rv_paired.fq.gz --out1 f_paired_trimq.fastq --out2 r_paired_trimq.fastq --qualified_quality_phred 20 --unqualified_percent_limit 10 -A -L >> Results_fragseq/report 2>&1
rm out_fw_paired.fq.gz
rm out_rv_paired.fq.gz
#| xargs -L2 bash -c 'fastp --in1 ${0} --in2 ${1} --out1 f_paired_trimq_a.fastq --out2 r_paired_trimq_a.fastq -Q -L $adapter_file'
mv fastp.html Results_fragseq/fastp_1.html
mv fastp.json Results_fragseq/fastp_1.json

echo "" >> Results_fragseq/report
fastp --in1 f_paired_trimq.fastq --in2 r_paired_trimq.fastq --out1 Results_fragseq/${noext_r1}_paired_trimq_a.fastq --out2 Results_fragseq/${noext_r2}_paired_trimq_a.fastq -Q -L $adapter_file >> Results_fragseq/report 2>&1
#--adapter_fasta adapters.fasta'
rm f_paired_trimq.fastq
rm r_paired_trimq.fastq
mv fastp.html Results_fragseq/fastp_2.html
mv fastp.json Results_fragseq/fastp_2.json

echo "" >> Results_fragseq/report
echo "" >> Results_fragseq/report
echo "-align reads-" >> Results_fragseq/report
echo "" >> Results_fragseq/report
echo "" >> Results_fragseq/report

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

#############################################################################

bwa mem -t 4 -k 20 ${f} Results_fragseq/${noext_r1}_paired_trimq_a.fastq Results_fragseq/${noext_r2}_paired_trimq_a.fastq | samtools view -S -b | tee Results_fragseq/inter.bam | java -jar jvarkit/dist/biostar84452.jar | samtools view -S -b | samtools view --threads 4 -b -f 2 > Results_fragseq/${noext}_${ref_pre}_mapped_0.bam

samtools view -F 2048 Results_fragseq/${noext}_${ref_pre}_mapped_0.bam -b > Results_fragseq/${noext}_${ref_pre}_mapped.bam 
#samtools view -F 2048 Results_fragseq/${noext}_${ref_pre}_mapped.bam -b > Results_fragseq/${noext}_${ref_pre}_uni_proper_mapped.bam


samtools view -bq 1 Results_fragseq/${noext}_${ref_pre}_mapped.bam > Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam


###########################
#remove reads with indels

samtools sort Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam -o Results_fragseq/${noext}_${ref_pre}_uni_mapped_sorted.bam
samtools index Results_fragseq/${noext}_${ref_pre}_uni_mapped_sorted.bam

python rm_indels_bam.py Results_fragseq/${noext}_${ref_pre}_uni_mapped_sorted.bam Results_fragseq/${noext}_${ref_pre}_umap_id_rm.bam 

rm Results_fragseq/${noext}_${ref_pre}_uni_mapped_sorted.bam
rm Results_fragseq/${noext}_${ref_pre}_mapped_0.bam
rm Results_fragseq/${noext}_${ref_pre}_uni_mapped_sorted.bam.bai


###########################



#in this uni_mapped file we can see that some reads were removed but their mates were saved
######## we need to remove reads without mates ##########
#extract ids of the reads in uni file

#samtools view Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam | cut -f 1 > Results_fragseq/ids_uni.txt


samtools view Results_fragseq/${noext}_${ref_pre}_umap_id_rm.bam | cut -f 1 > Results_fragseq/ids_uni.txt
python list_for_remove_reads.py Results_fragseq/ids_uni.txt Results_fragseq/ids_for_rm.txt

#samtools view -h Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam | grep -vf Results_fragseq/ids_for_rm.txt | samtools view -bS -o Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam -

#samtools view -h Results_fragseq/${noext}_${ref_pre}_umap_id_rm.bam | grep -vf Results_fragseq/ids_for_rm.txt | samtools view -bS -o Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam -
###
samtools view -N Results_fragseq/ids_for_rm.txt -U Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam -o /dev/null Results_fragseq/${noext}_${ref_pre}_umap_id_rm.bam

#####reports######

samtools flagstat Results_fragseq/inter.bam >> Results_fragseq/report 2>&1
echo "" >> Results_fragseq/report

echo "-only mapped reads--supplementary/secondary aligned reads were removed-" >> Results_fragseq/report
echo "" >> Results_fragseq/report

samtools flagstat Results_fragseq/${noext}_${ref_pre}_mapped.bam >> Results_fragseq/report 2>&1

echo "" >> Results_fragseq/report
echo "-uniquely properly mapped reads-" >> Results_fragseq/report
echo "" >> Results_fragseq/report

samtools flagstat Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam >> Results_fragseq/report 2>&1

echo "" >> Results_fragseq/report
echo "-REMOVE UNPAIRED READS after filtering the unique and removing indels-" >> Results_fragseq/report
echo "" >> Results_fragseq/report

samtools flagstat Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam >> Results_fragseq/report 2>&1


echo "" >> Results_fragseq/report
echo "-add pysam statistics about references-" >> Results_fragseq/report
echo "" >> Results_fragseq/report

####### plus pysam report ###########


samtools sort Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam -o Results_fragseq/${noext}_${ref_pre}_pr_map_sorted.bam
samtools index Results_fragseq/${noext}_${ref_pre}_pr_map_sorted.bam

python stats_pysam.py Results_fragseq/${noext}_${ref_pre}_pr_map_sorted.bam Results_fragseq/stat_pysam.txt 
cat Results_fragseq/stat_pysam.txt >> Results_fragseq/report 2>&1
echo "" >> Results_fragseq/report
echo "" >> Results_fragseq/report

rm Results_fragseq/inter.bam
rm Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam
#rm Results_fragseq/stat_pysam.txt

#Rscript --vanilla bam_to_csv.R Results_fragseq/${noext}_${ref_pre}_uni_mapped.bam Results_fragseq/${noext}_${ref_pre}_uni_mapped.csv
#Rscript --vanilla bam_to_csv.R Results_fragseq/Oct14_S10_R1_001_R2_kd403-real_uni_mapped.bam Results_fragseq/Oct14_S10_R1_001_R2_kd403-real_uni_mapped.csv


samtools view -h Results_fragseq/${noext}_${ref_pre}_proper_mapped.bam | \
  awk 'substr($0,1,1)=="@" || ($9>= 20 && $9<=40) || ($9<=-20 && $9>=-40)' | \
  samtools view -b > Results_fragseq/${noext}_${ref_pre}_proper_20_40.bam

samtools sort Results_fragseq/${noext}_${ref_pre}_proper_20_40.bam -o Results_fragseq/${noext}_${ref_pre}_proper_20_40_sorted.bam 
samtools index Results_fragseq/${noext}_${ref_pre}_proper_20_40_sorted.bam 
python stats_pysam.py Results_fragseq/${noext}_${ref_pre}_proper_20_40_sorted.bam Results_fragseq/stat_20_40_pysam.txt


####### plus pysam report after 20 40 lengths filter ###########
echo "-statistics for 20-40 nt -" >> Results_fragseq/report
samtools flagstat Results_fragseq/${noext}_${ref_pre}_proper_20_40.bam >> Results_fragseq/report 2>&1
echo "" >> Results_fragseq/report
echo "" >> Results_fragseq/report
cat Results_fragseq/stat_20_40_pysam.txt >> Results_fragseq/report 2>&1

#rm Results_fragseq/stat_20_40_pysam.txt




