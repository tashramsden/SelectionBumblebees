#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=5gb

# Uses BWA to align the sequences (output from process_radtags) to the B. terrestris genome. Then uses samtools to convert the output sam files to bam files and sort them.
# Requires: bwa and samtools software packages

module load anaconda3/personal
cd $HOME/bee_project/data/

echo -e "Start: `date`" >> stopwatch_align.txt

# file_index=1  # for local testing
file_index=$PBS_ARRAY_INDEX  # 1-204

file=$(ls process_radtags_outputs/assayed_results_204n | head -n $file_index | tail -n 1)
sample=$(ls process_radtags_outputs/assayed_results_204n | head -n $file_index | tail -n 1 | awk -F "." '{print $1}')

# echo $file
echo $sample >> stopwatch_align.txt

echo "Aligning sample to reference genome..."
# bwa mem to align the sequence to the terrestris reference genome
bwa mem align_seqs_bwa/indexed_terrestris/terrestris_index process_radtags_outputs/assayed_results_204n/$file > align_seqs_bwa/aligned_to_terr_sam/$sample.sam

echo "Converting to bam file and sorting..."
# samtools to convert the sam to bam files and sort them
samtools view -b align_seqs_bwa/aligned_to_terr_sam/$sample.sam | samtools sort -O bam > align_seqs_bwa/aligned_to_terr_bam/$sample.bam
# -b for view = output to bam format
# -O for sort = output to bam format

echo "Saving percentage of sequence mapped to ref genome..."
# samtools to find the percentage of sample mapped to the ref genome
samtools flagstats align_seqs_bwa/aligned_to_terr_bam/$sample.bam > align_seqs_bwa/flagstats_terr/$sample.txt
mapped=$(egrep -o '[0-9]+\.[0-9]+%' align_seqs_bwa/flagstats_terr/$sample.txt)
echo $sample $mapped >> align_seqs_bwa/percent_mapped_to_terrestris.csv

echo "Finished sample"

echo -e "Finish: `date`" >> stopwatch_align.txt
