#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=32:mem=15gb
module load anaconda3/personal
cd $HOME/bee_project/data/
touch ref_map_outputs/stopwatch.txt  # make a stopwatch
echo -e "178n terr Start date: `date`" >> ref_map_outputs/stopwatch.txt  # start time
ref_map.pl \
    -T 32 `# num CPUs` \
    --samples align_seqs_bwa/aligned_to_terr_bam/ `# location of BAM files` \
    --popmap metadata/popmap_178n_spp.csv `# population map file: sample names to be used and their species` \
    -o ref_map_outputs/refmap_terr_178n/ `# path to output file directory`
echo -e "Finish date: `date`" >> ref_map_outputs/stopwatch.txt  # finish time
