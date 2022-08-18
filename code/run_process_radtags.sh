#!/bin/bash
#PBS -l walltime=7:00:00
#PBS -l select=1:ncpus=1:mem=1gb
module load anaconda3/personal
cd $HOME/bee_project/data/
touch process_radtags_outputs/stopwatch.txt  # make a stopwatch
echo -e "Start date: `date`" >> process_radtags_outputs/stopwatch.txt  # start time
process_radtags \
    -p FASTQ_sequences `# directory containing g-zipped fastq files` \
    -o process_radtags_outputs/results `# output directory` \
    -e SgrAI `# restriction enzyme` \
    -c `# clean data, remove any read with an uncalled base` \
    -q  `# discard reads with low quality scores (phred < 10)` \
    -r `# rescue barcodes and RAD-Tags` \
    -D `# capture discarded reads to a file` \
    -i gzfastq `# input file type` \
    &> process_radtags_outputs/process_radtags.oe `# save outputs here`
echo -e "Finish date: `date`" >> process_radtags_outputs/stopwatch.txt  # finish time
