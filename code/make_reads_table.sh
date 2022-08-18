#!/bin/sh
# Takes the process_radtags log file and creates a reads_table.tsv 

log_file=../data/process_radtags_outputs/results/process_radtags.FASTQ_sequences.log

# start_table=$(grep -n "BEGIN per_file_raw_read_counts" $log_file | cut -d: -f1)
# end_table=$(grep -n "END per_file_raw_read_counts" $log_file | cut -d: -f1)
# # get line num of last line of file to include (ie end of the reads table)
# last_table_inc=$(expr $end_table - 2)  # not include table end line or control sample
# # get length of table
# table_length=$(expr $last_table_inc - $start_table)

# # trim file to just include the table and save to reads_table.tsv
# head -$last_table_inc $log_file | tail -$table_length > ../data/process_radtags_outputs/reads_table.tsv

stacks-dist-extract $log_file per_file_raw_read_counts > ../data/process_radtags_outputs/reads_table.tsv
