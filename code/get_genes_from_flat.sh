#!/bin/bash
# script to use within R scripts (pops_smoothed_analysis) to get genes from each outlier region from flat files (NCBI)

# smooth_path="../data/populations_smoothed"
# flat_locs=$smooth_path/*flat_genes*

# for folder in $flat_locs;
#     do

#         touch $folder/flat_genes.txt

#         for f in $folder/*.flat;
#             do
#                 grep "\gene=" $f >> $folder/flat_genes.txt
#             done

#     done
# exit;

# do each genomic region separately - can keep track of which genes from where
# smooth_path="../data/populations_smoothed"
flat_loc=../data/populations_smoothed/flat_genes_hort_h/one_at_a_time
touch $flat_loc/flat_genes.txt
for f in $flat_loc/*.flat;
    do
        grep "\gene=" $f >> $flat_loc/flat_genes.txt
    done
cat $flat_loc/flat_genes.txt
