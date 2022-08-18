#!/bin/sh
# Takes the ref_map gstacks.log.distribs file to extract key assay tables to separate tsv files. Repeats this for samples aligned to hort and terr ref genomes.

hort_path=../data/ref_map_outputs/refmap_hort_178n
terr_path=../data/ref_map_outputs/refmap_terr_178n

for path in $hort_path $terr_path;
    do

        echo $path

        log_file=$path/gstacks.log.distribs

        # verify proportion of retained reads
        echo "gstacks.log content:"
        start_table=$(grep -n "BAM records:" $path/gstacks.log | cut -d: -f1)
        end_table=$(grep -n "skipped some suboptimal" $path/gstacks.log | cut -d: -f1)
        table_length=$(expr $end_table - $start_table + 1)
        head -$end_table $path/gstacks.log | tail -$table_length

        # bam stats
        stacks-dist-extract $log_file bam_stats_per_sample > $path/bam_stats.tsv

        # coverage depth (remove first 2 lines - comment)
        stacks-dist-extract $log_file effective_coverages_per_sample | tail +3 > $path/coverage.tsv

        # phasing rates
        stacks-dist-extract $log_file phasing_rates_per_sample > $path/phasing.tsv

    done
exit;
