#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=2gb

module load anaconda3/personal
cd $HOME/bee_project/data/

num_chains=$PBS_ARRAY_INDEX  # to repeat each run of bayescan -> check convergence of chains

touch bayescan/stopwatch.txt


# ruderatus ----

# aligned to hortorum ref genome 
echo -e "Start rud $num_chains: `date`" >> bayescan/stopwatch.txt
rud_input=PGDSpider/rud_bayescan/pgds_vcf_bayescan_rud
rud_output_dir=bayescan/rud_sites
# default params
bayescan_2.1 $rud_input -od $rud_output_dir -o rud_sites$num_chains -threads 32 -out_pilot -out_freq
# higher thinning interval
bayescan_2.1 $rud_input -od $rud_output_dir -o rud_sites_thin50_$num_chains -threads 32 -out_pilot -out_freq -thin 50
# with pr_odds higher 100
# bayescan_2.1 $rud_input -od $rud_output_dir -o rud_sites_prodds100_$num_chains -threads 32 -out_pilot -out_freq -pr_odds 100
echo -e "Finish rud $num_chains: `date`" >> bayescan/stopwatch.txt

# # aligned to terrestris ref genome
# echo -e "Start rud alignt $num_chains: `date`" >> bayescan/stopwatch.txt
# rud_input_t=PGDSpider/rud_bayescan/pgds_vcf_bayescan_rud_alignt
# rud_output_dir_t=bayescan/rud_sites_alignt
# # default params
# bayescan_2.1 $rud_input_t -od $rud_output_dir_t -o rud_sites$num_chains -threads 32 -out_pilot -out_freq
# # higher thinning interval
# bayescan_2.1 $rud_input_t -od $rud_output_dir_t -o rud_sites_thin50_$num_chains -threads 32 -out_pilot -out_freq -thin 50
# # with pr_odds higher 100
# bayescan_2.1 $rud_input_t -od $rud_output_dir_t -o rud_sites_prodds100_$num_chains -threads 32 -out_pilot -out_freq -pr_odds 100
# echo -e "Finish rud $num_chains: `date`" >> bayescan/stopwatch.txt


# # hortorum ----

# # aligned to hortorum ref genome
# echo -e "Start hort $num_chains: `date`" >> bayescan/stopwatch.txt
# hort_input=PGDSpider/hort_bayescan/pgds_vcf_bayescan_hort
# hort_output_dir=bayescan/hort_sites
# # default params
# bayescan_2.1 $hort_input -od $hort_output_dir -o hort_sites$num_chains -threads 32 -out_pilot -out_freq
# # higher thinning interval
# bayescan_2.1 $hort_input -od $hort_output_dir -o hort_sites_thin50_$num_chains -threads 32 -out_pilot -out_freq -thin 50
# # with pr_odds higher 100
# bayescan_2.1 $hort_input -od $hort_output_dir -o hort_sites_prodds100_$num_chains -threads 32 -out_pilot -out_freq -pr_odds 100
# echo -e "Finish hort $num_chains: `date`" >> bayescan/stopwatch.txt

# # aligned to terrestris ref genome
# echo -e "Start hort alignt $num_chains: `date`" >> bayescan/stopwatch.txt
# hort_input_t=PGDSpider/hort_bayescan/pgds_vcf_bayescan_hort_alignt
# hort_output_dir_t=bayescan/hort_sites_alignt
# # default params
# bayescan_2.1 $hort_input_t -od $hort_output_dir_t -o hort_sites$num_chains -threads 32 -out_pilot -out_freq
# # higher thinning interval
# bayescan_2.1 $hort_input_t -od $hort_output_dir_t -o hort_sites_thin50_$num_chains -threads 32 -out_pilot -out_freq -thin 50
# # with pr_odds higher 100
# bayescan_2.1 $hort_input_t -od $hort_output_dir_t -o hort_sites_prodds100_$num_chains -threads 32 -out_pilot -out_freq -pr_odds 100
# echo -e "Finish hort $num_chains: `date`" >> bayescan/stopwatch.txt


# # terrestris ----

# echo -e "Start terr $num_chains: `date`" >> bayescan/stopwatch.txt
# terr_input=PGDSpider/terr_bayescan/pgds_vcf_bayescan_terr
# terr_output_dir=bayescan/terr_sites_alignt
# # default params
# bayescan_2.1 $terr_input -od $terr_output_dir -o terr_sites$num_chains -threads 32 -out_pilot -out_freq
# # higher thinning interval
# bayescan_2.1 $terr_input -od $terr_output_dir -o terr_sites_thin50_$num_chains -threads 32 -out_pilot -out_freq -thin 50
# # with pr_odds higher 100
# bayescan_2.1 $terr_input -od $terr_output_dir -o terr_sites_prodds100_$num_chains -threads 32 -out_pilot -out_freq -pr_odds 100
# echo -e "Finish terr $num_chains: `date`" >> bayescan/stopwatch.txt
