#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=15gb

module load anaconda3/personal
cd $HOME/bee_project/code

## Smoothed stats
## Without --write-single-snp to get full snp data:

path=../data/ref_map_outputs
smooth_path=../data/populations_smoothed

hort_popmap=../data/metadata/popmap_hort_67n_site.csv
rud_popmap=../data/metadata/popmap_rud_56n_site.csv
terr_popmap=../data/metadata/popmap_terr_55n_site.csv

touch $smooth_path/stopwatch.txt

sigmas=(100000 125000 150000 175000 200000)

## local testing
# test=5
# index=$(expr $test - 1)

# -J 1-5
index=$(expr $PBS_ARRAY_INDEX - 1)
sigma=${sigmas[$index]}

# sigma=150000		

echo -e "sigma = ${sigma}" >> $smooth_path/stopwatch.txt
label="${sigma:0:3}"

mkdir $smooth_path/terr_t_smooth_sigma$label

		
## terr - vcf - smoothed stats - first run - generate bootstrap archive
echo -e "start vcf terr t smooth: `date`" >> $smooth_path/stopwatch.txt
/rds/general/user/ter21/home/stacks-2.62/populations \
	--in-path $path/refmap_terr_178n 		`# input dir - stacks aligned to hort ref genome` \
	--popmap $terr_popmap					`# popmap` \
	--threads 32 							`# num threads` \
	--min-samples-per-pop 0.5				`# loci in min 50% samples for each pop` \
	--hwe 									`# calculate hwe for loci` \
	--smooth								`# enable kernel-smoothed pi, Fis, Fst stats` \
	--bootstrap-archive						`# archive stats for use in bootstrap resampling in subsequent run` \
	--sigma $sigma							`# standard deviation of the kernel smoothing weight distribution: half sliding window = 3*sigma` \
	--vcf-all									`# output results in VCF format` \
	--out-path $smooth_path/terr_t_smooth_sigma$label	`# output dir`
echo -e "end vcf terr t smooth: `date`" >> $smooth_path/stopwatch.txt

## second run - use bootstrap archive to bootstrap resample
echo -e "start vcf terr t smooth 2 w 10000 reps: `date`" >> $smooth_path/stopwatch.txt
/rds/general/user/ter21/home/stacks-2.62/populations \
	--in-path $path/refmap_terr_178n 		`# input dir - stacks aligned to hort ref genome` \
	--popmap $terr_popmap					`# popmap` \
	--threads 32 							`# num threads` \
	--min-samples-per-pop 0.5				`# loci in min 50% samples for each pop` \
	--smooth								`# enable kernel-smoothed pi, Fis, Fst stats` \
	--bootstrap								`# boostrap resampling for all kernel-smoothed stats` \
	--bootstrap-reps 10000 					`# increase num bootstrap reps` \
	--sigma $sigma							`# standard deviation of the kernel smoothing weight distribution: half sliding window = 3*sigma` \
	--vcf-all									`# output results in VCF format` \
	--out-path $smooth_path/terr_t_smooth_sigma$label	`# output dir`
echo -e "end vcf terr t smooth 2: `date`" >> $smooth_path/stopwatch.txt
