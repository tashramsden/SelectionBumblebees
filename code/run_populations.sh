#!/bin/bash

path=../data/ref_map_outputs

hort_popmap=../data/metadata/popmap_hort_67n_site.csv
rud_popmap=../data/metadata/popmap_rud_56n_site.csv
terr_popmap=../data/metadata/popmap_terr_55n_site.csv

rud_clumpak_subpop_popmap=../data/metadata/popmap_rud_56n_clumpak_subpops.csv

touch $path/pop_stopwatch.txt

# rud - aligned to hort ref genome - vcf for bayescan
echo -e "start vcf just rud: `date`" >> $path/pop_stopwatch.txt
populations \
	--in-path $path/refmap_hort_178n 		`# input dir - stacks aligned to hort ref genome` \
	--popmap $rud_popmap					`# popmap` \
	--threads 5 							`# num threads` \
	--min-samples-per-pop 0.5				`# loci in min 50% samples for each pop` \
	--min-maf 0.05							`# minimum minor allele frequency required to process a nucleotide site at a locus` \
	--hwe 									`# calculate hwe for loci` \
	--fstats								`# enable SNP-based F statistics` \
	--write-single-snp 						`# only retain first snp per locus` \
	--vcf									`# output results in VCF format` \
	--out-path $path/populations_vcf_rud_test	`# output dir`
echo -e "end vcf just rud: `date`" >> $path/pop_stopwatch.txt

# # hort - vcf
# echo -e "start vcf just hort: `date`" >> $path/pop_stopwatch.txt
# populations \
# 	--in-path $path/refmap_hort_178n 		`# input dir - stacks aligned to hort ref genome` \
# 	--popmap $hort_popmap					`# popmap` \
# 	--threads 5 							`# num threads` \
# 	--min-samples-per-pop 0.5				`# loci in min 50% samples for each pop` \
# 	--min-maf 0.05							`# minimum minor allele frequency required to process a nucleotide site at a locus` \
# 	--hwe 									`# calculate hwe for loci` \
# 	--fstats								`# enable SNP-based F statistics` \
# 	--write-single-snp 						`# only retain first snp per locus` \
# 	--vcf									`# output results in VCF format` \
# 	--out-path $path/populations_vcf_hort	`# output dir`
# echo -e "end vcf just hort: `date`" >> $path/pop_stopwatch.txt

# terr - vcf - aligned to terrestris genome
# echo -e "start vcf just terr: `date`" >> $path/pop_stopwatch.txt
# populations \
# 	--in-path $path/refmap_terr_178n 		`# input dir - stacks aligned to terr ref genome` \
# 	--popmap $terr_popmap					`# popmap` \
# 	--threads 5 							`# num threads` \
# 	--min-samples-per-pop 0.5				`# loci in min 50% samples for each pop` \
# 	--min-maf 0.05							`# minimum minor allele frequency required to process a nucleotide site at a locus` \
# 	--hwe 									`# calculate hwe for loci` \
# 	--fstats								`# enable SNP-based F statistics` \
# 	--write-single-snp 						`# only retain first snp per locus` \
# 	--vcf									`# output results in VCF format` \
# 	--out-path $path/pops_vcf_terr_t	`# output dir`
# echo -e "end vcf just terr: `date`" >> $path/pop_stopwatch.txt
