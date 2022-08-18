#!/bin/bash
# gets snp ids from vcf files in order passed to bayescan

dir=../data/ref_map_outputs/

# ruderatus
grep -v "#" $dir/populations_vcf_rud/populations.snps.vcf | cut -f 3 > $dir/populations_vcf_rud/rud_snp_ids.txt

# hortorum
grep -v "#" $dir/populations_vcf_hort/populations.snps.vcf | cut -f 3 > $dir/populations_vcf_hort/hort_snp_ids.txt


# ## for aligned to terr ref genome

# # ruderatus
# grep -v "#" $dir/pops_vcf_rud_t/populations.snps.vcf | cut -f 3 > $dir/pops_vcf_rud_t/rud_snp_ids_t.txt

# # hortorum
# grep -v "#" $dir/pops_vcf_hort_t/populations.snps.vcf | cut -f 3 > $dir/pops_vcf_hort_t/hort_snp_ids_t.txt

# terrestris
grep -v "#" $dir/pops_vcf_terr_t/populations.snps.vcf | cut -f 3 > $dir/pops_vcf_terr_t/terr_snp_ids.txt
