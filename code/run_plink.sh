#!/bin/bash

## vcf to bed format for PCAdapt

# ruderatus - aligned to hort
plink \
    --vcf ../data/ref_map_outputs/populations_vcf_rud/populations.snps.vcf \
    --make-bed `# convert to bed format` \
    --double-id `# plink wants to use _ to get family and sample id from sample name - tell it not to` \
    --allow-extra-chr `# ignore unrecognised chromosome codes` \
    --out ../data/plink/rud_sites_h

# # ruderatus - aligned to terr
# plink \
#     --vcf ../data/ref_map_outputs/pops_vcf_rud_t/populations.snps.vcf \
#     --make-bed \
#     --double-id \
#     --allow-extra-chr \
#     --out ../data/plink/rud_sites_t

# # hortorum - aligned to hort
# plink \
#     --vcf ../data/ref_map_outputs/populations_vcf_hort/populations.snps.vcf \
#     --make-bed \
#     --double-id \
#     --allow-extra-chr \
#     --out ../data/plink/hort_h

# # hortorum - aligned to terr
# plink \
#     --vcf ../data/ref_map_outputs/pops_vcf_hort_t/populations.snps.vcf \
#     --make-bed \
#     --double-id \
#     --allow-extra-chr \
#     --out ../data/plink/hort_t


# # terrestris
# plink \
#     --vcf ../data/ref_map_outputs/pops_vcf_terr_t/populations.snps.vcf \
#     --make-bed \
#     --double-id \
#     --allow-extra-chr \
#     --out ../data/plink/terr
