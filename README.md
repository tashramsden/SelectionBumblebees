# Signals of Selection in UK Bumblebees

**MSc research project**

**Languages**: Bash (version 5.0.17(1)), R (version 4.1.3), LaTeX (pdfTeX version 3.14159265-2.6-1.40.20)

**Dependencies**: `Stacks` (version 2.62), `BWA` (0.7.17), `SAMtools` (1.15.1), `PGDSpider` (2.1.1.5), `BayeScan` (2.1), `VCFtools` (0.1.16), `PLINK` (1.9)

R packages: `ggplot2`, `stringr`, `gridExtra`, `coda`, `dplyr`, `tidyr`, `viridis`, `readxls`, `cowplot`, `grDevices`, `viridis`, `scales`, `VennDiagram`, `ggpubr`, `pcadapt`, `ggsignif`



## Code files


### Bioinformatics Processing

#### Demultiplexing: Stacks process_radtags

* [**create_init_popmap.R**](code/create_init_popmap.R)
  * Extracts a starting popmap file (containing all 206 samples) from the original metadata.

* [**run_process_radtags.sh**](code/run_process_radtags.sh)
  * Bash script to be run on the cluster. Runs the stacks program process_radtags to demultiplex the raw RADseq data; samples with a phred score of less than 10 are discarded.

* [**make_reads_table.sh**](code/make_reads_table.sh)
  * Extracts read count information for each sample from process_radtags for assaying and saves to reads_table.tsv

* [**assay_process_radtags.R**](code/assay_process_radtags.R)
  * Assaying process_radtags. Takes process_radtags log file. Checks contribution of each sample to total library and assesses retained reads. Outputs an updated reads table, metadata, popmap and assay figures.

#### Aligning sequences to reference genomes: bwa and samtools

* [**index_reference_genomes.sh**](code/index_reference_genomes.sh)
  * Indexing the reference genomes to create searchable databases to be used by the aligner in the next step.

* [**align_hort_bam.sh**](code/align_hort_bam.sh) and [**align_terr_bam.sh**](code/align_terr_bam.sh)
  * Bash script to be run on the cluster. Aligns each sample (output from process_radtags) to the B. hortorum or B. terrestris reference genome using BWA mem, produces SAM files. Converts SAM to BAM files using Samtools. Creates csv file containing the percentage of each sample successfully mapped to the reference genome.

* [**assay_bwa_alignments.R**](code/assay_bwa_alignments.R)
  * Assaying the alignments from BWA and Samtools. Takes the file created by align_hort_bam.sh of the percentages of samples successfully mapped to the reference genomes. Produces figures of how well samples form each species align to the reference genomes, updates metadata.

#### Stacks ref_map

* [**make_final_popmap_metadata.R**](code/make_final_popmap_metadata.R)
  * Collates samples to be removed based on previous analyses and assaying; updates metadata; creates popmap for ref_map pipeline.

* [**run_ref_map_hort.sh**](code/run_ref_map_hort.sh) and [**run_ref_map_terr.sh**](code/run_ref_map_terr.sh)
  * Bash scripts to be run on the cluster. Runs the stacks ref_map pipeline to create a catalogue of genotyped RAD loci for each species, aligned to the B. hortorum or B. terrestris reference genome.

* [**make_refmap_assay_tables.sh**](code/make_refmap_assay_tables.sh)
  * Extracts information on retained alignments, coverage depth and phasing from the ref_map gstacks log file and saves to tsvs for assaying.

* [**assay_ref_map.R**](code/assay_ref_map.R)
  * Assaying ref_map pipeline. Takes tsv files created by make_refmap_assay_tables. Produces figures exploring retained alignments, coverage depth and phasing.


### Signals of Selection

#### Nucleotide diversity selection analyses

* [**run_VCFtools.sh**](code/run_VCFtools.sh)
  * Runs VCFtools' --SNPdensity to find the number of SNPs in windows of a given size; used to help determine minimum window size for kernel smoothing in Stacks.

* [**find_window_size.R**](code/find_window_size.R)
  * Takes output tsvs from run_VCFtools.sh to find the mean number of SNPs per window for different window sizes.

* [**run_psmooth_rud_h.sh**](code/run_psmooth_rud_h.sh), [**run_psmooth_hort_h.sh**](code/run_psmooth_hort_h.sh) and [**run_psmooth_terr_h.sh**](code/run_psmooth_terr_h.sh)
  * Bash scripts to be run on the cluster. Runs Stacks' populations programme to calculate kernel smoothed statistcis and calculate significance of nucleotide diversity outliers using bootstrap resampling.

* [**pops_smoothed_analysis_rud.R**](code/pops_smoothed_analysis_rud.R), [**pops_smoothed_analysis_hort.R**](code/pops_smoothed_analysis_hort.R) and [**pops_smoothed_analysis_terr.R**](code/pops_smoothed_analysis_terr.R)
  * Explores nucleotide diversity in each species and creates Manhattan-style plots of nucleotide diversity across the species' genomes with outlier regions highlighted. In pops_smoothed_analysis_rud.R Fst is additionally explored and plotted.

* [**get_genes_from_flat.sh**](code/get_genes_from_flat.sh)
  * Quick bash script to extract genes found in nucleotide diversity outlier regions from flat files (flat files manually extracted from NCBI Data Viewer); used in pops_smoothed_analysis scripts.

* [**gene_annotation_analysis.R**](code/gene_annotation_analysis.R)
  * Comparing nucleotide diversity among the three species and genes found in outlier regions of the three species; creates venn diagram and pie charts plots. 

#### Additional tests for B. ruderatus

* [**run_populations.sh**](code/run_populations.sh)
  * Runs the stacks populations programme to create VCF data for each species separately. Filters data according to requirements for downstream analysis. 

* [**run_PGDSpider.sh**](code/run_PGDSpider.sh)
  * Runs PGDSpider to convert files from Stacks' populations output (VCF) into the format ready for Bayescan.

* [**run_bayescan.sh**](code/run_bayescan.sh)
  * Bash script to be run on the cluster. Runs BayeScan to detect Fst outliers in B. ruderatus based on population differentiation between the geographic sites. Runs with default parameters and increased thinning interval.

* [**get_snp_ids.sh**](code/get_snp_ids.sh)
  * Quick bash script to extract SNP IDs from populations' VCF files. Used to find IDs of BayeScan and PCAdapt outlier SNPs, saves IDs to text files; outputs used in bayescan_analysis.R and run_PCAdapt.R.

* [**get_chr_bp_from_outliers**](code/get_chr_bp_from_outlier.sh)
  * Quick bash script to extract chromosome and basepair positions of outlier loci (used in bayescan_analysis.R and run_PCAdapt.R).

* [**bayescan_analysis.R**](code/bayescan_analysis.R)
  * Inspects BayeScan results, evaluating MCMC chains for convergence and efficiency. Saves data about alleles at outlier SNPs ready for making heatmaps.

* [**run_plink.sh**](code/run_plink.sh)
  * Runs PLINK to convert populations' VCF files into bed format ready for PCAdapt. 

* [**run_PCAdapt.R**](code/run_PCAdapt.R)
  * Runs PCAdapt to detect outlier SNPs in B. ruderatus based on population structure identified using principal component analysis. Saves data about alleles at outlier SNPs ready for making heatmaps.

* [**make_heatmaps.R**](code/make_heatmaps.R)
  * Calculates dominant allele frequency at BayeScan and PCAdapt outlier SNPs. Plots a heatmap of dominant allele frequency at these outliers across the geographic locations, alongside a table of how many tests detected each outlier. 



## Thesis folder

The LaTeX code files used to create my full report including supplementary material are available in the [thesis](thesis) folder, as well as the bibtex file (containing all citations, including the 58 used for gene annotation) and a bash script to compile the report. Figures for the report can be created using the code files above.



## Author

Natasha Ramsden | tash.ramsden21@imperial.ac.uk
