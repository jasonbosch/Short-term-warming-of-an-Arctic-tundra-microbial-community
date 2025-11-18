#Jason Bosch
#Analysis adapted from Gabriele Tosadori's scripts

#The data we are using are already demultiplexed and barcodes were already removed. So, the FASTQ files are ready to perform the DADA2 pipeline

# following the pipeline suggested here, for inferring ASVs
# https://benjjneb.github.io/dada2/tutorial.html

#####SET UP#####

#BioConductor libraries
library(dada2) #Denoising to ASVs

###Directory###

setwd("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Bacteria_ASV/")

#####LOAD DATA#####

# import functions
source("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/03_Analysis_Scripts/Gabriele/000_micro_functions_disko2013.R")

#get metadata
metadata <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/Disko2013-control_SW_metadata_for_Jason.csv")

# get forward and reverse reads
dna_fwds <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Bacteria/DNA/", pattern="_R1_", full.names = TRUE))
dna_revs <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Bacteria/DNA", pattern="_R2_", full.names = TRUE))
rna_fwds <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Bacteria/RNA/", pattern="_R1_", full.names = TRUE))
rna_revs <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Bacteria/RNA", pattern="_R2_", full.names = TRUE))

#Combine and filter to relevant files
fn_fwds <- c(dna_fwds, rna_fwds)
fn_revs <- c(dna_revs, rna_revs)
fn_fwds <- fn_fwds[grepl(paste(metadata$sampleID_DNA_bacteria,collapse = "|"),fn_fwds)|grepl(paste(metadata$sampleID_RNA_bacteria,collapse = "|"),fn_fwds)]
fn_revs <- fn_revs[grepl(paste(metadata$sampleID_DNA_bacteria,collapse = "|"),fn_revs)|grepl(paste(metadata$sampleID_RNA_bacteria,collapse = "|"),fn_revs)]

#####ANALYSIS#####

# extract sample names
sample_names <- vector()
# for each line in fn_fwds
for (ln in c(1:length(fn_fwds))) {
  # split the string and get the first occurrence. also remove the reference to the read direction, i.e. R1 in the case of fn_fwds
  sample_names <- append(sample_names, strsplit(strsplit(basename(fn_fwds[ln]), "\\.")[[1]][1], "_R1")[[1]][1])
}

# where filtered reads will be placed
filt_fwds <- paste0("filtered_forward/", paste0(sample_names, "_F_filt.fastq"))
filt_revs <- paste0("filtered_reverse/", paste0(sample_names, "_R_filt.fastq"))

#Check quality scores
f_qual <- plotQualityProfile(fn_fwds,aggregate = T)
svg("read_quality_forward.svg", width = 2*84*0.04, height = 2*84*0.04)
f_qual
dev.off()
r_qual <- plotQualityProfile(fn_revs,aggregate = T)
svg("read_quality_reverse.svg", width = 2*84*0.04, height = 2*84*0.04)
r_qual
dev.off()

# filter reads based on quality
# Manual selection of truncation positions
# Allow more expected errors on forward reads because of lower quality. Many more are kept.
out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncLen=c(220, 220), maxN=0, maxEE=c(4, 2), truncQ=2, rm.phix=FALSE, compress=TRUE, multithread=7)

# learning the error rates
err_fwd <- learnErrors(filt_fwds, multithread=7)
err_rev <- learnErrors(filt_revs, multithread=7)

# and plot the errors, to verify if everything is ok.
plot_err_fwd <- plotErrors(err_fwd, nominalQ=TRUE)
svg("errors_forward.svg", width = 2*84*0.04, height = 2*84*0.04)
plot_err_fwd
dev.off()
plot_err_rev <- plotErrors(err_rev, nominalQ=TRUE)
svg("errors_reverse.svg", width = 2*84*0.04, height = 2*84*0.04)
plot_err_rev
dev.off()

# perform sample inference
dada_fwds <- dada(filt_fwds, err=err_fwd, multithread=7, pool="pseudo")
dada_revs <- dada(filt_revs, err=err_rev, multithread=7, pool="pseudo")

# Merge paired end reads
mergers <- mergePairs(dada_fwds, filt_fwds, dada_revs, filt_revs, verbose=FALSE, returnRejects=FALSE, minOverlap=12, maxMismatch=0)

# create the sequence table
seq_tab <- makeSequenceTable(mergers)

# export seq_tab, before removing chimeras
saveRDS(seq_tab, "seq_tab.Rds")

# remove chimeras
seq_tab_no_chim <- removeBimeraDenovo(seq_tab, method="consensus", multithread=7, verbose=FALSE)

# get the table that summarises how many reads were discarded for each filtering step
get_names <- function(x) sum(getUniques(x))
filter_info <- cbind.data.frame(sample_names, out, sapply(dada_fwds, get_names), sapply(dada_revs, get_names), sapply(mergers, get_names), rowSums(seq_tab_no_chim))
colnames(filter_info) <- c("Sample","Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Non-chimera")
filter_info$Percent <- signif(filter_info$`Non-chimera`/filter_info$Input*100,4)
# export filtering information
write.table(filter_info, "dada_filtering_summary.csv", col.names = T, row.names = F, sep = ",")

# export the final dada2 table clean from chimeras
saveRDS(seq_tab_no_chim, "counts_table_ASV.Rds")

#set seed for reproducibility
set.seed(123456)

# run assignment for ASVs
#Using dada2 silva 138.1 sequences: https://zenodo.org/records/4587955
taxa_asv <- assignTaxonomy(seq_tab_no_chim, "~/PostDoc/04_Software/DADA2-classifier/silva_nr99_v138.1_train_set.fa.gz", minBoot=70, outputBootstraps=TRUE, tryRC=TRUE, multithread=7)

# export taxonomy table
saveRDS(taxa_asv[[1]], "taxonomy_table_boot_70_NBC_ASV.Rds")
# and bootstrap scores
saveRDS(taxa_asv[[2]], "taxonomy_table_boot_70_NBC_ASV_scores.Rds")

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: x86_64-pc-linux-gnu
# Running under: Arch Linux
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblas.so.0.3;  LAPACK version 3.12.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Prague
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dada2_1.32.0 Rcpp_1.0.12 
# 
# loaded via a namespace (and not attached):
# [1] SummarizedExperiment_1.34.0 gtable_0.3.5                hwriter_1.3.2.1             ggplot2_3.5.1               latticeExtra_0.6-30        
# [6] Biobase_2.63.1              lattice_0.22-6              vctrs_0.6.5                 tools_4.4.1                 bitops_1.0-7               
# [11] generics_0.1.3              stats4_4.4.1                parallel_4.4.1              tibble_3.2.1                fansi_1.0.6                
# [16] pkgconfig_2.0.3             Matrix_1.7-0                RColorBrewer_1.1-3          S4Vectors_0.41.7            RcppParallel_5.1.7         
# [21] lifecycle_1.0.4             GenomeInfoDbData_1.2.12     compiler_4.4.1              farver_2.1.2                stringr_1.5.1              
# [26] deldir_2.0-4                Rsamtools_2.20.0            Biostrings_2.72.1           munsell_0.5.1               codetools_0.2-20           
# [31] GenomeInfoDb_1.39.14        pillar_1.9.0                crayon_1.5.2                BiocParallel_1.38.0         DelayedArray_0.30.1        
# [36] ShortRead_1.62.0            abind_1.4-5                 tidyselect_1.2.1            stringi_1.8.4               dplyr_1.1.4                
# [41] reshape2_1.4.4              labeling_0.4.3              grid_4.4.1                  colorspace_2.1-0            cli_3.6.2                  
# [46] SparseArray_1.4.8           magrittr_2.0.3              S4Arrays_1.4.1              utf8_1.2.4                  withr_3.0.0                
# [51] scales_1.3.0                UCSC.utils_0.99.7           pwalign_1.0.0               XVector_0.43.1              httr_1.4.7                 
# [56] matrixStats_1.3.0           jpeg_0.1-10                 interp_1.1-6                png_0.1-8                   GenomicRanges_1.56.1       
# [61] IRanges_2.37.1              rlang_1.1.3                 glue_1.7.0                  BiocGenerics_0.49.1         rstudioapi_0.16.0          
# [66] jsonlite_1.8.8              R6_2.5.1                    plyr_1.8.9                  MatrixGenerics_1.16.0       GenomicAlignments_1.40.0   
# [71] zlibbioc_1.49.3       