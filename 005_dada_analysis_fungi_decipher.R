#Jason Bosch
#Analysis adapted from Gabriele Tosadori's scripts

#The data we are using are already demultiplexed and barcodes were already removed. So, the FASTQ files are ready to perform the DADA2 pipeline

# following the pipeline suggested here, for inferring ASVs
# https://benjjneb.github.io/dada2/tutorial.html

#####SET UP#####

#BioConductor libraries
library("dada2") #Denoising to ASVs
library("DECIPHER")
library("Biostrings")
library("phyloseq")
library("speedyseq")

###Directory###

setwd("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/")

#####LOAD DATA#####

# # import functions
# source("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/03_Analysis_Scripts/Gabriele/000_micro_functions_disko2013.R")

#get metadata
metadata <- read.csv("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/06_Miscellaneous/Disko2013-control_SW_metadata_for_Jason.csv")

# get forward and reverse reads
dna_fwds <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Fungi/DNA/", pattern="_R1_", full.names = TRUE))
dna_revs <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Fungi/DNA", pattern="_R2_", full.names = TRUE))
rna_fwds <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Fungi/RNA/", pattern="_R1_", full.names = TRUE))
rna_revs <- sort(list.files("~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/01_Raw_Data/Fungi/RNA", pattern="_R2_", full.names = TRUE))

#Combine and filter to relevant files
fn_fwds <- c(dna_fwds, rna_fwds)
fn_revs <- c(dna_revs, rna_revs)
fn_fwds <- fn_fwds[grepl(paste(metadata$sampleID_DNA_fungi,collapse = "|"),fn_fwds)|grepl(paste(metadata$sampleID_RNA_fungi,collapse = "|"),fn_fwds)]
fn_revs <- fn_revs[grepl(paste(metadata$sampleID_DNA_fungi,collapse = "|"),fn_revs)|grepl(paste(metadata$sampleID_RNA_fungi,collapse = "|"),fn_revs)]

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
# Using maxEE = 4 for max read retention.
out <- filterAndTrim(fwd=fn_fwds, filt=filt_fwds, rev=fn_revs, filt.rev=filt_revs, truncLen=c(220, 220), maxN=0, maxEE=c(4, 4), truncQ=2, rm.phix=FALSE, compress=TRUE, multithread=7)

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

# export the final dada2 table clean from chimeras
saveRDS(seq_tab_no_chim, "counts_table_ASV.Rds")

# export sequences as fasta file to extract ITS
itsxnm <- paste(">asv", seq(1, ncol(seq_tab_no_chim)), sep="_")
itsxas <- vector()

# loop through all the names
for (i in c(1:length(itsxnm))) {
  # and put together the table
  itsxas <- append(itsxas, rbind(itsxnm[i], colnames(seq_tab_no_chim)[i]))
}
# export the table
write.table(itsxas, "sequences_for_itsx.fasta", quote=F, col.names=F, row.names=F)

# ITS extraction
system("~/PostDoc/04_Software/ITSx_1.1.3/ITSx -i ~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/sequences_for_itsx.fasta -o ~/PostDoc/02_Projects/08_Disko_Island_Snow_Fence/02_Cleaned_Data/Fungi_ASV_OTU/ITS_extracted.fasta -t all -cpu 7")

# import the fasta file
its_extracted <- read.csv("./ITS_extracted.fasta.ITS2.fasta", header=F)

# get sequences only, from the fasta file
its_seqs <- its_extracted[seq(2, nrow(its_extracted), by=2), ]

# verify if all taxa have a match in the its_seqs file
all_matches <- vector()
col_names <- colnames(seq_tab_no_chim)
# loop through the its sequences
for (a_seq in its_seqs) {
	# find match
	all_matches <- append(all_matches, grep(a_seq, col_names))
}

# filter the seq_tab_no_chim table with ITSs sequences only
seq_tab_its <- seq_tab_no_chim[, unique(all_matches)]

# create phyloseq object out of seq_tab_its
phylo_tab <- phyloseq(otu_table(seq_tab_its, taxa_are_rows=F))

# get asv sequences that needs to be clustered and
# cast them into a format DECIPHER likes, i.e. DNAStringSet
dna <- Biostrings::DNAStringSet(taxa_names(phylo_tab))

# find clusters of ASVs to form the new OTUs
aligned_seqs <- AlignSeqs(dna, processors=6)

# compute distance matrix between the aligned sequences
dist_mat <- DistanceMatrix(aligned_seqs, processors=6)

# cluster together the sequences at 97% of similarity, i.e. cutoff=0.03
clusters <- TreeLine(myDistMatrix=dist_mat, method="complete", cutoff=0.03, type="clusters", processors=6)

# merge the ASVs into the OTUs that were found through ad-hoc function
# https://mikemc.github.io/speedyseq/reference/merge_taxa_vec.html
# https://github.com/mikemc/speedyseq/blob/main/NEWS.md#speedyseq-020
merged_seq_tab_as_phylo <- merge_taxa_vec(phylo_tab, group=clusters$cluster, tax_adjust=2)

# get the table that summarises how many reads were discarded for each filtering step
get_names <- function(x) sum(getUniques(x))
filter_info <- cbind.data.frame(sample_names, out, sapply(dada_fwds, get_names), sapply(dada_revs, get_names), sapply(mergers, get_names), rowSums(seq_tab_no_chim), rowSums(seq_tab_its), rowSums(merged_seq_tab_as_phylo))
colnames(filter_info) <- c("Sample","Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Non-chimera", "ITS Extracted", "Clustered 97%")
filter_info$Percent <- signif(filter_info$`Non-chimera`/filter_info$Input*100,4)
# export filtering information
write.table(filter_info, "dada_filtering_summary.csv", col.names = T, row.names = F, sep = ",")

# export table
saveRDS(merged_seq_tab_as_phylo, "counts_table_OTU_ITSx.Rds")

# get final counts table
counts_table <- otu_table(merged_seq_tab_as_phylo)

#set seed for reproducibility
set.seed(123456)

# run assignment for ASVs
#Using UNITE 10.0 All_Euk: https://doi.plutof.ut.ee/doi/10.15156/BIO/2959334
taxa_asv <- assignTaxonomy(counts_table, "~/PostDoc/04_Software/DADA2-classifier/sh_general_release_dynamic_all_04.04.2024.fasta", minBoot=70, outputBootstraps=TRUE, tryRC=TRUE, multithread=7)

# export taxonomy table
saveRDS(taxa_asv[[1]], "taxonomy_table_boot_70_NBC_OTU.Rds")
# and bootstrap scores
saveRDS(taxa_asv[[2]], "taxonomy_table_boot_70_NBC_OTU_scores.Rds")
