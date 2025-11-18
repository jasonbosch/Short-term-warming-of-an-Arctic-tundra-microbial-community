##Preparation##

All sequencing data associated with the project can be downloaded from https://www.ncbi.nlm.nih.gov/sra/PRJNA1264795

To demultiplex and prepare the files, please run the following scripts, available from https://github.com/gabrielet/tosadori_disko:  
000_parsing_FASTQs.py  
002_demultiplex_samples.R  
003_run_cutadapt.R  

You will also need the 000_micro_functions_disko2013.R script from https://github.com/gabrielet/tosadori_disko.

The preceeding scripts were described in Tosadori et al., under review.

Other required data (metadata, temperature recordings, fungal guild identification) and the final R workspace are available on Zenodo: https://doi.org/10.5281/zenodo.17610846

##Analysis##

Run the following scripts:  
005_dada2_analysis_bacteria.R  
005_dada_analysis_fungi_decipher.R  
006_Disko_Analysis_Final.R

##Note##

You will need to update the folder locations and pointers to external scripts.
