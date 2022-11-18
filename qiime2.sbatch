#!/bin/bash -l

#sbatch script for running qiime2 on Wahab HPC at ODU
#set all variables in config.qiime file
#comment in or out modules below as needed

#SBATCH --job-name=qiime2      
#SBATCH -o qiime2-%j.out
#SBATCH --time=00:00:00

enable_lmod
module load container_env qiime2/2022.08


#this will trim the fq.gz files using the settings in config.4.all
#it is assumed that the directory you run this script from has the 
#fq.gz files
bash qiime2.bash 1_setup_qiime2 config.qiime2

#this will assemble the fq.gz files in the mkREF and mkBAM directories that were
#created by the trim function.  We have to cd to the mkREF directory
#then execute dDocentHPC mkREF
#bash qiime2.bash 1a_feature_classifier_qiime2 config.qiime2

#this will use dDocent mkBAM to map the reads in the mkBAM directory to the reference*fasta created by 
#dDocent mkREF.  We have to change directories, copy the reference*fasta to mkBAM, then run dDocent mkBAM 
#bash qiime2.bash 2_file_import_qiime2 config.qiime2

#this will use dDocent fltrBAM to filter the BAM files
#bash qiime2.bash 3_preprocessing_qiime2 config.qiime2

#this will use freebayes to genotype the bam files and make a VCF
#bash qiime2.bash 4_filter_coremetrics_qiime2 config.qiime2

#bash qiime2.bash 5_univariate_stats_qiime2 config.qiime2
#bash qiime2.bash 6_multivarite stats_qiime2 config.qiime2
#bash qiime2.bash 7_ancom_qiime2 config.qiime2


#bash qiime2.bash 4_filter_coremetrics_qiime2 config.qiime2

