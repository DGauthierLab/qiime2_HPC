#!/bin/bash -l

#sbatch script for running qiime2 on Wahab HPC at ODU
#set all variables in config.qiime file
#comment in or out modules below as needed

#SBATCH --job-name=qiime2      
#SBATCH -o qiime2-%j.log
#SBATCH --time=00:00:00

enable_lmod
module load container_env qiime2/2022.08


#this will trim the fq.gz files using the settings in config.4.all
#it is assumed that the directory you run this script from has the 
#fq.gz files

#Qiime2 modules available are:
#1: directory setup
#2: feature classifier setup
#3: file import
#4: DADA2 preprocessing
#5: taxonomic filtering 
#6: coverage/metadata filtration and coremetrics
#7: univariate statistics
#8: multivariate statistics
#9: ANCOM
#10: core features

#Qiime2 module to run

MODULE=8

#If you change config files alter this line
bash qiime2_HPC/qiime2.bash $MODULE config.qiime2_HPC_RIPBIO_C

