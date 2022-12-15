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

#Qiime2 modules available are:
#1: directory setup
#2: feature classifier setup
#3: file import
#4: DADA2 preprocessing
#5: filtering and coremetrics
#5a: mass-filtering of multiple taxa
#6: univariate statistics
#7: multivariate statistics
#8: ANCOM


#Qiime2 module to run
MODULE=51

bash qiime2_HPC/qiime2.bash $MODULE qiime2_HPC/config.qiime2_HPC


