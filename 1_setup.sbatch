#!/bin/bash -l

#Folder setup for QIIME2 processing on Wahab HPC cluster

#SBATCH --job-name=folder_setup
#SBATCH -o qiime2-%j.out
#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=32

mkdir -p metadata
mkdir -p fastq

##Subsequent to running this file:
#place tab-delimited .txt metadata file in metadata folder
#place all .fastq files (.gz is OK) in fastq folder
