#!/bin/bash -l

#Preprocessing script for paired-end Illumina reads using dada2
#saves demux and joined files to folder called demux in top analysis directory.

#SBATCH --job-name=qiime2_import
#SBATCH -o logs/qiime2-%j.out
#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=32

enable_lmod
module load container_env qiime2

mkdir -p demux

##SET fastq folder path

FASTQ=/cm/shared/courses/gauthier/Microbiome/MM_microbiome/fastq



#.gz OK

#For demultiplexed MiSeq data, only change --input-path below

crun qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path FASTQ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux/demux-paired-end.qza

#Don't change
crun qiime demux summarize \
  --i-data demux/demux-paired-end.qza \
  --o-visualization demux/demux-paired-end.qzv

##Output from demux-paired-end.qzv should be examined to determine and set parameters for next step

#assembling pairs, assuming you have overlapping pairs
#this step not necessary if using dada2 for read processing

#crun qiime vsearch join-pairs \
#  --i-demultiplexed-seqs demux-paired-end.qza \
#  --o-joined-sequences demux-joined.qza \
#  --p-threads 4

#crun qiime demux summarize \
#  --i-data demux-joined.qza \
#  --o-visualization demux-joined.qzv

##Output from demux-joined.qzv should be examined to determine and set parameters for next step if not using dada2



#-------->3_proceed to preprocessing_qiime2.sbatch

