#!/bin/bash -l

#Preprocessing script for paired-end Illumina reads using dada2
#requires demux paired end files in folder called demux under top analysis directory

#SBATCH --job-name=qiime2_dada2_preprocessing
#SBATCH -o logs/qiime2-%j.out
#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=32

enable_lmod
module load container_env qiime2

##SET metadata path

METAPATH=metadata/metadata_19_21b.txt

##SET truncation of F and R reads

FTRUNC=275
RTRUNC=275

##SET left trims of F and R reads

FTRIM=17
RTRIM=21

##SET fold-overabundance for chimera detection

CHIMERA=2

mkdir -p dada2_${FTRIM}_${RTRIM}_${FTRUNC}_${RTRUNC}_p${CHIMERA}
cd dada2_${FTRIM}_${RTRIM}_${FTRUNC}_${RTRUNC}_p${CHIMERA}


##follows from 1_file_import_qiime2 script
#examine demux-paired-end.qzv file to estimate settings
#note min-fold-parent-over-abundance setting here (reduces # of detected chimeras)
crun qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ../demux/demux-paired-end.qza \
  --p-trunc-len-f ${FTRUNC} \
  --p-trunc-len-r ${RTRUNC} \
  --p-trim-left-f ${FTRIM} \
  --p-trim-left-r ${RTRIM} \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza \
  --p-min-fold-parent-over-abundance ${CHIMERA} \
  --p-n-threads 0

#Shouldn't need changed
crun qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

#Shouldn't need changed
crun qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

##Taxonomic analysis (generate taxonomy.qza file)

crun qiime feature-classifier classify-sklearn \
  --i-classifier /cm/shared/courses/gauthier/Microbiome/training_feature_classifiers/slv_ssu_138.1_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

crun qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

##Run alpha-rarefaction analysis on whole data
##Alpha Rarefaction

crun qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-steps 50 \
  --p-min-depth 1000 \
  --p-max-depth 50000 \
  --m-metadata-file $METAPATH \
  --o-visualization alpha-rarefaction.qzv

