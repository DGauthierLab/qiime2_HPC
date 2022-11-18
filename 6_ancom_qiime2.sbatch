#!/bin/bash -l

#This script is designed to perform post-filtration analyses 
#Sections are intended to be run sequentially, and should be commented out as needed

#SBATCH --job-name=qiime2_ancom
#SBATCH -o qiime2_ancom-%j.out
#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=32

enable_lmod
module load container_env qiime2

##SET metadata path, filtering conditions (must match filter_coremetrics script)

METAPATH=/cm/shared/courses/gauthier/Microbiome/metadata/metadata_19_21b.txt
FILTCOV=10000
FILTMETA=nocontrol_wild
FILTTAX=

##SET metadata column for beta comparisons

BETACOMPVAR=Rp_status

##SET diversity measure

DMALPHA=shannon
DMBETA=bray_curtis

##SET taxonomic level at which to collapse (optional)
#1-Kingdom, 2-Phylum, 3-Class, 4-Order, 5-Family, 6-Genus, 7-Species (no collapse)

COLV=5

##AUTOSET, Don't change
if [ -z $FILTCOV ]
then
echo "coverage filter not specified"
else
FILTCOV=_$FILTCOV
fi

if [ -z $FILTMETA ]
then
echo "metadata filter not specified"
else
FILTMETA=_$FILTMETA
fi

if [ -z $FILTTAX ]
then
echo "taxonomy filter not specified"
else
FILTTAX=_$FILTTAX
fi

COREMETRICS=core-metrics-results${FILTCOV}${FILTMETA}${FILTTAX}

##RUNSCRIPTS, don't change

crun qiime taxa collapse \
  --i-table table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
  --i-taxonomy taxonomy.qza \
  --p-level ${COLV} \
  --o-collapsed-table collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qza

##ANCOM analysis

crun qiime composition add-pseudocount \
--i-table table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
--o-composition-table comp-table${FILTCOV}${FILTMETA}${FILTTAX}.qza

crun qiime composition ancom \
--i-table comp-table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
--m-metadata-file $METAPATH \
--m-metadata-column $BETACOMPVAR \
--o-visualization ancom-$BETACOMPVAR-${FILTCOV}${FILTMETA}${FILTTAX}.qzv

#ANCOM analysis at specific taxonomic level

crun qiime taxa collapse \
--i-table table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
--i-taxonomy taxonomy.qza \
--p-level ${COLV} \
--o-collapsed-table collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qza 

crun qiime composition add-pseudocount \
--i-table collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
--o-composition-table comp-collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qza 

crun qiime composition ancom \
--i-table comp-collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qza \
--m-metadata-file $METAPATH \
--m-metadata-column $BETACOMPVAR \
--o-visualization comp-collapseLevel${COLV}-table${FILTCOV}${FILTMETA}${FILTTAX}.qzv
