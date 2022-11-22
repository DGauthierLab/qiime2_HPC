#!/bin/bash -l

VERSION=1.1
#This script is a bash wrapper to perform various functions of qiime2 analysis 

#Prints out title and contact info
echo; echo -e "\n* qiime2_HPC v$VERSION by Dave Gauthier * \n"
echo -e "Contact dgauthie@odu.edu with any problems \n\n "

#with apologies to John Puritz and Chris Bird as this code is pilfered from dDocent

if [ -n "$1" ]; then
        #getting functions from command line
        FUNKTION=$(echo $1)
        echo; echo "Running script $FUNKTION..."
else
        echo ""; echo `date` "ERROR:    qiime2_HPC must be run with 2 arguments, "
        echo "                  qiime2.bash [function] [config file]"
        exit
fi


if [ -n $2 ]; then
        echo " "
        echo `date` "Files output to: " $(pwd)
        echo ""; echo `date` "Reading config file... "
        CONFIG=$2
        FASTQ=$(grep '.fastq files:' $CONFIG | awk '{print $3;}')
        METAPATH=$(grep 'metadata files:' $CONFIG | awk '{print $3;}')
	SILVAVERSION=$(grep 'qiime rescript get-silva-data --p-version <integer>' $CONFIG | awk '{print $1;}') 
	SILVATRGET=$(grep 'qiime rescript get-silva-data --p-target <string>' $CONFIG | awk '{print $1;}')
	FPRIMER=$(grep 'qiime feature-classifier extract-reads --p-f-primer' $CONFIG | awk '{print $1;}')
	RPRIMER=$(grep 'qiime feature-classifier extract-reads --p-r-primer' $CONFIG | awk '{print $1;}')
	MAXAMP=$(grep 'qiime feature-classifier extract-reads --p-max-length' $CONFIG | awk '{print $1;}')
	MINAMP=$(grep 'qiime feature-classifier extract-reads --p-min-length' $CONFIG | awk '{print $1;}')
       	FTRUNC=$(grep 'qiime dada2 denoise-paired --p-trunc-len-f <integer>' $CONFIG | awk '{print $1;}')
        RTRUNC=$(grep 'qiime dada2 denoise-paired --p-trunc-len-r <integer>' $CONFIG | awk '{print $1;}')
        FTRIM=$(grep 'qiime dada2 denoise-paired --trim-len-f <integer>' $CONFIG | awk '{print $1;}')
        RTRIM=$(grep 'qiime dada2 denoise-paired --trim-len-r <integer>' $CONFIG | awk '{print $1;}')
        CHIMERA=$(grep 'qiime dada2 denoise-paired --p-min-fold-parent-over-abundance <integer>' $CONFIG | awk '{print $1;}')
        ALPHARMIN=$(grep 'minimum sampling depth for alpha rarefaction' $CONFIG | awk '{print $1;}')
        ALPHARMAX=$(grep 'qiime diversity alpha-rarefaction --p-max-depth <integer>' $CONFIG | awk '{print $1;}')
        ALPHARSTEPS=$(grep 'qiime diversity alpha-rarefaction --p-steps <integer>' $CONFIG | awk '{print $1;}')
        FILTCOV=$(grep 'qiime feature-table filter-samples --p-min-frequency <integer>' $CONFIG | awk '{print $1;}')
        FILTTAX=$(grep 'FILTTAX filename variable <string> ' $CONFIG | awk '{print $1;}')
        FILTMETA=$(grep 'FILTMETA filename variable <string>' $CONFIG | awk '{print $1;}')
        TFILT1=$(grep 'taxonomy filtration argument 1 <exclude/include> <"string">' $CONFIG | awk '{print $1, $2;}')
        TMODE1=$(grep 'taxonomy filtration argument 1 <exact/contains> <"string">' $CONFIG | awk '{print $1, $2;}')
        TFILT2=$(grep 'taxonomy filtration argument 2 <exclude/include> <"string">' $CONFIG | awk '{print $1, $2;}')
        TMODE2=$(grep 'taxonomy filtration argument 2 <exact/contains> <"string">' $CONFIG | awk '{print $1, $2;}')
        TFILT3=$(grep 'taxonomy filtration argument 3 <exclude/include> <"string">' $CONFIG | awk '{print $1, $2;}')
        TMODE3=$(grep 'taxonomy filtration argument 3 <exact/contains> <"string">' $CONFIG | awk '{print $1, $2;}')
        MFILT=$(grep 'metadata filtration argument 1  <mySQL where statement>' $CONFIG | gawk '{if (match($0,/where.*"/,m)) print m[0]}')
        DMALPHA=$(grep 'alpha diversity metric <string>' $CONFIG | awk '{print $1;}')
        DMBETA=$(grep 'beta diversity metric <string>' $CONFIG | awk '{print $1;}')
        BETACOMPVAR=$(grep 'comparison variable for univariate beta diversity comparison <string> ' $CONFIG | awk '{print $1;}')
        PERMFORM=$(grep 'qiime diversity adonis --p-formula <string>  ' $CONFIG | awk '{print $1;}')
        COLV=$(grep 'qiime taxa collapse --p-level <integer>' $CONFIG | awk '{print $1;}')
fi

echo "Done reading config file.  Proceeding to module $1"

if [[ $FUNKTION == 1 ]];
	then
		mkdir -p metadata
		mkdir -p fastq
		echo "place tab-delimited .txt metadata file in metadata folder"
		echo "place all .fastq files (.gz is OK) in fastq folder.  Ensure that only .fastq files are present in the folder"

elif [[ $FUNKTION == 2 ]];
	then
		echo "Running 2_feature_classifier_qiime2"
		echo "Silva version is $SILVAVERSION"
		echo "Silva target is $SILVATARGET"
		echo "Forward sequencing primer: $FPRIMER"
		echo "Reverse sequencing primer: $RPRIMER"
		echo "Min expected amplicon length: $MINAMP"
		echo "Max expected amplicon length: $MAXAMP"

#This script will train the feature classifier for use in taxonomy analysis.  
#see: https://docs.qiime2.org/2021.4/tutorials/feature-classifier/ for details
 
mkdir training_feature_classifiers

crun qiime rescript get-silva-data \
--p-version '138.1' \
--p-target 'SSURef_NR99' \
--p-include-species-labels \
--p-download-sequences \
--o-silva-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs.qza \
--o-silva-taxonomy training_feature_classifiers/silva-138.1-ssu-nr99-tax.qza

crun qiime rescript cull-seqs \
--i-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs.qza \
--o-clean-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned.qza

crun qiime rescript filter-seqs-length-by-taxon \
--i-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned.qza \
--i-taxonomy training_feature_classifiers/silva-138.1-ssu-nr99-tax.qza \
--p-labels Archaea Bacteria Eukaryota \
--p-min-lens 900 1200 1400 \
--o-filtered-seqs training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned_filt.qza \
--o-discarded-seqs training_feature_classifiers/silva-138-ssu-nr99-seqs_discard.qza

crun qiime rescript dereplicate \
--i-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned_filt.qza \
--i-taxa training_feature_classifiers/silva-138.1-ssu-nr99-tax.qza \
--p-rank-handles 'silva' \
--p-mode 'uniq' \
--o-dereplicated-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned_filt_derep.qza \
--o-dereplicated-taxa training_feature_classifiers/silva-138.1-ssu-nr99-tax_derep.qza

#Seems not to be necessary if previous steps run.  If no cleaning was done, this is necessary.

#crun qiime rescript reverse-transcribe \
#--i-rna-sequences silva-138.1-ssu-nr99-seqs_cleaned_filt_derep.qza \
#--o-dna-sequences silva-138.1-ssu-nr99-seqs_cleaned_filt_derep_dna.qza

crun qiime feature-classifier extract-reads \
--i-sequences training_feature_classifiers/silva-138.1-ssu-nr99-seqs_cleaned_filt_derep.qza \
--p-f-primer CCTACGGGNGGCWGCAG \
--p-r-primer GACTACHVGGGTATCTAATCC \
--p-min-length 350 \
--p-max-length 550 \
--o-reads training_feature_classifiers/ref-seqs_silva138_NR99.qza

crun qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads training_feature_classifiers/ref-seqs_silva138_NR99.qza \
--i-reference-taxonomy training_feature_classifiers/silva-138.1-ssu-nr99-tax_derep.qza \
--o-classifier training_feature_classifiers/slv_ssu_138.1_classifier.qza

	echo "Module 2 completed successfully"

elif [[ $FUNKTION == 3 ]];
        then
                echo "Running 3_file_import_qiime2" 
		echo "Finding .fastq files at $FASTQ"

#Preprocessing script for paired-end Illumina reads using dada2
#saves demux and joined files to folder called demux in top analysis directory.

mkdir -p demux

#For demultiplexed MiSeq data, only change --input-path below

crun qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${FASTQ} \
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

	echo "Module 3 completed successfully"

elif [[ $FUNKTION == 4 ]];
        then
                echo "Running 4_DADA2_preprocessing_qiime2" 
		echo "Forward truncation at $FTRUNC bp"
		echo "Reverse truncation at $RTRUNC bp"
		echo "Forward trimming $FTRIM bp"
		echo "Reverse trimming $RTRIM bp"
		echo "Chimera fold overabundance: $CHIMERA"
		echo -e "\nAlpha rarefaction performed from $ALPHARMIN to $ALPHARMAX X coverage in $ALPHARSTEPS X steps"

#Preprocessing script for paired-end Illumina reads using dada2
#requires demux paired end files in folder called demux under top analysis directory

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
  --i-classifier ../training_feature_classifiers/slv_ssu_138.1_classifier.qza \
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
  --p-steps $ALPHARSTEPS \
  --p-min-depth $ALPHARMIN \
  --p-max-depth $ALPHARMAX \
  --m-metadata-file $METAPATH \
  --o-visualization alpha-rarefaction.qzv


elif [[ $FUNKTION == 5 ]];
        then
                echo "Running 5_filtering_and_coremetrics_qiime2" 
		echo "Minimum coverage set at $FILTCOV"
		echo "Taxonomic filtration is $FILTTAX"
		echo "Metadata filtration is $FILTMETA"


if [ -z $FILTTAX ]
        then
                echo "taxonomic filter not specified"
        else
                FILTTAX=_$FILTTAX

                if [ -z $TFILT1 ]
                        then
                                echo "taxonomic filter 1 not specified"
                        else
                                if [ -z $TMODE1 
                                        then
                                                echo "taxonomic mode 1 not specified"

                                        else
                                                crun qiime taxa filter-table \
                                                --i-table table.qza \
                                                --i-taxonomy taxonomy.qza \
                                                --p-${TMODE1} \
                                                --p-${TFILT1} \
                                                --o-filtered-table table${FILTTAX}.qza 

                if [ -z $TFILT2 ]

                        then
                                echo "taxonomic filter 2 not specified"
                        else
                                if [ -z $TMODE2 ]
                                        then
                                                echo "taxonomic mode 2 not specified"

                                        else

                                                crun qiime taxa filter-table \
                                                --i-table table.qza \
                                                --i-taxonomy taxonomy.qza \
                                                --p-${TMODE2} \
                                                --p-${TFILT2} \
                                                --o-filtered-table table${FILTTAX}.qza

                if [ -z $TFILT3 ]

                        then
                                echo "taxonomic filter 3 not specified"

                        else
                                if [ -z $TMODE3 ]
                                        then
                                                echo "taxonomic mode 3 not specified"

                                        else
                                                crun qiime taxa filter-table \
                                                --i-table table${FILTTAX}.qza \
                                                --i-taxonomy taxonomy.qza \
                                                --p-${TMODE3} \
                                                --p-${TMODE3} \
                                                --o-filtered-table table${FILTTAX}.qza 

fi

if [ -z $FILTCOV ]
        then
                echo "coverage filter not specified"
        else
                COV=$FILTCOV
                FILTCOV=_$FILTCOV
                        crun qiime feature-table filter-samples \
                          --i-table table${FILTTAX}.qza \
                          --p-min-frequency $COV \
                          --o-filtered-table table${FILTTAX}${FILTCOV}.qza
fi

if [ -z $FILTMETA ]
        then
                echo "metadata filter not specified"
        else
                FILTMETA=_$FILTMETA
                crun qiime feature-table filter-samples \
                  --i-table table${FILTTAX}${FILTCOV}.qza \
                  --m-metadata-file $METAPATH \
                  --p-${MFILT} \
                  --o-filtered-table table${FILTTAX}${FILTCOV}${FILTMETA}.qza
fi

##Summarize the new filtered table
crun qiime feature-table summarize \
  --i-table table${FILTTAX}${FILTCOV}${FILTMETA}.qza \
  --o-visualization table${FILTTAX}${FILTCOV}${FILTMETA}.qzv \
  --m-sample-metadata-file $METAPATH

##Alpha and Beta diversity metrics
crun qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table${FILTTAX}${FILTCOV}${FILTMETA}.qza \
  --p-sampling-depth $COV \
  --m-metadata-file $METAPATH \
  --output-dir core-metrics-results${FILTTAX}${FILTCOV}${FILTMETA}

##change filenames to reflect filtering parameters
cd core-metrics-results${FILTTAX}${FILTCOV}${FILTMETA} &&
for f in * ;
do mv "$f" $(echo "$f" | sed "s/\(.*\).qz\([av]\)/\1${FILTTAX}${FILTCOV}${FILTMETA}.qz\2/g") ;
done
cd ..

##Following assumes taxonomic classification has already been run
crun qiime taxa barplot \
--i-table table${FILTTAX}${FILTCOV}${FILTMETA}.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file $METAPATH \
--o-visualization barplot${FILTTAX}${FILTCOV}${FILTMETA}.qzv

	echo "Module 5 completed successfully"

elif [[ $FUNKTION == 6 ]];
        then
                echo "Running 6_univariate_stats_qiime2" 

elif [[ $FUNKTION == 7 ]];
        then
                echo "Running 7_multivariate_stats_qiime2" 

elif [[ $FUNKTION == 8 ]];
        then
                echo "Running 8_ANCOM_feature_classifier_qiime2" 

fi

