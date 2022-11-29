# wahab_qiime2

being scripts and files for running basic qiime2 functions on Wahab HPC cluster at ODU

Three critical files present:
	
qiime2_HPC.sbatch
	sbatch script for running qiime modules on cluster.  Only module number needs to be changed according to the legend given therein.  Specific functions of each modules subject of next README revision.

qiime2_HPC.bash
	bash wrapper script containing qiime commands and variable setting functions.  Ask DG for modifications/bugfixes.

config.qiime2_HPC
	configuration script for variables desired during qiime2 analysis.  Recommend setting up entire config file before running anything.  Variables can be changed and filenames will change accordingly.


Getting started:

1) Create a uniquely named folder for your top directory.  Clone or move the qiime2_HPC folder therein.  
2) All commands will be issued from within the top directory.  Pay attention to filepaths.
3) %nano qiime2_HPC/qiime2_HPC.sbatch
4) change module number to desired module and save/exit.  Start with Module 1 
5) %nano qiime2_HPC/config.qiime2_HPC
6) adjust variables according to taste.  Save and close. 
7) %sbatch qiime2_HPC/qiime2_HPC.sbatch
8) logfiles will appear in top directory as qiime2-xxxxxxx.out.  Consult them for information and possible run errors.
9) module 1 will create fastq and metadata directories in the top directory.  Transfer all fastq files to the fastq directory and metadata file to metadata directory.  If you wish to use a different location for these files, modify the config.qiime2_HPC file accordingly.
10) proceed with modules


	
