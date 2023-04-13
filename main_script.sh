#!/bin/bash

#SBATCH --job-name=protOrth        # job name
#SBATCH --nodes=1                            # nodes
#SBATCH -p long                              # queue
#SBATCH -A jgu-cbdm
#SBATCH -c 20       # nb of cores/task
#SBATCH --mem=34000M                          # memory
#SBATCH --time=2-00:00:00                      # time
#SBATCH --ramdisk=30G   # ramdisk
##SBATCH --error=tmp.err                 # error file name
##SBATCH --output=tmp.out                # output file name


# run with: sbatch launch_Rscript_new.sh filename


# filenames are: 
# 

# recommandation: sbatch launch_Rscript.sh 20 # 20 is nb of cores to use / optimal value 


# Modules loading

#module load lang/R

module load lang/R/3.4.3-intel-2018.02-X11-20171023

#module load bio/BLAST+/2.6.0-foss-2017a-Python-2.7.13  # proper versions to avoid gcc compiler conflict
#module load lang/R/3.4.1-foss-2017a

# Global variables

JOBDIR=/localscratch/$SLURM_JOB_ID 
RAMDISK=$JOBDIR/ramdisk # folder associated to the ramdisk
#cores=$1  # nb of cores to use in parallelization
cores=$1  # nb of cores to use per R script

# Main script

# step1: retrieving orthology 1:1 information from biomaRt

srun Rscript ./biomaRt.R

# step2: filling putative missing data by doing reciprocal searches

# transfering the sequence databases to the folder associated with the ramdisk:

#cp -r /home/dfournie/cbdm/prot_evol/orthologues_finding/blastp/db $RAMDISK
cp -r /home/dfournie/dfournie/prot_evol/orthologues_finding/blastp/db $RAMDISK # temporary working directory til fs6 back to normal

cd /home/dfournie/dfournie/prot_evol/orthologues_finding

mkdir tmp1; mkdir tmp2
tail -n +2 orthos_1to1_absence_presence.txt | split -500   # spits file into chunks of 500 bits - also removes the header
mv x* tmp1/
awk '{print $6}' orthos_1to1_metadata.txt | paste  orthos_1to1_gene_names.txt - > test.txt
tail -n +2 test.txt | split -500  
rm test.txt
mv x* tmp2/

# launch of the main R script:

> Rreport.txt # report of R steps
for file in tmp2/x*; do
filename="${file##*/}"  # extracts file name from full path
#for file in x*; do 
echo "processing file ${filename}"
sleep 2
srun Rscript ./blast.R ${RAMDISK} $cores $filename  # passing ramdisk directory to R script
done

# cleaning - better do it manually

# cd ..
#rm -r tmp1; rm -r tmp2
