#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --job-name=RNAseq_sacCer3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lv38@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out
#SBATCH --error=/scratch/%u/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#
# Dependencies:
#     - S288C sacCer3 genome reference: .fasta and .gff
#     - Bowtie2 index
#       (Assumed to be in '~/Library/Bowtie2_sacCer3')

### RUN AS:
# sbatch --export EXPID="for_output_files",RUNDIR="path/to/dir",\
# FQ="absolute/path/to/reads.fastq",\
# ~/Pipeline/Slurm/RNA-seq-sacCer3.sh

### EXAMPLE:
# sbatch --export EXPID="AH119_3h",RUNDIR="scratch/lv38",\
# T_MAP="scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq" \
# ~/Pipeline/Slurm/RNA-seq-sacCer3.sh

#------------------------------------------------------------------------------#

date
cd ${RUNDIR}

# Abort if output directory already exists
if [ -d "${EXPID}" ]
then
    echo "Directory already exists"; exit 1;
fi

mkdir ${EXPID}/
cd ${EXPID}/

#--------------------------------------------#
# Map reads to reference genome with TopHat2 #
#--------------------------------------------#

mkdir ${EXPID}_sacCer3_TopHat2-nnjuncs/

echo "Map reads with TopHat2..."
module load tophat/intel/2.1.1
module load bowtie2/intel/2.2.9
module load samtools/intel/1.3.1

tophat2 -p 8 \
    --no-novel-juncs \
    --library-type fr-firststrand \
    -o ${EXPID}_sacCer3_TopHat2-nnjuncs/ \
    -G ~/Library/Bowtie2_sacCer3/S288C_chr_copy.gff \
    ~/Library/Bowtie2_sacCer3/sacCer3 ${FQ}


#---------------------------------#
# Sort and index BAM file for IGV #
#---------------------------------#

echo "Sort and index bam file for IGV..."
module purge
module load samtools/intel/1.3.1

samtools sort -o ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}-accepted_hits_s.bam \
    ${EXPID}_sacCer3_TopHat2-nnjuncs/accepted_hits.bam
samtools index ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}-accepted_hits_s.bam

#--------------------------------------------------#
# Count reads with featureCounts (Subread package) #
#--------------------------------------------------#

echo "Count reads with featureCounts..."
module purge
module load subread/intel/1.5.1

featureCounts -s 2 \
    -t CDS \
    -g Name \
    -a ~/Library/Bowtie2_sacCer3/S288C_chr_copy.gff \
    -o ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}_featureCounts.txt \
    ${EXPID}_sacCer3_TopHat2-nnjuncs/accepted_hits.bam


echo "Done!"
date
