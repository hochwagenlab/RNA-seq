#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --job-name=RNAseq_analysis
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lv38@nyu.edu
#SBATCH --output=/scratch/%u/%x_%j.out
#SBATCH --error=/scratch/%u/%x_%j.err

#------------------------------------------------------------------------------#
#                                INSTRUCTIONS                                  #
#------------------------------------------------------------------------------#
# Dependencies:
#     - Genome reference: .fasta and .gff
#     - Bowtie2 index
#       (Assumed to be in '~/Library/Bowtie2_sacCer3')

### Argument options:
# EXPID     Custom ID for output files
# RUNDIR    Path to directory to run script and save output in
# FQ        Absolute path to input fastq file
# GENDIR    Absolute path to directory containing reference genome files.
#           Must include:
#               FASTA file
#               Matching GFF file
#           An existing Bowtie2 index with a basename ("bt2_base") matching the
#           FASTA file name is used if found in the same directory; otherwise a
#           new index is built

### EXAMPLE:
# sbatch --export EXPID="AH119_3h",RUNDIR="/scratch/lv38",\
# FQ="/scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq",\
# GENDIR="/home/lv38/Library/SK1Yue" ~/Pipeline/Slurm/RNA-seq-sacCer3.sh

#------------------------------------------------------------------------------#

echo ">>>>> Started pipeline:"
date

# Find reference genome files
FA=$(find $GENDIR -iname "*.fa*")
GFF=$(find $GENDIR -iname "*.gff")

# Exit if not found
if [ -z "$FA" ] || [ -z "$GFF" ]
then
    echo ">>>>> Could not find reference genome files"; exit 1;
fi

# Find Bowtie2 index (is there a file named as "fasta_base_name.1.bt2"?)
IX=$(basename $FA)                              # Drop path to file
IX=${IX%.*}                                     # Drop extension
CHECKIX=$(find $GENDIR -iname "${IX}.1.bt2")    # Search file
#IX=$(basename $IX | cut -d '.' -f 1)

# Build Bowtie2 index if not found
if [ -z "$CHECKIX" ]
then
    echo ">>>>> Building Bowtie2 index"
    module purge
    module load bowtie2/intel/2.2.9
    # Build index
    cd $GENDIR
    bowtie2-build -f $FA $IX
fi


#--------------------------------------------#
# Map reads to reference genome with TopHat2 #
#--------------------------------------------#
cd $RUNDIR

# Abort if output directory already exists
if [ -d "$EXPID" ]
then
    echo ">>>>> Directory already exists"; exit 1;
fi

mkdir ${EXPID}/
cd ${EXPID}/
mkdir ${EXPID}_sacCer3_TopHat2-nnjuncs/

echo ">>>>> Map reads with TopHat2..."
module purge
module load tophat/intel/2.1.1
module load bowtie2/intel/2.2.9
module load samtools/intel/1.3.1

tophat2 -p 8 \
    --no-novel-juncs \
    --library-type fr-firststrand \
    -o ${EXPID}_sacCer3_TopHat2-nnjuncs/ \
    -G $GFF \
    $GENDIR/$IX $FQ


#---------------------------------#
# Sort and index BAM file for IGV #
#---------------------------------#

echo ">>>>> Sort and index bam file for IGV..."
module purge
module load samtools/intel/1.3.1

samtools sort -o ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}-accepted_hits_s.bam \
    ${EXPID}_sacCer3_TopHat2-nnjuncs/accepted_hits.bam
samtools index ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}-accepted_hits_s.bam

#--------------------------------------------------#
# Count reads with featureCounts (Subread package) #
#--------------------------------------------------#

echo ">>>>> Count reads with featureCounts..."
module purge
module load subread/intel/1.5.1

featureCounts -s 2 \
    -t CDS \
    -g Name \
    -a $GFF \
    -o ${EXPID}_sacCer3_TopHat2-nnjuncs/${EXPID}_featureCounts.txt \
    ${EXPID}_sacCer3_TopHat2-nnjuncs/accepted_hits.bam


echo ">>>>> Completed pipeline:"
date
