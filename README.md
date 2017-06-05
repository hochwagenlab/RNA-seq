# RNA-seq
RNA-seq experiment analysis code

## Fastq analysis pipeline: `RNA-seq_sacCer3_slurm_job.sh`

### Dependencies:

* S288C sacCer3 genome reference: .fasta and .gff
* Bowtie2 index

(must be in '~/Library/Bowtie2_sacCer3')

### Run as:
<pre>
sbatch --export EXPID="for_output_files",RUNDIR="path/to/dir",\
 FQ="absolute/path/to/reads.fastq",\
 ~/Pipeline/Slurm/RNA-seq-sacCer3.sh
</pre>

### Example:

<pre>
sbatch --export EXPID="AH119_3h",RUNDIR="scratch/lv38",\
 T_MAP="scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq" \
 ~/Pipeline/Slurm/RNA-seq-sacCer3.sh
</pre>

## License
This project is licensed under the terms of the MIT license. See [LICENSE](LICENSE) file for details.
