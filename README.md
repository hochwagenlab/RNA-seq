# RNA-seq
RNA-seq experiment analysis code

## Fastq analysis pipeline:

__`RNA-seq_slurm_job.sh`__

### Dependencies:

* Genome reference files (`FASTA` and matching `GFF`)
* `FASTQ` file

#### Argument options:

* __EXPID__     Custom ID for output files
* __RUNDIR__    Path to directory to run script and save output in
* __FQ__        Absolute path to input fastq file
* __GENDIR__    Absolute path to directory containing reference genome files.
                Must include:
                    * `FASTA` file
                    * matching `GFF` file
                An existing Bowtie2 index with a basename (`bt2_base`) matching
                the `FASTA` file name is used if found in the same directory;
                otherwise a new index is built

### Example job submission:

```
sbatch --export EXPID="AH119_3h",RUNDIR="/scratch/lv38",\
FQ="/scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq",\
GENDIR="/home/lv38/Library/SK1Yue" ~/Pipeline/RNA-seq_slurm_job.sh
```

## License
This project is licensed under the terms of the MIT license. See [LICENSE](LICENSE) file for details.
