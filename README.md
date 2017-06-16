# RNA-seq
RNA-seq experiment analysis code

## Fastq analysis pipeline:

__`RNA-seq_slurm_job.sh`__

* Aligns `FASTQ` to specified reference genome
* Sorts and indexes output `BAM` file to use in genome browser
* Summarizes read counts per gene using `featureCounts`

#### Argument options:

* __EXPID__     Custom ID for output files.
* __RUNDIR__    Path to directory to run script and save output in.
* __FQ__        Absolute path to input fastq file.
* __GENNAME__   Reference genome file basename preceded by absolute path to directory
                containing the reference genome files. The directory must include:
     * `FASTA` file
     * matching `GFF` file.

Both files must use the same basename. If an existing Bowtie2 index with a matching basename
(Bowtie2's `bt2_base`) is found in the same directory it will be used;
otherwise a new index is built.
* __FEAT__      `GFF` feature type (featureCounts argument `-t`). Suggested values:
     * `FEAT="gene"` for SK1Yue
     * `FEAT="CDS"` for sacCer3

* __ATTR__      `GFF` attribute type used to group features (featureCounts argument `-g`). Suggested values:
     * `ATTR="ID"` for SK1Yue
     * `ATTR="Name"` for sacCer3

### Example job submission:

```
sbatch --export EXPID="AH119_3h_SK1Yue",RUNDIR="/scratch/lv38",\
FQ="/scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq.gz",\
GENNAME="/home/lv38/Library/SK1Yue/Yue.SK1.genome.nuclear.mito.2micr",\
FEAT="gene",ATTR="ID" ~/Pipeline/RNA-seq/RNA-seq_slurm_job.sh
```

## License
This project is licensed under the terms of the MIT license. See [LICENSE](LICENSE) file for details.
