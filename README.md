# RNA-seq
RNA-seq experiment analysis code

## Fastq analysis pipeline:

__`RNA-seq_slurm_job.sh`__

* Aligns `FASTQ` to specified reference genome
* Sorts and indexes output `BAM` file to use in genome browser
* Summarizes read counts per gene using `featureCounts`

#### Argument options:

* __EXPID__     Custom ID for output files
* __RUNDIR__    Path to directory to run script and save output in
* __FQ__        Absolute path to input fastq file
* __GENDIR__    Absolute path to directory containing reference genome files.
Must include:
      * `FASTA` file
      * matching `GFF` file.

If an existing Bowtie2 index with a basename (`bt2_base`) matching the `FASTA` file name is
found in the same directory it will be used; otherwise a new index is built
* __FEAT__      `GFF` feature type (featureCounts arg `-t`).
                Suggested values:
                
                * `FEAT="gene"` for SK1Yue
                * `FEAT="CDS"` for sacCer3 

* __ATTR__      `GFF` attribute type used to group features (featureCounts arg `-g`).
                Suggested values:
                
                * `ATTR="ID"` for SK1Yue
                * `ATTR="Name"` for sacCer3

### Example job submission:

```
sbatch --export EXPID="AH119_3h_SK1Yue",RUNDIR="/scratch/lv38",\
FQ="/scratch/lv38/C8C2NACXX_l03n01_ah119-3-030316.3510000004e291.fastq.gz",\
GENDIR="/home/lv38/Library/SK1Yue",FEAT="exon",ATTR="Name" \
~/Pipeline/RNA-seq/RNA-seq_slurm_job.sh
```

## License
This project is licensed under the terms of the MIT license. See [LICENSE](LICENSE) file for details.
