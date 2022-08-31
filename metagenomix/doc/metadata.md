# Metadata

The sample metadata is a mandatory [configuration file](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configuration-files).
It must be a simple tab-separated file with a header row containing at least 
one column.

The variable name for this first column can be anything, such as 
`sample_name`, `SampleID`, or any other character string such as 
`my_samples`.

This column must contain the sample names written **exactly** as they appear in 
the basename of the input fastq files that contain the sequence data for 
these samples.

For example, a user may provide these two folders as input for illumina and 
nanopore reads (using options `-i illumina_fastqs -k nanopore_fastqs`):

```
.
├── illumina_fastqs
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   └── sample2_2.fastq.gz
└── nanopore_fastqs
    ├── sample1.fastq.gz
    ├── sample2.fastq.gz
    └── sampleX.fastq.gz
```

Then, for these samples to be integrated in the analyses of the 
configured [pipeline](), it is crucial that the metadata looks like:

```
sample_name variable_A   variable_B
sample1 value_A1    value_B1
sample2 value_A2    value_B2
sample3 value_A2    value_B2
```

In this example, the long Nanopore reads in `sampleX` will be ignored and the 
metadata sample `sample3` will be ignored since there is no data for it.
