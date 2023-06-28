[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
![maintenance-status](https://img.shields.io/badge/maintenance-passively--maintained-yellowgreen.svg)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


# Clonal evolution reconstruction

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)

## Overview
Clonal evolution describes the accumulation of genetic changes within cancer cells as they undergo successive rounds of cell division. Initially, a single normal cell acquires a specific initiating mutation, known as a driver mutation, which provides a selective advantage and leads to the formation of a clonal population of cells. Over time, these cells accumulate additional genetic alterations, including both driver and passenger mutations, due to ongoing genomic instability and selective pressures within the tumor microenvironment.

Advancements in technology, such as whole exome sequencing (WES), now enable the reconstruction of a tumor's evolutionary history by analyzing mutational patterns across different regions of the tumor or longitudinally sampled specimens.

WES produces a list of mutations that requires functional annotation. To facilitate this process, a pipeline has been developed using nextflow. This pipeline builds upon the output of [variantalker](https://github.com/zhanyinx/variantalker) and allows for the reconstruction of clonal phylogeny. Clonal inference is performed using (https://github.com/Roth-Lab/pyclone-vi), while phylogeny reconstruction is achieved with [citup](https://github.com/amcpherson/citup). The resulting phylogeny can be visualized using [timescape](https://github.com/shahcompbio/timescape) through fishplot visualization.


## Installation
Clone the repo

```bash
git clone git@github.com:zhanyinx/clonal_evolution.git
```

Download the necessary databases (ASCAT databases for WES)

```bash
wget -r -N --no-parent -nH --cut-dirs=2 -P public_databases https://bioserver.ieo.it/repo/dima/ascat_wes_files/
```

Please update the "ascat_genome_basedir" parameter by replacing it with the absolute path of the downloaded databases.


## Usage


```bash
nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --outdir outdir
```

## Input

clonal_evolution takes as input a csv samplesheet with 9 columns (very similar to [nf-core/sarek pipeline](https://nf-co.re/sarek/3.2.3/usage))

__IMPORTANT: HEADER is required__ 

| patient        | sex | status | sample   | cram or bam | crai or bai  | cellularity  | maf       |
| -------------- | --- | ------ | -------- | ----------- | ------------ | ------------ | --------- |
| patient1       | XX  | 1      | sampleid | path2/cram  | path2/crai   | 0.1          | path2/maf |
| .....          | ... | ...... | ........ | ......      | .....        | .....        | ...       |



- patient: patient id

- sex: sex of the patient (XX or XY)

- status: tumor (1) or normal (0)

- sample: sample id

- cram (bam): cram (or bam) file from the sample

- crai (bai): index. of cram (bam)

- cellularity: if available, put the cellularity (number between 0 and 1) from histology. Alternatively leave it empty.

- maf: annotated vcf, output of variantalker

## Output


variantalker outputs for each sample multiple files organised in 3 folders


Output structure:

```
params.outdir

|-- patient
|   `-- pyclone
|       |-- patient.pyclone.tsv (input file for pyclone)
|       |-- patient.pyclone.output.tsv (cellular prevalence estimates from pyclone)
|   `-- citup
|       |-- patient.h5 (citup model)
|   `-- timescape
|       |-- patient.timescape.html (timescape result output)
```





