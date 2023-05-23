#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowMain.initialise(workflow, params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.input, params.fasta, 
    params.target, 
    ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLONAL_EVOLUTION } from './workflows/clonal_evolution'

process help {
    label "help"

    script:
    """
    echo "Usage: nextflow run main.nf -c nextflow.config"
    echo ""
    echo "This pipeline performs whole exome sequencing clinical annotation (snp, indel, cnv)."
    echo ""
    echo "Options:"
    echo "  --input samplesheet.csv                      CSV file with 4 columns: patient, tumor_tissue, sample_file, sample_type." 
    echo "                                               Available sample_type: somatic, germline and cnv. For somatic and germline" 
    echo "                                               vcf.gz file is accepted for sample_file. For CNV, either cnr from CNVkit or" 
    echo "                                               vcf.gz file from dragen is accepted." 
    echo "                                               For available tumor type, checkout https://github.com/zhanyinx/variantalker/blob/main/README.md" 
    echo "  --pipeline <pipeline>                        The pipeline used to generate the input data (options: 'Sarek', 'DRAGEN'"
    echo "                                               or 'Iontorrent' (only SNP/INDEL), default: 'Sarek')." 
    echo ""
    """
}

workflow {
   CLONAL_EVOLUTION()
}