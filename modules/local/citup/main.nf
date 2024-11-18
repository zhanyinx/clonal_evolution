process convert2citup{
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    fair true
    input:
        tuple val(patient), file(pyclone)
    output:
        tuple val(patient), file("*/citup_frequencies.tsv"), file("*/citup_clusters.tsv")
        file("*/mut_sample_clus.txt")
    script:
        """
        name="\$(basename $pyclone | sed 's/\\.pyclone\\.output\\.tsv//g')"
        pyclone2citup.R -p $pyclone -o \$name
        """
}


process citup{
    cpus 2
    maxRetries = 3
    memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    publishDir "${params.outdir}/${patient}/citup", mode: "copy"
    fair true
    input:
        tuple val(patient), file(frequency), file(clusters)
    output:
        tuple val(patient), file("*.h5")
    script:
        """
        run_citup_qip.py $frequency $clusters ${patient}.h5 --submit local
        """
}
