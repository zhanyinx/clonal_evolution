process convert2citup{
    input:
        file(pyclone)
    output:
        tuple file("*/citup_frequencies.tsv"), file("*/citup_clusters.tsv")
        file("*/mut_sample_clus.txt")
    script:
        """
        name="\$(basename $pyclone | sed 's/\\.pyclone\\.output\\.tsv//g')"
        pyclone2citup.R -p $pyclone -o \$name
        """
}


process citup{
    cpus 10
    maxRetries = 2
    memory { 5.GB * task.attempt }
    input:
        tuple file(frequency), file(clusters)
    output:
        file("*.h5")
    script:
        """
        outname="\$(dirname $frequency)"
        # create output folder
        if ! [ -d ${launchDir}/${params.outdir}/${params.date}/\${outname}/citup/ ]; then
            mkdir -p ${launchDir}/${params.outdir}/${params.date}/\${outname}/citup/
        fi
        run_citup_qip.py $frequency $clusters \$outname.h5 --submit local
        cp \$outname.h5 ${launchDir}/${params.outdir}/${params.date}/\$outname/citup
        """
}