
process ASCAT {
    cpus '4'
    errorStrategy 'retry'
    maxRetries = 2
    memory { 40.GB * task.attempt }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0':
        'quay.io/biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:dfe5aaa885de434adb2b490b68972c5840c6d761-0' }"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path(allele_files)
    path(loci_files)
    path(bed_file)  // optional
    path(fasta)     // optional
    path(gc_file)   // optional
    path(rt_file)   // optional

    output:
    tuple val(meta.patient), val(meta.tumor_id), val(meta.sex), val(meta.maf), val(meta.cellularity), path("${input_tumor}"), path("${index_tumor}"), path("*purityploidy.txt"),  path("*cnvs.txt"),          emit: cram
    path("*.purity_fit.tsv"), emit: purity_fit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args        ?: ''
    def prefix         = task.ext.prefix      ?: "${meta.id}"
    def gender         = args.gender          ?  "$args.gender" :        "NULL"
    def genomeVersion  = args.genomeVersion   ?  "$args.genomeVersion" : "NULL"
    def purity         = args.purity          ?  "$args.purity" :        "NULL"
    def ploidy         = args.ploidy          ?  "$args.ploidy" :        "NULL"
    def gc_input       = gc_file              ?  "$gc_file" :            "NULL"
    def rt_input       = rt_file              ?  "$rt_file" :            "NULL"

    def minCounts_arg                    = args.minCounts                     ?  ",minCounts = $args.minCounts" : ""
    def bed_file_arg                     = bed_file                           ?  ",BED_file = '$bed_file'": ""
    def chrom_names_arg                  = args.chrom_names                   ?  ",chrom_names = $args.chrom_names" : ""
    def min_base_qual_arg                = args.min_base_qual                 ?  ",min_base_qual = $args.min_base_qual" : ""
    def min_map_qual_arg                 = args.min_map_qual                  ?  ",min_map_qual = $args.min_map_qual" : ""
    def fasta_arg                        = fasta                              ?  ",ref.fasta = '$fasta'" : ""
    def skip_allele_counting_tumour_arg  = args.skip_allele_counting_tumour   ?  ",skip_allele_counting_tumour = $args.skip_allele_counting_tumour" : ""
    def skip_allele_counting_normal_arg  = args.skip_allele_counting_normal   ?  ",skip_allele_counting_normal = $args.skip_allele_counting_normal" : ""

    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    #build prefixes: <abspath_to_files/prefix_chr>
    allele_path = normalizePath("$allele_files")
    allele_prefix = paste0(allele_path, "/", "$allele_files", "_chr")

    loci_path = normalizePath("$loci_files")
    loci_prefix = paste0(loci_path, "/", "$loci_files", "_chr")

    #prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "$input_tumor",
        normalseqfile = "$input_normal",
        tumourname = paste0("$prefix", ".tumour"),
        normalname = paste0("$prefix", ".normal"),
        allelecounter_exe = "alleleCounter",
        alleles.prefix = allele_prefix,
        loci.prefix = loci_prefix,
        gender = "$gender",
        genomeVersion = "$genomeVersion",
        nthreads = $task.cpus
        $minCounts_arg
        $bed_file_arg
        $chrom_names_arg
        $min_base_qual_arg
        $min_map_qual_arg
        $fasta_arg
        $skip_allele_counting_tumour_arg
        $skip_allele_counting_normal_arg
    )


    #Load the data
    ascat.bc = ascat.loadData(
        Tumor_LogR_file = paste0("$prefix", ".tumour_tumourLogR.txt"),
        Tumor_BAF_file = paste0("$prefix", ".tumour_tumourBAF.txt"),
        Germline_LogR_file = paste0("$prefix", ".tumour_normalLogR.txt"),
        Germline_BAF_file = paste0("$prefix", ".tumour_normalBAF.txt"),
        genomeVersion = "$genomeVersion",
        gender = "$gender"
    )

    #Plot the raw data
    ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".before_correction."))

    # optional LogRCorrection
    if("$gc_input" != "NULL") {
        gc_input = paste0(normalizePath("$gc_input"), "/", "$gc_input", ".txt")

        if("$rt_input" != "NULL"){
            rt_input = paste0(normalizePath("$rt_input"), "/", "$rt_input", ".txt")
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = rt_input)
            #Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc_rt."))
        }
        else {
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = $rt_input)
            #Plot raw data after correction
            ascat.plotRawData(ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc."))
        }
    }

    #Segment the data
    ascat.bc = ascat.aspcf(ascat.bc)

    #Plot the segmented data
    ascat.plotSegmentedData(ascat.bc)

    #Fit. Purity (rho) is forced: the samplesheet cellularity when given, else the VAF
    #purity of the unconstrained fit (see vaf_purity). Ploidy (psi) is then the best
    #goodness of fit over a grid at that fixed rho. A failed fit stops the run.
    fit_free = tryCatch(ascat.runAscat(ascat.bc, gamma=1, img.dir=tempdir()), error=function(e) NULL)
    free_ok = isTRUE(!is.null(fit_free) && length(fit_free[["rho"]]) > 0 && !is.na(fit_free[["rho"]][1]))

    vaf_purity = function(segs, maf_path) {
        #2 x median VAF of the top VAF cluster (above the largest gap) of autosomal MAF SNVs
        #in ASCAT 1/1 segments, depth >= 30, alt >= 5, VAF <= 0.6; top cluster >= 3 SNVs.
        lines = readLines(maf_path)
        maf = read.delim(text=lines[!startsWith(lines, "#")], stringsAsFactors=FALSE, quote="", check.names=FALSE)
        snv = maf[maf[["Variant_Type"]] == "SNP", ]
        chr = sub("^chr", "", snv[["Chromosome"]]); pos = as.numeric(snv[["Start_Position"]])
        alt = as.numeric(snv[["t_alt_count"]]); dp = alt + as.numeric(snv[["t_ref_count"]])
        if (is.null(segs)) {
            in_dip = rep(TRUE, length(chr))  #no fit: every autosomal SNV is taken as copy neutral
        } else {
            segs[["chr"]] = sub("^chr", "", segs[["chr"]])
            dip = segs[segs[["nMajor"]] == 1 & segs[["nMinor"]] == 1, ]
            in_dip = mapply(function(c, p) any(dip[["chr"]] == c & dip[["startpos"]] <= p & dip[["endpos"]] >= p), chr, pos)
        }
        ok = chr %in% as.character(1:22) & dp >= 30 & alt >= 5 & in_dip
        v = sort((alt / dp)[ok]); v = v[v <= 0.6]
        n = length(v)
        if (n < 3) stop(paste("VAF purity: only", n, "heterozygous SNVs in diploid segments; set cellularity in the samplesheet"))
        top = v
        if (n > 3) {
            g = diff(v); cand = which(n - seq_along(g) >= 3)
            i = cand[which.max(g[cand])]; top = v[(i + 1):n]
        }
        list(purity = min(1, 2 * median(top)), n = n, n_top = length(top))
    }

    purity_given = "${purity}"
    if (purity_given != "NULL") {
        rho = as.numeric(purity_given); purity_source = "samplesheet"; n_snv = NA; n_top = NA
    } else {
        vp = vaf_purity(if (free_ok) fit_free[["segments"]] else NULL, "${meta.maf}")
        rho = vp[["purity"]]; n_snv = vp[["n"]]; n_top = vp[["n_top"]]
        purity_source = if (free_ok) "vaf" else "vaf_no_ascat_fit"
    }
    if (is.na(rho) || rho <= 0 || rho > 1) stop(paste("Purity", rho, "is not in (0, 1]; stopping"))

    best = NULL; best_gof = -Inf; best_psi = NA
    psi_grid = if ("${ploidy}" != "NULL") as.numeric("${ploidy}") else seq(1.5, 5.5, by=0.05)
    for (ps in psi_grid) {
        f = tryCatch(ascat.runAscat(ascat.bc, gamma=1, rho_manual=rho, psi_manual=ps, img.dir=tempdir()), error=function(e) NULL)
        if (!is.null(f) && !is.na(f[["goodnessOfFit"]][1]) && f[["goodnessOfFit"]][1] > best_gof) {
            best = f; best_gof = f[["goodnessOfFit"]][1]; best_psi = ps
        }
    }
    if (is.null(best)) stop(paste("ASCAT fit at forced purity", rho, "failed for every ploidy"))
    #re-run the chosen fit to write its plots to the work dir
    ascat.output = ascat.runAscat(ascat.bc, gamma=1, rho_manual=rho, psi_manual=best_psi)

    fit = data.frame(sample="${meta.tumor_id}", purity=rho, purity_source=purity_source,
                     vaf_purity_snvs=n_snv, vaf_purity_top=n_top, psi=best_psi, goodness_of_fit=best_gof,
                     ploidy=ascat.output[["ploidy"]][1],
                     free_purity=if (free_ok) fit_free[["aberrantcellfraction"]][1] else NA,
                     free_ploidy=if (free_ok) fit_free[["ploidy"]][1] else NA,
                     free_goodness_of_fit=if (free_ok) fit_free[["goodnessOfFit"]][1] else NA)
    write.table(fit, file=paste0("$prefix", ".purity_fit.tsv"), sep="\t", quote=F, row.names=F)

    #Extract metrics from ASCAT profiles
    QC = ascat.metrics(ascat.bc,ascat.output)

    #Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("$prefix", ".segments.txt"), sep="\t", quote=F, row.names=F)

    #Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("$prefix",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #Write out purity and ploidy info (the forced purity)
    summary = matrix(c(rho, ascat.output[["ploidy"]][1]), ncol=2, byrow=TRUE)
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("$prefix",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    write.table(QC, file=paste0("$prefix", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # version export
    f <- file("versions.yml","w")
    alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
    ascat_version = sessionInfo()\$otherPkgs\$ASCAT\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    alleleCounter:", alleleCounter_version), f)
    writeLines(paste("    ascat:", ascat_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    echo stub > ${prefix}.before_correction.test.tumour.germline.png
    echo stub > ${prefix}.before_correction.test.tumour.tumour.png
    echo stub > ${prefix}.cnvs.txt
    echo stub > ${prefix}.metrics.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.purityploidy.txt
    echo stub > ${prefix}.purity_fit.tsv
    echo stub > ${prefix}.segments.txt
    echo stub > ${prefix}.tumour.ASPCF.png
    echo stub > ${prefix}.tumour.sunrise.png
    echo stub > ${prefix}.tumour_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.tumour_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.tumour_normalBAF.txt
    echo stub > ${prefix}.tumour_normalLogR.txt
    echo stub > ${prefix}.tumour_tumourBAF.txt
    echo stub > ${prefix}.tumour_tumourLogR.txt

    echo "${task.process}:" > versions.yml
    echo ' alleleCounter: 4.3.0' >> versions.yml
    echo ' ascat: 3.0.0' >> versions.yml

    """


}

process generate_ascat_loci {
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    label 'process_single'
    input:
        file(loci)
    output:
        path("${loci.baseName}")
    script:
    """
        mkdir ${loci.baseName}
        cd ${loci.baseName}
        unzip ../$loci
        cd ..
    """
}

process generate_ascat_alleles {
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    label 'process_single'
    input:
        file(alleles)
    output:
        path("${alleles.baseName}")
    script:
    """
        mkdir ${alleles.baseName}
        cd ${alleles.baseName}
        unzip ../$alleles
        cd ..
    """
}


process generate_ascat_rt {
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    label 'process_single'
    input:
        file(rt)
    output:
        path("${rt.baseName}")
    script:
    """
        mkdir ${rt.baseName}
        cd ${rt.baseName}
        unzip ../$rt
        cd ..
    """
}


process generate_ascat_gc {
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    label 'process_single'
    input:
        file(gc)
    output:
        path("${gc.baseName}")
    script:
    """
        mkdir ${gc.baseName}
        cd ${gc.baseName}
        unzip ../$gc
        cd ..
    """
}
