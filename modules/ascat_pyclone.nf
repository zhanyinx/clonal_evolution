/*
    input: bam file and maf file
    run: ascat + pyclone-vi
    output: pyclone clones
*/


process ascat_pyclone {
    tag "ascat_pyclone"
    input:
        file(samp)
        file(conf)
    script:
    """     
        header="\$(awk -F ',' 'BEGIN{getline; printf "%s", \$1; for(i=2;i<=NF-2;i++) printf ",%s", \$i; print "";}' ${samp})" # get header for sarek input file

        # loop over patients
        for patient in `awk -F ',' 'BEGIN{getline}{print \$1}' ${samp} | sort | uniq`; do 
            
            # get all unique mutations in same patient
            allmaf=()
            for patientmaf in `awk -F ',' 'BEGIN{getline; nrow=NF}{if(\$1 == "'"\$patient"'" && \$3 == 1) print \$nrow}' ${samp} | sort | uniq`; do
                allmaf+=(\$patientmaf)
                awk -F '\\t' '{if(\$10=="SNP") print \$5"-"\$6"-"\$7"-"\$11"-"\$13"-"\$1"-"\$9"-"\$10}' \$patientmaf >> appo
            done

            cat appo | sort | uniq > list_unique_mutations
            sed -i 's/-/ /g' list_unique_mutations
            line=`wc -l list_unique_mutations | awk '{print \$1}'`

            echo "mutation_id sample_id ref_counts alt_counts normal_cn major_cn minor_cn tumour_content" > \${patient}.pyclone.tsv
            sed -i 's/ /\\t/g' \${patient}.pyclone.tsv
            
            for selectedmaf in "\${allmaf[@]}"; do    
                # create patient specific input file and run nf sarek with ascat 
                echo \${header} > tmp.csv
                awk -F ',' '{if(\$1 == "'"\$patient"'") {if(\$3==0 || \$NF == "'"\$selectedmaf"'") {printf "%s", \$1; for(i=2;i<=NF-2;i++) printf ",%s", \$i; print "";}}}' ${samp} >> tmp.csv 

                if [[ ${params.build} == "hg19" ]]; then
                    cp ${params.target} intervals.bed
                    sed -i 's/chr//g' intervals.bed
                    nextflow run nf-core/sarek -r 3.1.1 \
                    -profile singularity \
                    --input tmp.csv \
                    -resume \
                    -c ${conf} \
                    --genome "GATK.GRCh37" \
                    --step ${params.biomarkers_ascat_step} \
                    --intervals intervals.bed 
                else
                    nextflow run nf-core/sarek -r 3.1.1 \
                    -profile singularity \
                    --input tmp.csv \
                    -resume \
                    -c ${conf} \
                    --genome "GATK.GRCh38" \
                    --step ${params.biomarkers_ascat_step} \
                    --intervals ${params.target} 
                fi

                # extract cellularity, tumor/normal sample name and name from dragen
                cellularity="\$(awk 'BEGIN{getline; nrow=NF-1}{if(\$1 == "'"\$patient"'" && \$NF == "'"\$selectedmaf"'") cellularity = \$nrow}END{print cellularity}' ${samp})"
                tumor="\$(awk -F ',' '{if(\$1 == "'"\$patient"'" && \$3 == 1 && \$NF == "'"\$selectedmaf"'") out = \$4}END{print out}' ${samp})"
                normal="\$(awk -F ',' '{if(\$1 == "'"\$patient"'" && \$3 == 0) out = \$4}END{print out}' ${samp})"
                tumor_bam="\$(awk -F ',' '{if(\$1 == "'"\$patient"'" && \$3 == 1 && \$NF == "'"\$selectedmaf"'") out = \$6}END{print out}' ${samp})"
                name=\$tumor
                
                # create pyclone input file
                ascat_file="\$(ls results/variant_calling/ascat/\${tumor}_vs_\${normal}/*cnvs.txt)"

                # create output folder
                if ! [ -d ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/ ]; then
                    mkdir -p ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/
                fi

                # if ascat files do not exist, can't calculate clonal TMB, set it to NA
                if [ -z \${ascat_file} ]; then
                    echo "Clonal TMB: NA" > \${name}.clonalTMB.txt
                    cp \${name}.clonalTMB.txt ${launchDir}/${params.output}/${params.date}/biomarkers/\${name}/
                    continue
                fi
                
                # get number of normal and tumor allele counts
                echo "Chromosome Start_Position End_Position Tumor_Sample_Barcode t_ref_count t_alt_count Genome_Change" > tmp_maf
                tumor_sample=`basename \${selectedmaf}`
                for idx in `seq 1 \$line`; do
                    strings=`awk '{if(NR=="'"\$idx"'"+0.) print \$1":"\$2"-"\$3}' list_unique_mutations`
                    start_end_chr=`awk '{if(NR=="'"\$idx"'"+0.) print \$1,\$2,\$3}' list_unique_mutations`
                    alt=`awk '{if(NR=="'"\$idx"'"+0.) print \$5}' list_unique_mutations`
                    altcounts=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} |  cut -f 5 | tr '[a-z]' '[A-Z]' | fold -w 1 | sort | uniq -c | awk 'BEGIN{ee=0}{if(\$2=="'"\$alt"'") {print \$1; ee=99}}END{if(ee==0) print "0"}'`

                    ref=`awk '{if(NR=="'"\$idx"'"+0.) print \$4}' list_unique_mutations`
                    ref_detected=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} | awk '{print \$3}'`
                    refcounts=`samtools mpileup -f $params.fasta -r \$strings \${tumor_bam} | awk '{print \$4 - "'"\$altcounts"'"+0.}'`
                    echo "\$start_end_chr \${tumor_sample} \$refcounts \$altcounts \$strings.\$ref.\$alt" >> tmp_maf

                done

                sed -i 's, ,\\t,g' tmp_maf
                if ! [ -z \$cellularity ]; then
                    create_input4pyclone.py -as \${ascat_file} -c \${cellularity} -m tmp_maf -o \${patient}.pyclone.tsv --skip_filter
                else
                    ascat_tumor_estimate="\$(ls results/variant_calling/ascat/\${tumor}_vs_\${normal}/*.purityploidy.txt)"
                    create_input4pyclone.py -as \${ascat_file} -m tmp_maf -o \${patient}.pyclone.tsv -ac \${ascat_tumor_estimate} --skip_filter
                fi
                
            done

            # run pyclone
            pyclone-vi fit -i \${patient}.pyclone.tsv -o \${patient}.pyclone.h5 -c 40 -d beta-binomial -r 20
            pyclone-vi write-results-file -i \${patient}.pyclone.h5 -o \${patient}.pyclone.output.tsv

            cp \${patient}.pyclone.output.tsv ${launchDir}/${params.output}/${params.date}/biomarkers/\$patient/
        done
    """

}