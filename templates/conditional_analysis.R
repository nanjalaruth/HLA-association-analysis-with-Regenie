process ExtractAllele {
    tag "Extract allele with smallest p-value"
    input:
    tuple val(dataset), val(pheno_label), path(saige_output_file)

    output:
    tuple val(dataset), val(pheno_label), path(output)

    script: 
    output = "${pheno_label}_extracted_allele.txt"
    template "extract_allele.py"
}

process Extractdosage {
    input:
    tuple val(dataset), path(genotype_data_file), path(index)

    output:
    tuple val(dataset), path("${dataset}_id_geno_samp.txt") 

    script:
    """
    # Use bcftools to extract genotypes for the specific SNP
    bcftools query -f '%ID[\\t%SAMPLE=%GT]\\n' ${genotype_data_file} > ${dataset}_id_geno_samp.txt
    """
}

process ExtractGenotypes {
    input:
    tuple val(dataset), val(pheno_label), path(extracted_allele_file), path(dosage)

    output:
    tuple val(dataset), val(pheno_label), path(output) 

    script:
    output = "${pheno_label}_id_geno_samp.txt"

    """
    # Extract SNP ID
    SNP_ID=\$(cut -f3 ${extracted_allele_file})

    #Extract dosage for that SNP
    awk -v snp_id="\$SNP_ID" -F '\\t' '\$1 == snp_id {print}' ${dosage} > sig_geno

    #remove id name and transpose row to column
    cut -f 2- sig_geno | cat | tr '\\t' '\\n' > sig_geno_edited 

    #separate sample from genotype
    awk -F '=' '{print \$1 "\\t" \$2}' sig_geno_edited > high_sig_geno

    #replace genotype information
    awk '{if (\$2 == "0|1") \$2 = 1; else if (\$2 == "0|0") \$2 = 0; \
    else if(\$2 == "1|1") \$2 = 2; else \$2 = 1; print}' high_sig_geno | \
    awk -F ' ' '{print \$1 "\\t" \$2}' > high_sig_geno_edited

    #add header
    echo -e "IID\\t${pheno_label}_dosage" > header.tsv
    
    #add header
    cat header.tsv high_sig_geno_edited > ${output}
    """
}

process ModifiedCovariate {
    input:
        tuple val(pheno_label), path(high_sig_geno_edited_file), path(covar_file)

    output:
        tuple val(pheno_label), path("${pheno_label}_updated_covar.txt")
    
    script:
        template "merge_cov.R"

}

process step_1_cond {
    tag "Run logistic regression on ${dataset}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), path(bed), path(bim), path(fam), val(pheno_label), path(phenotype), path(covariate)

    output:
       tuple val(dataset), val(pheno_label), file("${out}_regenie-step1_pred.list"), path(phenotype), path(covariate)

    script:
        out = "${dataset}_${pheno_label}"
        base = bed.baseName

        """
        regenie --step 1 \\
        --bt \\
        --cv 3 \\
        --bed ${base} \\
        --covarFile ${covariate} \\
        --phenoFile ${phenotype} \\
        --phenoCol ${pheno_label} \\
        --bsize 1000 \\
        --lowmem \\
        --lowmem-prefix ${out}_regenie-step1_tmp_rg \\
        --threads 8 \\
        --out ${out}_regenie-step1
        """
}

process step_2_cond {
    tag "Run logistic regression on ${dataset}"
    publishDir "${params.outdir}/Regression_results", mode: 'copy', overwrite: true
    label "bigmem"
    container "/apps/singularity/saige_1.1.6.3.sif"

    input:
        tuple val(dataset), path(bed), path(bim), path(fam), val(pheno_label), path(regenie_step1_pred_list), path(phenotype), path(covariate)

    output:
       tuple val(dataset), val(pheno_label), path("${out}_regenie-step2-bin-out-firth*")

    script:
        out = "${dataset}_${pheno_label}"
        base = bed.baseName

        """
        regenie --step 2 \\
        --bt \\
        --bed ${base} \\
        --covarFile ${covariate} \\
        --phenoFile ${phenotype} \\
        --phenoCol ${pheno_label} \\
        --bsize 1000 \\
        --pred ${regenie_step1_pred_list} \\
        --minMAC 1 \\
        --out ${out}_regenie-step2-bin-out-firth \\
        --threads 1
        """
}

process plot_manhattan_cond {
    publishDir "${params.outdir}/plots", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), val(pheno_label), path(step_2_out) 

    output:
        tuple val(dataset), val(pheno_label), path("${dataset}_${pheno_label}_phewas.pdf")

    script:
        template "manhattan.R"
}

process top_alleles_cond {
    tag "extracting top alleles on ${pheno_label}"
    publishDir "${params.outdir}/top_alleles", mode: 'copy', overwrite: true

    input:
	    tuple val(dataset), val(pheno_label), path(step_2_out) 

    output:
	    tuple val(dataset), val(pheno_label), path("clean_${pheno_label}.tsv")

    script:
        template "top_alleles.R"
}

process combineFiles_cond {
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true
    input:
        path(input_files) 

    output:
        path("combined_file.tsv")

    script:
        template "combine.R"
}


