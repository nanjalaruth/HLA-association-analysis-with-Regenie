process step_1 {
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

process step_2 {
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
        --threads 8
        """
}

process plot_manhattan {
    publishDir "${params.outdir}/plots", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), val(pheno_label), path(step_2_out) 

    output:
        tuple val(dataset), val(pheno_label), path("${dataset}_${pheno_label}_phewas.pdf")

    script:
        template "manhattan.R"
}

process top_alleles {
    tag "extracting top alleles on ${pheno_label}"
    publishDir "${params.outdir}/top_alleles", mode: 'copy', overwrite: true

    input:
	    tuple val(dataset), val(pheno_label), path(step_2_out) 

    output:
	    tuple val(dataset), val(pheno_label), path("clean_${pheno_label}.tsv")

    script:
        template "top_alleles.R"
}


process combineFiles {
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true
    input:
        path(input_files) 

    output:
        path("combined_file.tsv")

    script:
        template "combine.R"
}


