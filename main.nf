nextflow.enable.dsl=2

include { step_1; step_2; 
        plot_manhattan; combineFiles;  top_alleles
} from './modules/regenie_analysis'


workflow{

     //step 1
     cov_pheno_ch = Channel.fromList(params.cov_pheno)
     plink_ch = Channel.fromList(params.plink_files)
     input = plink_ch
        .combine(cov_pheno_ch)
     //input.view()
    step_1(input)
  
     //step2
    step1_out = step_1.out
    // step1_out.view()
    vcf_ch = Channel.fromList(params.plink_files)
        .combine(step1_out, by:0)
    // vcf_ch.view()
    step_2(vcf_ch)

    // step 3
    in = step_2.out
        .map { dataset, pheno_label, files -> [ dataset, pheno_label, files[1]]}
    //n.view()
    plot_manhattan(in) 

//     // step 4 Combine files
    top_in = step_2.out
        .map { dataset, pheno_label, files -> [ dataset, pheno_label, files[1]]}
    top_alleles(top_in)

    input = top_alleles.out
        .map { dataset, pheno_label, files -> [files]}
        .collect()
    // input.view()
    combineFiles(input)


// //------------------------------------------------------
//     // Step B
// //     // Conditional analysis round one
//     input = top_alleles.out
//     //input.view()
//     ExtractAllele(input)

    // geno = Channel.fromList(params.whole_ckb_vcf)
//     //geno.view()
//     Extractdosage(geno)

//     allele = ExtractAllele.out
//     dosage = Extractdosage.out
//     in = allele
//         .combine(dosage, by:0)
//     //in.view()
//     ExtractGenotypes(in)

//     dosage = ExtractGenotypes.out
//         .map{dataset, id, geno -> [id, geno]}
//     cov = Channel.fromList(params.cov_pheno)
//     input = dosage
//         .combine(cov, by:0)
//     //input.view()
//     ModifiedCovariate(input)

// //     // Run SAIGE after conditional analysis
//      cov_pheno = ModifiedCovariate.out
//      input = saige_spop_ch
//         .combine(cov_pheno)
//         .combine(sparse_in)
//     //input.view()
//     saige_logistic_cond_1(input)

//     step1_out = saige_logistic_cond_1.out
//     reg_2_ch = vcf_ch.combine(ids_ch)
//           .combine(step1_out, by:0)
//           .combine(sparse_in)
//     //reg_2_ch.view()
//     saige_logistic_cond_2(reg_2_ch)

// //     plot_in = saige_logistic_cond_2.out
// //     plot_cond_manhattan(plot_in) 

//     // Combine files
//     cond_in = saige_logistic_cond_2.out
//     top_cond_alleles(cond_in)

//     inpt_cond = top_cond_alleles.out
//      .map{dataset, pheno, file -> [file]}
//      .collect()
//     //inpt.view()
//     combineFiles_cond(inpt_cond)

// // //---------------------------------------
// // // Conditional analysis round 2
//     input = top_cond_alleles.out
//     //input.view()
//     ExtractAllele_1(input)

// //     //geno = Channel.fromList(params.whole_ckb_vcf)
// //     //geno.view()
// //     //Extractdosage_1(geno)

//     allele_1 = ExtractAllele_1.out
//     dosage_1 = Extractdosage.out
//     in_1 = allele_1
//         .combine(dosage_1, by:0)
//     //in.view()
//     ExtractGenotypes_1(in_1)

//     dosage = ExtractGenotypes_1.out
//         .map{dataset, id, geno -> [id, geno]}
//     cov = ModifiedCovariate.out
//     //cov.view()
//     input = dosage
//         .combine(cov, by:0)
//     //input.view()
//     ModifiedCovariate_1(input)

//     // Run SAIGE after conditional analysis
//     cov_pheno_2 = ModifiedCovariate_1.out
//     in_2 = saige_spop_ch
//         .combine(cov_pheno_2)
//         .combine(sparse_in)
//     //in_2.view()
//     saige_logistic_cond_1_1(in_2)

//     step1_1_out = saige_logistic_cond_1_1.out
//     reg_2_1_ch = vcf_ch.combine(ids_ch)
//           .combine(step1_1_out, by:0)
//           .combine(sparse_in)
//     //reg_2_ch.view()
//     saige_logistic_cond_2_1(reg_2_1_ch)

// //     //Plot
// //     plot_in = saige_logistic_cond_2_1.out
// //     plot_cond_manhattan_1(plot_in) 

//     // Combine files
//     cond_in = saige_logistic_cond_2_1.out
//     top_cond_alleles_1(cond_in)

//     inpt_cond = top_cond_alleles_1.out
//      .map{dataset, pheno, file -> [file]}
//      .collect()
// //     //inpt.view()
//     combineFiles_cond_1(inpt_cond)
}