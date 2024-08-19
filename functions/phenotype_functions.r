chisquared = function(tumor, selection, search_term = "neuronal system", seq_data = NULL, nmf = NULL, clusterprofiler = NULL, pheno_data = NULL) {
    # function: performs a chi-squared test for patients of the given tumor dataset. The patients are sorted based on the selection term and the search term in the enrich_terms list. If the maximum exposure of a patient is in at least one signature that contains the search term, the patient is considered to be exposed to the search term and put in the patients vector. 

    # Input:
        # tumor: the name of the tumor in the tcga_df DataFrame
        # selection: the column name of the pheno_data DataFrame that should be used for the chi-squared test
        # search_term: the phrase that should be searched for in the enrich_terms, used for sorting the patients
        # seq_data: Optional: the original sequencing data which contains patients as columns. The names of the columns need to be the TCGA patient sample ids. Default: NULL. When no sequencing data is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.
        # nmf: Optional: the nmf object. Default: NULL. When no nmf object is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.
        # clusterprofiler: Optional: the result from the call_clusterprofiler() function. Default: NULL. When no clusterprofiler output object is given, the function calls it based on the tumor argument from the tcga_df defined for the thesis.
        # pheno_data: Optional: the phenotype data. Default: NULL. When no phenotype data is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.


    # Output:
        # test: the result of the chi-squared test for the given selection term and search term. test is a htest object


    dataset = tumor

    if(is.null(seq_data)) {seq_data = get(tcga_df[tcga_df$Dataset == dataset, "seq_data"])}
    if(is.null(pheno_data)) {pheno_data = get(tcga_df[tcga_df$Dataset == dataset, "pheno_data"])}
    if(is.null(nmf)) {nmf = get(tcga_df[tcga_df$Dataset == dataset, "nmf"])}
    if(is.null(clusterprofiler)) {cp = get(tcga_df[tcga_df$Dataset == dataset, "cp"])}


    if((selection %in% colnames(pheno_data)) == FALSE) {
        return(data.frame(p.value = 1))
    }

    specific_selection = c(which(colnames(pheno_data) == "submitter_id.samples"), which(colnames(pheno_data) == selection))

    pheno_data = pheno_data[, specific_selection]

    missing_data = pheno_data[, selection] %>% {which(. == "")}
    if(length(missing_data) != 0) {
        pheno_data = pheno_data[-missing_data,]
    }
    
    patients = seq_data %>% colnames() %>% {gsub("\\.", "-", .)}
    pheno_data = pheno_data$submitter_id.samples %in% patients %>% which() %>% pheno_data[.,]

    neuronal_patients = sort_patients(enrich_terms = cp, search_term = search_term)
    pheno_data$neuronal = ifelse(pheno_data$submitter_id.samples %in% neuronal_patients, 1, 0)

    test = chisq.test(table(pheno_data$neuronal, pheno_data[,selection]))
    return(test)
}