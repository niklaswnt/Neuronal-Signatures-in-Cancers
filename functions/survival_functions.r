generateSurvivalplot = function(tumor, search_term, seq_data = NULL, nmf = NULL, clusterprofiler = NULL, surv_data = NULL) {
    # function: generates a Kaplan-Meier plot for the survival data of a given tumor dataset. The patients are sorted based on the search term in the enrich_terms list. If the maximum exposure of a patient is in at least one signature that contains the search term, the patient is considered to be exposed to the search term and put in the patients vector.

    # Input:
        # tumor: the name of the tumor in the tcga_df DataFrame
        # search_term: the phrase that should be searched for in the enrich_terms, used for sorting the patients
        # seq_data: Optional: the original sequencing data which contains patients as columns. The names of the columns need to be the TCGA patient sample ids. Default: NULL. When no sequencing data is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.
        # nmf: Optional: the nmf object. Default: NULL. When no nmf object is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.
        # clusterprofiler: Optional: the result from the call_clusterprofiler() function. Default: NULL. When no clusterprofiler output object is given, the function calls it based on the tumor argument from the tcga_df defined for the thesis.
        # surv_data: Optional: the survival data. Default: NULL. When no survival data is given, the functions calls it based on the tumor argument from the tcga_df defined for the thesis.

    # Output:
        # survplot: a Kaplan-Meier plot for the survival data of the given tumor dataset, split into patients that are exposed to the search term and patients that are not exposed to the search term


    dataset = tumor

    if(is.null(seq_data)) {seq_data = get(tcga_df[tcga_df$Dataset == dataset, "seq_data"])}
    if(is.null(surv_data)) {surv_data = get(tcga_df[tcga_df$Dataset == dataset, "surv_data"])}
    if(is.null(nmf)) {nmf = get(tcga_df[tcga_df$Dataset == dataset, "nmf"])}
    if(is.null(clusterprofiler)) {cp = get(tcga_df[tcga_df$Dataset == dataset, "cp"])}

    patients = seq_data %>% colnames() %>% {gsub("\\.", "-", .)}
    surv_data = surv_data$sample %in% patients %>% which() %>% surv_data[.,]

    neuronal_patients = sort_patients(enrich_terms = cp, search_term = search_term)
    surv_data$neuronal = ifelse(surv_data$sample %in% neuronal_patients, 1, 0)
            
    sfit = survfit(Surv(OS.time, OS) ~ neuronal, data = surv_data)
    
    survplot = ggsurvplot(sfit, pval=TRUE, conf.int=TRUE, risk.table=TRUE, 
            legend.labs=c("No", "Yes"), legend.title="Neuronal",  
            palette=c("dodgerblue2", "orchid2"), 
            title=paste0("Kaplan-Meier Curve for ", dataset, " Survival"))

    return(survplot)
}