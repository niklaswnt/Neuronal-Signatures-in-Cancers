get_gene_ids = function(nmf_data, gene_id_column, k = 3, k_i = 1) {
    # function: returns the gene ids of the genes that are specific for the given signature k_i of the given rank k, based on the feature_selection function.

    # Input:
        # nmf_data: the nmf object from which to extract the features
        # gene_id_column: the column that contains the gene ids. The order of the gene ids must match the order of the genes in the nmf object
        # k: the factorization rank to extract the features for
        # k_i: the signature to extract the gene ids for

    # Output:
        # signature_gene_ids: the gene ids of the genes that are specific for the given signature k_i of the given rank k.

    signature_gene_ids = feature_selection(nmf = nmf_data, k = k)$features %>%
                            {which(rowSums(.) == 1 & .[,k_i] == 1)} %>%
                            gene_id_column[.] %>%
                            str_remove("\\..*") %>%
                            unique()

    return(signature_gene_ids)
}

call_clusterprofiler = function(nmf, gene_id_column, ranks = 3:20, pvalue = 0.05, skip = NULL, geneType = "ENSEMBL") {
    # function: calls the enrichPathway function from the clusterProfiler package for a range of ranks k and signatures k_i for a specific nmf object

    # Input:
        # nmf: the nmf object from which to extract the features
        # gene_id_column: the column that contains the gene ids. The order of the gene ids must match the order of the genes in the nmf object
        # ranks: the factorization ranks to extract the features for
        # pvalue: the pvalue cutoff for the enrichPathway function. Default: 0.05
        # skip: a list of indices to skip for each rank. The indices are the signatures that should be skipped for the given rank. Input in the format of list("k1" = c(k_i1, k_i2), "k2" = c(k_i3, k_i4), ...). Used to skip signatures that crash the enrichPathway function (usually when they are too short). Default: NULL
        # geneType: the type of the gene ids. Needed for the bitr function of the clusterProfiler package. Default: "ENSEMBL"

    # Output:
        # terms: a list of lists of the enriched pathways for each signature of each rank. The list is structured as list("k1" = list("Signature_1" = enrichResult object, "Signature_2" = enrichResult object), "k2" = list("Signature_1" = enrichResult object, "Signature_2" = enrichResult object), ...)

    terms = lapply(ranks, function(k) {
        
        signatures = 1:k
        
        if (!is.null(skip)) {
            skip_indices = skip[[as.character(k)]]
            if (!is.null(skip_indices)) {
                signatures = signatures[!signatures %in% skip_indices]
            }
        }

        rank_terms = lapply(signatures, function(k_i) {
            print(paste0("k",k, " Sig", k_i))
            gene_list = get_gene_ids(nmf = nmf, gene_id_column = gene_id_column, k = k, k_i = k_i)
            entrez_ids = bitr(gene_list, fromType = geneType, toType = "ENTREZID", OrgDb=org.Hs.eg.db)
            gene_ids = entrez_ids$ENTREZID
            reactome_enrich = enrichPathway(gene = gene_ids, 
                                            pvalueCutoff = pvalue, 
                                            readable = TRUE, 
                                            organism = "human")
            return(reactome_enrich)
        })

        names(rank_terms) = paste0("Signature_", signatures)
        
        return(rank_terms)
    })

    names(terms) = paste0("k", ranks)

    return(terms)
}



calculate_enrichment_score = function(enrich_terms, ranks = 3:(length(enrich_terms)+2), search_term, mute_print = FALSE) {
    # function: calculates the enrichment score for a given search term of a call_clusterprofiler output

    # Input:
        # enrich_terms: the output of the call_clusterprofiler function
        # ranks: the ranks for which the enrichment score should be calculated. Default: 3:(length(enrich_terms)+2). The default is used to calculate the enrichment score for all ranks in the enrich_terms list
        # search_term: the phrase that should be searched for in the enrich_terms and for which the enrichment score should be calculated.
        # mute_print: if TRUE, the function will not print the enrichment score and description for every match with the search term. Default: FALSE

    # Output:
        # This function prints the enrichment score. Output clarifies for which signature and rank
        # This function prints the description of the first match for the search term
        # total_score: the total enrichment score for the given search term. Enrichment Score for each signature is based on the -log10 of the p.adjust value of the first match for the search term in the enrich_terms
        
    
    total_score = 0

    for (k in ranks) {
        for(k_i in 1:length(enrich_terms[[k-2]])) {
            if (dim(as.data.frame(enrich_terms[[k-2]][[k_i]]))[1] == 0) {next}

            enrichResult = enrich_terms[[k-2]][[k_i]] %>% as.data.frame() %>% arrange(p.adjust)
            
            matches = enrichResult %>% .$Description
            which_neuro = sapply(matches, function(x) grepl(search_term, x, ignore.case = TRUE)) %>% which()
            
            if(length(which_neuro > 0)) {
                first_neuro = which_neuro[1]
                score = enrichResult[first_neuro, "p.adjust"] %>% -log10(.)
                total_score = total_score + score
                
                if(!mute_print) {
                    print(paste0("Score for ", gsub("_", " ", names(enrich_terms[[k-2]])[k_i]), " for k = ", k, " is ", score))
                    print(paste0("Description: ", names(first_neuro)))
                }
            }
        }
    }

    return(total_score)
}



sort_patients = function(enrich_terms, ranks = 3:(length(enrich_terms)+2), search_term = "neuronal system") {
    # function: sorts the patients based on the search term in the enrich_terms list. If the maximum exposure of a patient is in at least one signature that contains the search term, the patient is considered to be exposed to the search term and put in the patients vector.

    # Input:
        # enrich_terms: the output of the call_clusterprofiler function
        # ranks: the ranks for which the search term should be checked. Default: 3:(length(enrich_terms)+2). The default is used to check for the search term in all ranks in the enrich_terms list
        # search_term: the phrase that should be searched for in the enrich_terms. Default: "neuronal system"

    # Output:
        # patients: a vector of the patients that are exposed to the search term. The patients are identified by the maximum exposure to a signature that contains the search term at least once.
    
    
    patients = c()

    for(k in ranks) {
        for(k_i in 1:length(enrich_terms[[k-2]])) {
            if (dim(as.data.frame(enrich_terms[[k-2]][[k_i]]))[1] == 0) {next}

            enrichResult = enrich_terms[[k-2]][[k_i]] %>% as.data.frame() %>% arrange(p.adjust)
            
            matches = enrichResult %>% .$Description
            contains_search_term = any(sapply(matches, function(x) grepl(search_term, x, ignore.case = TRUE)))
            
            if(contains_search_term == TRUE) {
                
                signature = gsub("Signature_", "", names(enrich_terms[[k-2]])[k_i]) %>% as.numeric()

                H_Matrix = normalize_H(nmf = nmf, k = k)

                patients = c(patients, colnames(H_Matrix)[which(H_Matrix[signature,] == 1)])
            }
        }
    }

    patients = unique(patients) %>% {gsub("\\.", "-", .)}

    return(patients)
}