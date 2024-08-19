call_ensembl = function(gene_ids, biomart = "ensembl", biomart_dataset_name = "hsapiens_gene_ensembl") {
    # function: calls the canonical transcript length of genes via the ensembl biomart

    # Input:
        # gene_ids: a vector of ensembl gene or transcript ids for which to get the transcript length
        # biomart: the biomart to use (default: "ensembl")
        # biomart_dataset_name: the dataset to use (default: "hsapiens_gene_ensembl")

    # Output:
        # unwanted_ids: a vector of the rows with genes that are not in the canonical ensembl biomart
        # transcript_length: a vector of the transcript lengths of the given genes. Genes without a canonical transcript are removed from the output


    gene_ids = str_remove(gene_ids, "\\..*")
    mart = useMart(biomart = biomart, dataset = biomart_dataset_name, host="https://useast.ensembl.org")

    biomart_result = getBM(attributes = c("ensembl_gene_id", "transcript_length", "transcript_is_canonical"),
                          filters = "ensembl_gene_id",
                          values = gene_ids,
                          mart = mart)

    canonical_genes = biomart_result[which(biomart_result$transcript_is_canonical == "1"),]

    retired_genes = which(!(gene_ids %in% biomart_result$ensembl_gene_id))

    duplicated_gene_ids = duplicated(gene_ids) %>% which()

    unwanted_rows = c(duplicated_gene_ids, retired_genes) %>% sort()

    cleaned_genes = gene_ids[-unwanted_rows]

    transcript_length = c()
    for (i in 1:length(cleaned_genes)) {
        gene_id = str_remove(cleaned_genes[i], "\\..*")

        length = which(canonical_genes$ensembl_gene_id == gene_id) %>% {canonical_genes$transcript_length[.]}
        transcript_length = c(transcript_length, length)
    }

    return(list(unwanted_ids = unwanted_rows, transcript_length = transcript_length))
}



feature_selection = function(nmf, k = 3) {
    # function: selects which genes are specific for a single signature of the given rank k and returns it ("specific_W"), as well as all features for the given k ("features")
    
    # Input:
        # nmf: the nmf object from which to extract the features
        # k: the factorization rank to extract the features for

    # Output:
        # features: a matrix with rows as genes and columns as signatures, with 1 if the gene is belonging to the signature and 0 if not
        # specific_W: a matrix with rows as genes and columns as signatures, with the value of the gene in each signature divided by the genes maximum value (normalized)


    features = SignatureSpecificFeatures(nmf, k = k, return_all_features = TRUE)
    colnames(features) = paste0("Signature", 1:k)

    specific = which(rowSums(features) == 1)
    specific_W = WMatrix(nmf, k = k)[specific,]
    specific_W = specific_W/matrixStats::rowMaxs(specific_W)
    colnames(specific_W) = paste0("Signature", 1:k)

    return(list(features = features, W = specific_W))
}



normalize_H = function(nmf, k) {
    # function: normalizes the H matrix of the given nmf object for the given rank k

    # Input:
        # nmf: the nmf object from which to extract the H matrix
        # k: the factorization rank to extract the H matrix for

    # Output:
        # normalized_H: the normalized H matrix


    H_Matrix = HMatrix(nmf, k = k_i)
    colmax = matrixStats::colMaxs(H_Matrix)
    H_Matrix = H_Matrix / rep(colmax, each = nrow(H_Matrix))
    rownames(H_Matrix) = paste0("Signature_", 1:k_i)
    return(H_Matrix)
}
