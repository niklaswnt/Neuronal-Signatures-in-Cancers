WHeatmap = function(W, name = "W Matrix") {
    # function: creates a heatmap for the W matrix. Optimal use with specific_W Matrix obtained from the feature_selection function. Used for all W Heatmaps in the thesis.

    # Input:
        # W: the W matrix to create the heatmap for
        # name: the name of the heatmap, displayed above the legend. Default: "W Matrix"

    # Output:
        # WHeatmap: a Heatmap object (from ComplexHeatmap package) for the W matrix

    WHeatmap = Heatmap(W,
                        col = inferno(100), 
                        name = name,
                        clustering_distance_columns = 'pearson',
                        show_column_dend = TRUE,
                        show_column_names = TRUE,
                        show_row_names = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = FALSE)
    return(WHeatmap)
}



HHeatmap = function(H_Matrix, name = "Exposure") {
    # function: creates a heatmap for the H matrix. Optimal use with H_Matrix Matrix obtained from the normalize_H function. Used for all H Heatmaps in the thesis.

    # Input:
        # H_Matrix: the H matrix to create the heatmap for
        # name: the name of the heatmap, displayed above the legend. Default: "Exposure"

    # Output:
        # h_heatmap: a Heatmap object (from ComplexHeatmap package) for the H matrix


    h_heatmap = Heatmap(H_Matrix,
                        col = viridis(100),
                        name = name,
                        clustering_distance_columns = 'pearson',
                        show_column_dend = TRUE,
                        show_column_names = FALSE,
                        show_row_names = TRUE,
                        cluster_rows = FALSE)
    return(h_heatmap)
}