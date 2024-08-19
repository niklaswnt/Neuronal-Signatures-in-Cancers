# Identification of Neuronal Signatures in 16 different Cancers using Non-Negative Matrix Factorization
This GitHub Repository contains the code and results of the Bachelor´s thesis of Niklas Markus Winata. It allows you to recreate all results and use it on additional tumor data or to identify different signatures.

## 1 Acknowledgements
My thesis was writen at the CROmLab at University of Heidelberg. Many thanks to Prof. Dr. Carl Herrmann and his group for their support in this journey.

## 2 Data availability
All data used in this analysis is found in the `data` folder. The contents of the `tumor_tissue`, `phenotype` and `survival` folder were obtained from the [GDC Xena Hub](https://gdc.xenahubs.net/), while the content of the `healthy_tissue` folder was obtained from the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression).
Additionally, all RObjects and figures resulting from the analysis is uploaded in the `analysis` folder.

## 3 Installation
### 3.1 Conda Environment
The conda environment used for all analyses was exported as a yaml file. This can be used to replicate the work environment. The actual analysis is to be performed in R.

### 3.2 ButchR
The non-negative matrix factorization (NMF) is based on the [ButchR](https://github.com/wurst-theke/ButchR) package by [Quintero _et al._](https://doi.org/10.1093/biomethods/bpaa022). In order to properly function, the ButchR package needs specific versions of `reticulate`and `tensorflow`, which are specified in the conda environment.

## 4 Loading necessary libraries, functions and data
The `libraries_and_data.r` file loads all necessary libraries needed for the analysis, as well as all functions that were created during the thesis and the data that was used in the analysis. Therefore, when starting the analysis, one should start by loading this file in R:
``` r 
source("/path/to/file/libraries_and_data.r")
```

However, the `libraries_and_data.r` file needs to be adjusted to fit the path for all data that it loads.

### 4.1 custom Functions
All functions that were created for this project are in the `functions` folder. They will be used in this analysis and are automatically available when loading the `libraries_and_data.r` file. There is a description for every function at the start of its definition in their respective functions file. 

## 5 Prepare Data
The data used in this project is either TCGA FPKM gene expression RNAseq data for 16 tumor tissues obtained from the [GDC Xena Hub](https://gdc.xenahubs.net/), while gene count healthy tissue data was obtained from the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression). While the FPKM gene expression data is ready to use in the `run_NMF_tensor()` function of the ButchR package, the gene count data from GTEx first had to be transformed to match the log2(FPKM + 1) transformation of the tumor data. \\
In my thesis, the gene count data was transformed into log2(FPKM + 1) data by using the `ensembl` database to get the transcript length of each gene (using the custom `call_ensembl()` function). For the transformation, the first three columns are ignored in the following code, as these columns are gene ids and descriptions in GTEx gene count data. If other gene count data is used, all non-numeric columns need to be ommited. If the healthy tissue data is already in FPKM format, the transformation is not necessary.

``` r 
# Read GCT (Example)
tissue_data = read.delim(file="/path/to/file/gene_reads_tissue.gct", skip=2)

ensembl = call_ensembl(gene_ids = tissue_data$Name)
transcript_length = ensembl$transcript_length
tissue_normal = tissue_data[-ensembl$unwanted_ids,]

# Transformation
fpkm_tissue = calculateFPKM(tissue_normal[,-(1:3)], transcript_length)
fpkm_tissue = log2(fpkm_tissue + 1)
```

## 6 Perform NMF
The NMF is performed using the `run_NMF_tensor()` function of the ButchR package. For that, the input data needs to be a matrix only containing numerics. For the parameters, the ranks 2 to 20 were chosen and the number of initializations was set to 30. All other parameters are the default values of the `run_NMF_tensor()` function.

### 6.1 For TCGA FPKM gene expression data
The TCGA FPKM data is already ready for running NMF, but the gene id column needs to be removed before starting the analysis.

``` r 
tumor_nmf = run_NMF_tensor(X = TCGA_FPKM_data[,-1], ranks = 2:20, method = "NMF", n_initializations = 30, extract_features = TRUE)
```

### 6.2 For healthy tissue data
For healthy tissue, the `fpkm_tissue` data obtained in section 5 is used. As it only contains numeric values, no column needs to be removed.

``` r 
tissue_nmf = run_NMF_tensor(X = fpkm_tissue, ranks = 2:20, method = "NMF", n_initializations = 30, extract_features = TRUE)
```

### 6.3 Normalize W Matrix
Before starting with further analysis, the W matrices of each NMF should be normalized using the `normalizeW()` function of the ButchR package:

``` r 
nmf = normalizeW(nmf)
```

## 7 Heatmaps for W and H matrices
To take a first look at the product of the NMF, heatmaps were created for all W and H matrices. These can be found in the `heatmaps` folder.

### 7.1 W Heatmaps
The W heatmaps show which gene belongs to which signature. The heatmaps only show signature-specific genes, so genes that belong to only one signature. They are sorted based on the `SignatureSpecificFeatures()` function of the ButchR package, and the custom `feature_selection()` and `WHeatmap()` functions are based on ButchR as well. 

``` r 
for (k_i in ranks) {
    featureselection = feature_selection(nmf = nmf, k = k_i)
    specific_W = featureselection[[2]]
    assign(paste0("WHeatmap_", k_i), WHeatmap(specific_W, name = paste0("k = ", k_i)))
}
```

The result should look like the following heatmap `(for ..., k = )`
![W](analysis/W_Heatmaps/breast_w_heatmap.png)

### 7.2 H Heatmaps
The H heatmaps show which patient is influenced by which signature. For that, the exposure is normalized by the maximum exposure per patient (column). The `HHeatmap()` function is based on ButchR.

``` r 
for (k_i in ranks) {
    H_Matrix = normalize_H(nmf, k = k_i)
    assign(paste0("HHeatmap_", k_i), HHeatmap(H_Matrix, name = paste0("k = ", k_i)))
}
```

The result should look like the following heatmap `(for ..., k = )`

## 8 Riverplots for signature stability
Using the `riverplot` package and the custom `scaleRiverplot()` function, riverplots are generated to look at the stability of the signatures for different ranks.

``` r
rp = generateRiverplot(nmf, ranks = 3:20)
rp = scaleRiverplot(rp, min_rank = 3, max_rank = 20)
riverplot(rp)
```

The result should look like the following riverplot `(for ...)`

## 9 Pathway Enrichment Analysis
In order to identify the biological relevance of each signature, pathway enrichment analysis is performed using the `Reactome` database via the `clusterProfiler` package. A list is created for every NMF, which contains an enrichResult for every signature and rank. This is achieved by using the custom `call_clusterprofiler()` function. 

``` r 
enriched_terms = call_clusterprofiler(nmf, gene_id_column, ranks = 3:20, pvalue = 0.05, skip = NULL, geneType = "ENSEMBL")
```

The result is in the following format:
``` r 
list("k1" = list("Signature_1" = enrichResult object, "Signature_2" = enrichResult object), "k2" = list("Signature_1" = enrichResult object, "Signature_2" = enrichResult object), ...)
```

For my thesis, signatures were checked for neuronal enrichment. These neuronal signatures were defined as signatures with enrichment (p-value < 0.05) for the "Neuronal System" superpathway of the Reactome database.

## 10 Comparing tumor patients with and without neuronal enrichment
The enrichment results were used to compare survival and phenotypic characteristics between tumor patients with neuronal enrichment and those without. The patients were divided depending on their maximum exposure to the signatures (as seen in the H heatmaps). If a patient had a maximum exposure for a neuronal signature at least once, they were counted as a neuronal patient, otherwise as a non-neuronal patient. These two groups were then compared.

### 10.1 Survival Analysis
The two groups are compared for their survival in all tumors. For this, Kaplan-Meier plots are created. The custom `generateSurvivalplot()` function creates such a plot after sorting the patients into two groups. The type of enrichment can be specified with the `search_term` argument, which allows for comparison of patients with enrichment different from neuronal enrichment. The data needed can be specified manually, or the function will call it based the `tcga_df` object defined for this thesis. This automatic call only works for the 16 tumors used in the thesis, unless the `tcga_df` object is manually expanded.
 
``` r 
# Example for automatic call
generateSurvivalplot(tumor = "BRCA", search_term = "neuronal system", seq_data = NULL, nmf = NULL, clusterprofiler = NULL, surv_data = NULL)
```

The result should look like the following Kaplan-Meier plot `(for ...)`

### 10.2 Phenotype comparison
Several phenotypic characteristics were compared as well between the two groups. Chi-squared tests are performed for a specified characteristic. This can be achieved using the custom `chisquare()` function. The input is identical to the `generateSurvivalplot()` function, but instead of the `surv_data` argument, it has the `pheno_data` argument. The output is a `htest` object, which includes the p-value and a contingency table.

``` r 
chisquared(tumor = dataset, selection = "tumor_stage.diagnoses")
chisquared(tumor = dataset, selection = "tumor_stage.diagnoses")$p.value    # Returns the p-value of the test
chisquared(tumor = dataset, selection = "tumor_stage.diagnoses")$observed   # Returns the contingency table of the test
```


If the chi-squared test should be performed for every tumor in the `tcga_df` object and checked for every signficiant difference, the following code can be used:

``` r 
selection = "primary_diagnosis.diagnoses"

p_values = sapply(tcga_df$Dataset, function(x) {
                    dataset = x
                    p = chisquared(tumor = dataset, selection)$p.value
                    return(p)
            })

sgnfct_p = which(p_values < 0.05)

contingency_tables = sapply(tcga_df$Dataset[sgnfct_p], function(x){
                        dataset = x
                        p = chisquared(tumor = dataset, selection = selection)$observed
                        return(p)
})

p_values[sgnfct_p]
contingency_tables
```
