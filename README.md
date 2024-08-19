# Identification of Neuronal Signatures in 16 different Cancers using Non-Negative Matrix Factorization
This GitHub Repository contains the code and results of the Bachelor´s thesis of Niklas Markus Winata. It allows you to recreate all results and use it on additional tumor data or to identify different signatures.

## 1 Acknowledgements
My thesis was writen at the CROmLab at University of Heidelberg. Many thanks to Prof. Dr. Carl Herrmann and his group for their support in this journey.

## 2 Data availability
All data used in this analysis is found in the `data` folder. The contents of the `tumor_tissue`, `phenotype` and `survival` folder were obtained from the `[GDC Xena Hub]`, while the content of the `healthy_tissue` folder was obtained from the `[GTEx Portal]`.
Additionally, all RObjects and figures resulting from the analysis is uploaded in the `analysis` folder.

## 3 Installation
### 3.1 Conda Environment
The conda environment used for all analyses was exported as a yaml file. This can be used to replicate the work environment. The actual analysis is to be performed in R.

### 3.2 ButchR
The non-negative matrix factorization (NMF) is based on the [ButchR](https://github.com/wurst-theke/ButchR) package by [Quintero _et al._](https://doi.org/10.1093/biomethods/bpaa022). In order to properly function, the ButchR package needs specific versions of `reticulate`and `tensorflow`, which are specified in the conda environment.

## 4 Loading necessary libraries, functions and data
The `libraries_and_data.r` file loads all necessary libraries needed for the analysis, as well as all functions that were created during the thesis and the data that was used in the analysis. Therefore, when starting the analysis, one should start by loading this file in R:
```{r}
source("/path/to/file/libraries_and_data.r")
```

### 4.1 Costum Functions
All functions that were created for this project are in the `functions` folder. They will be used in this analysis and are automatically available when loading the `libraries_and_data.r` file. There is a description for every function at the start of its definition in their respective functions file. 

## 5 Prepare Data
The data used in this project is either TCGA FPKM gene expression RNAseq data for 16 tumor tissues obtained from the `[GDC Xena Hub]`, while gene count healthy tissue data was obtained from the `[GTEx Portal]`. While the FPKM gene expression data is ready to use in the `run_NMF_tensor()` function of the ButchR package, the gene count data from GTEx first had to be transformed to match the log2(FPKM + 1) transformation of the tumor data.
In my thesis, the gene count data was transformed into log2(FPKM + 1) data by using the `ensembl` database to get the transcript length of each gene (using the costum `call_ensembl()` function). For the transformation, the first three columns are ignored in the following code, as these columns are gene ids and descriptions in `GTEx` gene count data. If other gene count data is used, all non-numeric columns need to be ommited. If the healthy tissue data is already in FPKM format, the transformation is not necessary.

```{r}
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

```{r}
tumor_nmf = run_NMF_tensor(X = TCGA_FPKM_data[,-1], ranks = 2:20, method = "NMF", n_initializations = 30, extract_features = TRUE)
```

### 6.2 For healthy tissue data
For healthy tissue, the `fpkm_tissue` data obtained in section 5 is used. As it only contains numeric values, no column needs to be removed.

```{r}
tissue_nmf = run_NMF_tensor(X = fpkm_tissue, ranks = 2:20, method = "NMF", n_initializations = 30, extract_features = TRUE)
```

### 6.3 Normalize W Matrix
Before starting with further analysis, the W matrices of each NMF should be normalized using the `normalizeW()` function of the ButchR package:

```{r}
nmf = normalizeW(nmf)
```

