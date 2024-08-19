library(BiocStyle)
library(knitr)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)
library(dplyr)

reticulate::use_condaenv("butchr", conda="/path/to/conda/miniconda3/bin/conda",required = TRUE)
library(reticulate)
library(ButchR)

library(stringr)
library(biomaRt)
library(enrichR)
library(magick)
library(scuttle)

library(riverplot)
library(nnls)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

library(survival)
library(survminer)

library(ggplot2)
library(gridExtra)

library(writexl)

source("/path/to/functions/functions.r")
source("/path/to/functions/matrix_functions.r")
source("/path/to/functions/heatmap_functions.r")
source("/path/to/functions/riverplot_functions.r")
source("/path/to/functions/clusterprofiler_functions.r")
source("/path/to/functions/survival_functions.r")
source("/path/to/functions/phenotype_functions.r")

# Load all Normal Tissue-related Data
load("/path/to/healthy_tissue/pancreas_normal.RData")
load("/path/to/healthy_tissue/prostate_normal.RData")
load("/path/to/healthy_tissue/colon_sigmoid_normal.RData")
load("/path/to/healthy_tissue/colon_transverse_normal.RData")
load("/path/to/healthy_tissue/lung_normal.RData")

load("/path/to/healthy_tissue/bladder_normal.RData")
load("/path/to/healthy_tissue/blood_normal.RData")
load("/path/to/healthy_tissue/breast_normal.RData")
load("/path/to/healthy_tissue/ectocervix_normal.RData")
load("/path/to/healthy_tissue/endocervix_normal.RData")
load("/path/to/healthy_tissue/liver_normal.RData")
load("/path/to/healthy_tissue/ovary_normal.RData")
load("/path/to/healthy_tissue/skin_no_sun_normal.RData")
load("/path/to/healthy_tissue/skin_sun_normal.RData")
load("/path/to/healthy_tissue/stomach_normal.RData")
load("/path/to/healthy_tissue/testis_normal.RData")
load("/path/to/healthy_tissue/uterus_normal.RData")



load("/path/to/analysis/nmf/pancreas_nmf.RData")
pancreas_nmf = normalizeW(pancreas_nmf)
load("/path/to/analysis/nmf/prostate_nmf.RData")
prostate_nmf = normalizeW(prostate_nmf)
load("/path/to/analysis/nmf/colon_sigmoid_nmf.RData")
colon_sigmoid_nmf = normalizeW(colon_sigmoid_nmf)
load("/path/to/analysis/nmf/colon_transverse_nmf.RData")
colon_transverse_nmf = normalizeW(colon_transverse_nmf)
load("/path/to/analysis/nmf/lung_nmf.RData")
lung_nmf = normalizeW(lung_nmf)

load("/path/to/analysis/nmf/bladder_nmf.RData")
bladder_nmf = normalizeW(bladder_nmf)
load("/path/to/analysis/nmf/blood_nmf.RData")
blood_nmf = normalizeW(blood_nmf)
load("/path/to/analysis/nmf/breast_nmf.RData")
breast_nmf = normalizeW(breast_nmf)
load("/path/to/analysis/nmf/ectocervix_nmf.RData")
ectocervix_nmf = normalizeW(ectocervix_nmf)
load("/path/to/analysis/nmf/endocervix_nmf.RData")
endocevix_nmf = normalizeW(endocervix_nmf)
load("/path/to/analysis/nmf/liver_nmf.RData")
liver_nmf = normalizeW(liver_nmf)
load("/path/to/analysis/nmf/ovary_nmf.RData")
ovary_nmf = normalizeW(ovary_nmf)
load("/path/to/analysis/nmf/skin_no_sun_nmf.RData")
skin_no_sun_nmf = normalizeW(skin_no_sun_nmf)
load("/path/to/analysis/nmf/skin_sun_nmf.RData")
skin_sun_nmf = normalizeW(skin_sun_nmf)
load("/path/to/analysis/nmf/stomach_nmf.RData")
stomach_nmf = normalizeW(stomach_nmf)
load("/path/to/analysis/nmf/testis_nmf.RData")
testis_nmf = normalizeW(testis_nmf)
load("/path/to/analysis/nmf/uterus_nmf.RData")
uterus_nmf = normalizeW(uterus_nmf)


# Load all Brain-related Data
load("/path/to/healthy_tissue/amygdala_normal.RData")
load("/path/to/healthy_tissue/hippocampus_normal.RData")
load("/path/to/data/fpkm/TCGA_GBM_seq.RData")

load("/path/to/analysis/nmf/amygdala_nmf.RData")
amydala_nmf = normalizeW(amygdala_nmf)
load("/path/to/analysis/nmf/hippocampus_nmf.RData")
hippocampus_nmf = normalizeW(hippocampus_nmf)
load("/path/to/analysis/nmf/GBM_nmf.RData")
GBM_nmf = normalizeW(GBM_nmf)

# Load all Tumor-related Data
load("/path/to/data/fpkm/TCGA_PRAD_seq.RData")
load("/path/to/data/fpkm/TCGA_PAAD_seq.RData")
load("/path/to/data/fpkm/TCGA_COAD_seq.RData")
load("/path/to/data/fpkm/TCGA_LUAD_seq.RData")

load("/path/to/data/fpkm/TCGA_BLCA_seq.RData")
load("/path/to/data/fpkm/TCGA_BRCA_seq.RData")
load("/path/to/data/fpkm/TCGA_CESC_seq.RData")
load("/path/to/data/fpkm/TCGA_HNSC_seq.RData")
load("/path/to/data/fpkm/TCGA_LAML_seq.RData")
load("/path/to/data/fpkm/TCGA_LIHC_seq.RData")
load("/path/to/data/fpkm/TCGA_READ_seq.RData")
load("/path/to/data/fpkm/TCGA_SKCM_seq.RData")
load("/path/to/data/fpkm/TCGA_STAD_seq.RData")
load("/path/to/data/fpkm/TCGA_TGCT_seq.RData")
load("/path/to/data/fpkm/TCGA_UCEC_seq.RData")


load("/path/to/analysis/nmf/PRAD_nmf.RData")
load("/path/to/analysis/nmf/PAAD_nmf.RData")
load("/path/to/analysis/nmf/COAD_nmf.RData")
load("/path/to/analysis/nmf/LUAD_nmf.RData")

load("/path/to/analysis/nmf/blca_nmf.RData")
blca_nmf = normalizeW(blca_nmf)
load("/path/to/analysis/nmf/brca_nmf.RData")
brca_nmf = normalizeW(brca_nmf)
load("/path/to/analysis/nmf/cesc_nmf.RData")
cesc_nmf = normalizeW(cesc_nmf)
load("/path/to/analysis/nmf/hnsc_nmf.RData")
hnsc_nmf = normalizeW(hnsc_nmf)
load("/path/to/analysis/nmf/laml_nmf.RData")
laml_nmf = normalizeW(laml_nmf)
load("/path/to/analysis/nmf/lihc_nmf.RData")
lihc_nmf = normalizeW(lihc_nmf)
load("/path/to/analysis/nmf/read_nmf.RData")
read_nmf = normalizeW(read_nmf)
load("/path/to/analysis/nmf/skcm_nmf.RData")
skcm_nmf = normalizeW(skcm_nmf)
load("/path/to/analysis/nmf/stad_nmf.RData")
stad_nmf = normalizeW(stad_nmf)
load("/path/to/analysis/nmf/tgct_nmf.RData")
tgct_nmf = normalizeW(tgct_nmf)
load("/path/to/analysis/nmf/ucec_nmf.RData")
ucec_nmf = normalizeW(ucec_nmf)

# clusterprofiler results
load("/path/to/analysis/enriched_terms/cp_paad.RData")
load("/path/to/analysis/enriched_terms/cp_prad.RData")
load("/path/to/analysis/enriched_terms/cp_coad.RData")
load("/path/to/analysis/enriched_terms/cp_luad.RData")
load("/path/to/analysis/enriched_terms/cp_amy.RData")
load("/path/to/analysis/enriched_terms/cp_hippo.RData")
load("/path/to/analysis/enriched_terms/cp_gbm.RData")
load("/path/to/analysis/enriched_terms/cp_pancreas.RData")
load("/path/to/analysis/enriched_terms/cp_prostate.RData")
load("/path/to/analysis/enriched_terms/cp_colon_sigmoid.RData")
load("/path/to/analysis/enriched_terms/cp_colon_transverse.RData")
load("/path/to/analysis/enriched_terms/cp_lung.RData")

load("/path/to/analysis/enriched_terms/cp_blca.RData")
load("/path/to/analysis/enriched_terms/cp_brca.RData")
load("/path/to/analysis/enriched_terms/cp_cesc.RData")
load("/path/to/analysis/enriched_terms/cp_hnsc.RData")
load("/path/to/analysis/enriched_terms/cp_laml.RData")
load("/path/to/analysis/enriched_terms/cp_lihc.RData")
load("/path/to/analysis/enriched_terms/cp_read.RData")
load("/path/to/analysis/enriched_terms/cp_skcm.RData")
load("/path/to/analysis/enriched_terms/cp_stad.RData")
load("/path/to/analysis/enriched_terms/cp_tgct.RData")
load("/path/to/analysis/enriched_terms/cp_ucec.RData")

load("/path/to/analysis/enriched_terms/cp_bladder.RData")
load("/path/to/analysis/enriched_terms/cp_blood.RData")
load("/path/to/analysis/enriched_terms/cp_breast.RData")
load("/path/to/analysis/enriched_terms/cp_ectocervix.RData")
load("/path/to/analysis/enriched_terms/cp_endocervix.RData")
load("/path/to/analysis/enriched_terms/cp_liver.RData")
load("/path/to/analysis/enriched_terms/cp_ovary.RData")
load("/path/to/analysis/enriched_terms/cp_skin_no_sun.RData")
load("/path/to/analysis/enriched_terms/cp_skin_sun.RData")
load("/path/to/analysis/enriched_terms/cp_stomach.RData")
load("/path/to/analysis/enriched_terms/cp_testis.RData")
load("/path/to/analysis/enriched_terms/cp_uterus.RData")


load("/path/to/analysis/enriched_terms/cp_result.RData")
load("/path/to/analysis/enriched_terms/cp_result_n_system.RData")

# Load all survival-related Data
blca_surv_data = read.delim("/path/to/data/survival/TCGA-BLCA.survival.tsv", header = TRUE)
brca_surv_data = read.delim("/path/to/data/survival/TCGA-BRCA.survival.tsv", header = TRUE)
cesc_surv_data = read.delim("/path/to/data/survival/TCGA-CESC.survival.tsv", header = TRUE)
coad_surv_data = read.delim("/path/to/data/survival/TCGA-COAD.survival.tsv", header = TRUE)
gbm_surv_data = read.delim("/path/to/data/survival/TCGA-GBM.survival.tsv", header = TRUE)
hnsc_surv_data = read.delim("/path/to/data/survival/TCGA-HNSC.survival.tsv", header = TRUE)
laml_surv_data = read.delim("/path/to/data/survival/TCGA-LAML.survival.tsv", header = TRUE)
lihc_surv_data = read.delim("/path/to/data/survival/TCGA-LIHC.survival.tsv", header = TRUE)
luad_surv_data = read.delim("/path/to/data/survival/TCGA-LUAD.survival.tsv", header = TRUE)
paad_surv_data = read.delim("/path/to/data/survival/TCGA-PAAD.survival.tsv", header = TRUE)
prad_surv_data = read.delim("/path/to/data/survival/TCGA-PRAD.survival.tsv", header = TRUE)
read_surv_data = read.delim("/path/to/data/survival/TCGA-READ.survival.tsv", header = TRUE)
skcm_surv_data = read.delim("/path/to/data/survival/TCGA-SKCM.survival.tsv", header = TRUE)
stad_surv_data = read.delim("/path/to/data/survival/TCGA-STAD.survival.tsv", header = TRUE)
tgct_surv_data = read.delim("/path/to/data/survival/TCGA-TGCT.survival.tsv", header = TRUE)
ucec_surv_data = read.delim("/path/to/data/survival/TCGA-UCEC.survival.tsv", header = TRUE)

# Load all Phenotype-related Data
blca_pheno = read.delim("/path/to/data/phenotype/TCGA-BLCA.GDC_phenotype.tsv", header = TRUE)
brca_pheno = read.delim("/path/to/data/phenotype/TCGA-BRCA.GDC_phenotype.tsv", header = TRUE)
cesc_pheno = read.delim("/path/to/data/phenotype/TCGA-CESC.GDC_phenotype.tsv", header = TRUE)
coad_pheno = read.delim("/path/to/data/phenotype/TCGA-COAD.GDC_phenotype.tsv", header = TRUE)
gbm_pheno = read.delim("/path/to/data/phenotype/TCGA-GBM.GDC_phenotype.tsv", header = TRUE)
hnsc_pheno = read.delim("/path/to/data/phenotype/TCGA-HNSC.GDC_phenotype.tsv", header = TRUE)
laml_pheno = read.delim("/path/to/data/phenotype/TCGA-LAML.GDC_phenotype.tsv", header = TRUE)
lihc_pheno = read.delim("/path/to/data/phenotype/TCGA-LIHC.GDC_phenotype.tsv", header = TRUE)
luad_pheno = read.delim("/path/to/data/phenotype/TCGA-LUAD.GDC_phenotype.tsv", header = TRUE)
paad_pheno = read.delim("/path/to/data/phenotype/TCGA-PAAD.GDC_phenotype.tsv", header = TRUE)
prad_pheno = read.delim("/path/to/data/phenotype/TCGA-PRAD.GDC_phenotype.tsv", header = TRUE)
read_pheno = read.delim("/path/to/data/phenotype/TCGA-READ.GDC_phenotype.tsv", header = TRUE)
skcm_pheno = read.delim("/path/to/data/phenotype/TCGA-SKCM.GDC_phenotype.tsv", header = TRUE)
stad_pheno = read.delim("/path/to/data/phenotype/TCGA-STAD.GDC_phenotype.tsv", header = TRUE)
tgct_pheno = read.delim("/path/to/data/phenotype/TCGA-TGCT.GDC_phenotype.tsv", header = TRUE)
ucec_pheno = read.delim("/path/to/data/phenotype/TCGA-UCEC.GDC_phenotype.tsv", header = TRUE)

tcga_df = data.frame(Dataset = c("BLCA", "BRCA", "CESC", "COAD", "GBM", "HNSC", "LAML", "LIHC", "LUAD", "PAAD", "PRAD", "READ",
                                 "SKCM", "STAD", "TGCT", "UCEC"),
                 seq_data = c("TCGA_BLCA_seq", "TCGA_BRCA_seq", "TCGA_CESC_seq", "TCGA_COAD_seq",  "TCGA_GBM_seq",
                              "TCGA_HNSC_seq", "TCGA_LAML_seq", "TCGA_LIHC_seq", "TCGA_LUAD_seq", "TCGA_PAAD_seq", 
                              "TCGA_PRAD_seq", "TCGA_READ_seq", "TCGA_SKCM_seq", "TCGA_STAD_seq", "TCGA_TGCT_seq", 
                              "TCGA_UCEC_seq"),
                 surv_data = c("blca_surv_data", "brca_surv_data", "cesc_surv_data", "coad_surv_data", "gbm_surv_data",
                               "hnsc_surv_data", "laml_surv_data", "lihc_surv_data", "luad_surv_data", "paad_surv_data",
                               "prad_surv_data", "read_surv_data", "skcm_surv_data", "stad_surv_data", "tgct_surv_data",
                               "ucec_surv_data"),
                pheno_data = c("blca_pheno", "brca_pheno", "cesc_pheno", "coad_pheno", "gbm_pheno", "hnsc_pheno", "laml_pheno",
                               "lihc_pheno", "luad_pheno", "paad_pheno", "prad_pheno", "read_pheno", "skcm_pheno", "stad_pheno",
                               "tgct_pheno", "ucec_pheno"),
                 nmf = c("blca_nmf", "brca_nmf", "cesc_nmf", "COAD_nmf", "GBM_nmf", "hnsc_nmf", "laml_nmf", "lihc_nmf",
                         "LUAD_nmf", "PAAD_nmf", "PRAD_nmf", "read_nmf", "skcm_nmf", "stad_nmf", "tgct_nmf", "ucec_nmf"),
                 cp = c("cp_blca", "cp_brca", "cp_cesc", "cp_coad", "cp_gbm", "cp_hnsc", "cp_laml", "cp_lihc", "cp_luad", 
                        "cp_paad", "cp_prad", "cp_read", "cp_skcm", "cp_stad", "cp_tgct", "cp_ucec"))

normal_df = data.frame(
                Dataset = c("amygdala", "bladder", "blood", "breast", "colon_sigmoid", "colon_transverse","hippocampus",
                            "liver", "lung", "ovary", "pancreas", "prostate", "skin_no_sun", "skin_sun", "stomach",
                            "testis", "uterus"),
                nmf = c("amygdala_nmf", "bladder_nmf", "blood_nmf", "breast_nmf", "colon_sigmoid_nmf", "colon_transverse_nmf",
                        "hippocampus_nmf", "liver_nmf", "lung_nmf", "ovary_nmf", "pancreas_nmf", "prostate_nmf",
                        "skin_no_sun_nmf", "skin_sun_nmf", "stomach_nmf", "testis_nmf", "uterus_nmf"),
                seq_data = c("amygdala_normal", "bladder_normal", "blood_normal", "breast_normal", "colon_sigmoid_normal",
                            "colon_transverse_normal", "hippocampus_normal", "liver_normal", "lung_normal", "ovary_normal",
                            "pancreas_normal", "prostate_normal", "skin_no_sun_normal", "skin_sun_normal", "stomach_normal",
                            "testis_normal", "uterus_normal"))
