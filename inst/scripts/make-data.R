# All the original data were downloaded from https://portal.gdc.cancer.gov/
# mainly using GDC data transfer tool (https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)


# For this package, we preprocessed gene expression quantification data (including HTSeq-Counts and HTSeq-FPKM),
# miRNA expression quantification data and DNA methylation data
# (processed using Illumina Human Methylation 450 platform) for 2582 patients with cancer in five primary sites:
# adrenal_gland, lung, uterus, kidney and colorectal.
# (Note: in the paper (https://arxiv.org/abs/1708.07136), we exclude 389 cases with colorectal cancer for analysis)


# We used HPC cluster to bulk process 137438 open-access files downloaded from GDC data portal.
# It takes a long time to process 1.43 TB original data files.
# The content of `Wall` of this package is derived from 10328 files, four files for each of the 2582 patients.


# We organize gene expression, miRNA expression, and DNA methylation data into matrices,
# with each row representing a gene/miRNA/CpG site, and each column represent a sample (an aliquot indeed)
# For each of the five primary sites, we have five matrices.
# Take adrenal_gland as an example:
mat.htseq_counts.adrenal_gland # Gene Expression HTSeq Counts
mat.fpkm.adrenal_gland # Gene expression FPKM values
mat.mirnas_counts.adrenal_gland # read count for miRNAs
mat.mirnas_normalized.adrenal_gland #normalized count in reads-per-million-miRNA-mapped
mat.methy450.adrenal_gland # DNA methylation beta values (Illumina HumanMethylation450 platform)
# Note: mat.mirans_couts.adrenal_gland and mat.mirans_normalized.adrenal_gland
# are derived from the same files with file name suffix "mirnas.quantification.txt"

# The rownames of these matrices are case_id (or the first 12 letters of TCGA barcode)
# The colnames are full TCGA barcode (28 letters)
# Here are details about TCGA barcode: https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode




# We omit the script for processing raw data files to produce these matrices
# In the following will only show how we computed `Wall` in this package
# Since these matrices are large, in the following we do not actually run the script.






#***************************** How to produce `Wall` from expression/methylation matrices ***************************
# We have 2582 tumor samples from five primary sites. For each primary site, the processing procedure
# is exactly the same.

# `Wall` is a list of five lists (each corresponding to a primary sites)
names(Wall)
# In the following, we use primary site adrenal_gland for illutration.


# Using different feature selection and transformation techniques, we produce six features:
names(Wall$adrenal_gland)

#***********What are `raw.all` and `log.all`**************
# We used the same procedure to process mat.htseq_counts.adrenal_gland and mat.mirnas_counts.adrenal_gland
# In the following, we only show code for processing mat.htseq_counts.

# First we use all raw counts to compute pairwise sample distance (Euclidean distance)
dist.htseq_counts.adrenal_gland <- as.matrix(dist(t(mat.htseq_counts.adrenal_gland)))
# we calculate affinity matrix W from D using function `affinity_matrix` in ANF package (https://github.com/BeautyOfWeb/ANF)
library(ANF)
W.htseq_counts.adrenal_gland <- affinity_matrix(dist.htseq_counts.adrenal_gland, k=15)
# Now: W.htseq_counts.adrenal_gland == Wall$adrenal_gland$`raw.all`$fpkm


# Besides using raw counts, we also used log2 transformation on raw counts
dist.log2.htseq_counts.adrenal_gland <- as.matrix(dist(log2(1 + t(mat.htseq_counts.adrenal_gland))))
W.log2.htseq_counts.adrenal_gland <- affinity_matrix(dist.log2.htseq_counts.adrenal_gland, k=15)
# Now: W.htseq_counts.adrenal_gland == Wall$adrenal_gland$`log.all`$fpkm




#***********What are `raw.sel`, `log.sel` and `vst.sel`**************
# We used DESeq2 (v1.16.1) for differential expression analysis for
# both HTSeq_counts and miRNA counts.
# Since there are hundreds of samples (most of which are tumor samples), it takes quite a while
# to run DESeq2. We used default parameter settings of DESeq2 software.
library(DESeq2)
# Now suppose we combine tumor samples and normal samples into a count matrix called `cnt`
# colnames(cnt) is TCGA barcode. The first 12 letters correspond to case_id,
# while letter 14-15 are either '01' (solid tumor) or '11' (solid normal)
# Normal samples for these five primary sites were not included in this package,
# but can be downloaded from GDC data portal.
coldata = data.frame(subject=substr(colnames(cnt),1,12),type=substr(colnames(cnt),14, 15))
dds = DESeqDataSetFromMatrix(cnt, coldata, design = ~subject+type)
dds = DESeq(dds)
if (data_type == "fpkm") {
  res <- results(dds, contrast = c("type","01","11"), lfcThreshold = 1, alpha = 0.05)
  degs=rownames(res)[!is.na(res$padj)&res$padj<0.05]
}
if (data_type == "mirnas") {
  res <-results(dds, contrast = c("type","01","11"))
  degs=rownames(res)[!is.na(res$padj)&res$padj<0.1]
}

# Now we get differentially expressed genes/mirnas (i.e., `degs` in the above code)
# We can use these diffentially expressed genes/mirnas as selected features for calculating pairwise sample distance
# Again we only show an example for HTSeq-Counts. The same procedure had been applied to miRNA counts
dist.sel.htseq_counts.adrenal_gland <- as.matrix(dist(t(mat.htseq_counts.adrenal_gland[degs, ])))
W.sel.htseq_counts.adrenal_gland <- affinity_matrix(dist.sel.htseq_counts.adrenal_gland, k=15)
# Now: W.sel.htseq_counts.adrenal_gland == Wall$adrenal_gland$`raw.sel`$fpkm

dist.log2.sel.htseq_counts.adrenal_gland <- as.matrix(dist(log2(1 + t(mat.htseq_counts.adrenal_gland[degs, ]))))
W.log2.sel.htseq_counts.adrenal_gland <- affinity_matrix(dist.log2.sel.htseq_counts.adrenal_gland, k=15)
# Now: W.log2.sel.htseq_counts.adrenal_gland == Wall$adrenal_gland$`log.sel`$fpkm

# We also calculate variance stabilizing transformation of raw counts for differentially expressed
# genes/mirnas using DESeq2
vsd <- varianceStabilizingTransformation(dds[degs, ], blind = F)
W.vsd.htseq_counts.adrenal_gland <- affinityMatrix(as.matrix(dist(t(assay(vsd[, tumor.sampleIDs])))), K = 10)
# Now: W.vsd.htseq_counts.adrenal_gland == Wall$adrenal_gland$`vst.sel`$fpkm




#***********What are `normalized`*********************
# For gene and mirna counts data, we have described five features: `raw.all`, `raw.sel`, `log.all`, `log.sel` and `vst.sel`
# For genes, we also have FPKM values;
# For miRNAs, we also have normalized count in reads-per-million-miRNA-mapped
# For DNA methylation data, we treat beta values as normalized values, too.
# We process these "normalized" features as follows.
# We use DNA methylation beta values as an examples. FPKM and miRNA normalized expression had been treated the same.
dist.methy450.adrenal_gland = 1 - cor(mat.methy450.adrenal_gland)
W.methy450.adrenal_gland <- affinity_matrix(dist.methy450.adrenal_gland, 10)
# Now: W.methy450.adrenal_gland == Wall$adrenal_gland$normalized$methy450



# Note: since DNA methylation data do not contain raw counts data,
# we actually do not have corresponding `raw.all`, `raw.sel`, `log.all`, `log.sel` and `vst.sel` features
# But for the convenience of using ANF, we also added `methy450` matrix to the above five features.
# That's to say:
# Wall$adrenal_gland$`raw.all`$methy450 == Wall$adrenal_gland$`raw.sel`$methy450 ==
# Wall$adrenal_gland$`log.all`$methy450 == Wall$adrenal_gland$`log.sel`$methy450 ==
# Wall$adrenal_gland$`vst.sel`$methy450 == Wall$adrenal_gland$`normalized`$methy450





# ******************* File UUIDs of 10328 files **************
# For the convenience of users, we also provide file UUIDs of these 10328 files
# Users can download raw files from GDC data portal.
# Note: miRNA expression data was from GDC data release 7.0, while GDC data release 8.0 has updated updated miRNA files

# The following code shows how we generate the data.frame `meta`
# It cannot be run by users since "file2aliquot.all.RData" has not been made publicly available, which maps file ID to TCGA barcode (i.e., aliquot ID).
data("Wall")
data("project_ids")
# file2aliquot.all.RData is stored locally
load("F:/TCGA/processed_RData/file2aliquot.all.RData")

gen_metadata <- function(aliquotIDs, fileUUIDs, meta, data_type_name, project_ids) {
  for (name in names(aliquotIDs)) {
    file_ids <- names(fileUUIDs[[name]])[match(aliquotIDs[[name]],
                                               fileUUIDs[[name]])]
    meta <- rbind(meta, data.frame(file_id=file_ids,
                                   aliquot_id=aliquotIDs[[name]],
                                   data_type=data_type_name,
                                   case_id=substr(aliquotIDs[[name]], 1, 12),
                                   primary_site=name,
                                   project_id=project_ids[substr(aliquotIDs[[name]], 1, 12)]))
  }
  return(meta)
}

meta <- data.frame()

# HTSeq-Counts:
aliquotIDs <- sapply(Wall, function(x) colnames(x$raw.all$fpkm))
fileUUIDs <- sapply(column2aliquot.all,
                    function(x) x$`htseq.counts.gz$`)
data_type_name <- "Gene Expression HTSeq-Counts"
meta <- gen_metadata(aliquotIDs, fileUUIDs, meta, data_type_name, project_ids)

# FPKM:
aliquotIDs <- sapply(Wall, function(x) colnames(x$normalized$fpkm))
fileUUIDs <- sapply(column2aliquot.all,
                    function(x) x$`FPKM.txt.gz$`)
data_type_name <- "Gene Expression FPKM Values"
meta <- gen_metadata(aliquotIDs, fileUUIDs, meta, data_type_name, project_ids)

# miRNA expression:
aliquotIDs <- sapply(Wall, function(x) colnames(x$raw.all$mirnas))
fileUUIDs <- sapply(column2aliquot.all,
                    function(x) x$`mirnas.quantification.txt$`)
data_type_name <- "miRNA Expression Quantifications"
meta <- gen_metadata(aliquotIDs, fileUUIDs, meta, data_type_name, project_ids)

# DNA methylation beta values
aliquotIDs <- sapply(Wall, function(x) colnames(x$raw.all$methy450))
fileUUIDs <- sapply(column2aliquot.all,
                    function(x) x$HumanMethylation450)
data_type_name <- "DNA Methylation Beta Values"
meta <- gen_metadata(aliquotIDs, fileUUIDs, meta, data_type_name, project_ids)

write.table(meta, file = "inst/extdata/fileUUIDs.csv", sep = ',', quote = FALSE, row.names = FALSE)






#***************************** What are `project_ids` and `surv.plot`********************

#`project_ids` is a named character mapping case_id to TCGA project ID (one project corresponds to one disease type),
# which can be directly derived from https://portal.gdc.cancer.gov/repository?facetTab=cases&searchTableTab=cases
head(project_ids)



# `surv.plot` (overall survival plot data) were downloaded from
# https://portal.gdc.cancer.gov/exploration?searchTableTab=genes
# It is used for patient suvival analysis and plotting survival plot.
