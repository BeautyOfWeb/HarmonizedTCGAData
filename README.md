# Processed high-quality harmonized TCGA data of five cancer types

**If you use the data from this package in published research, please cite:**

> Tianle Ma, Aidong Zhang,
> Integrate Multi-omic Data Using Affinity Network Fusion (ANF) for Cancer Patient Clustering, 
> https://arxiv.org/abs/1708.07136

This package contains three R objects: `Wall`, `project_ids` and `surv.plot`:

`Wall` contains lists inside list. In fact, `Wall` a list (five cancer type) of list (six feature normalization types: `raw.all`, `raw.sel`, `log.all`, `log.sel`,  `vst.sel`, `normalized`) of list (three feature spaces or views: `fpkm`, `mirna`, and `methy450`) of matrices. The rownames of each matrix is the submitter_id (can be seen as a patient id), and the column names of each matrix is the aliquot ID (which contains the submitter_id as prefix). Based on these aliquot ID, users can download original data from https://portal.gdc.cancer.gov/repository .

`project_ids` is a named character vector, that maps the submitter_id (represent a patient) to project_id (one-to-one correspond to disease type). This is used for evaluating clustering results, such as calculating NMI and Adjusted Rand Index (ARI).

`surv.plot` is a data.frame containing patient survival data for survival analysis, providing an "indirect" way to evaluate clustering results.

See paper https://arxiv.org/abs/1708.07136 for more explanation.
