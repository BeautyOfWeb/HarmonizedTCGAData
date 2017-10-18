meta <- data.frame(
  Title = "Processed harmonized TCGA data of five cancer types for patient clustering",
  Description = paste0("Pre-computed affinity matrices ",
                       "of 2582 patients from nine TCGA projects ",
                       "using gene expression, miRNA expression, ",
                       "and DNA methylation data"),
  BiocVersion = "3.6",
  Genome = "GRCh38",
  SourceType = "TXT, TSV",
  SourceUrl = "https://portal.gdc.cancer.gov/repository",
  SourceVersion = "June 29 2017",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "Harmonized Cancer Datasets Genomic Data Commons Data Portal",
  Maintainer = "Tianle Ma <tianlema@buffalo.edu>",
  RDataClass = "List",
  DispatchClass = "Rda",
  ResourceName = "Wall.rda",
  RDataPath = "HarmonizedTCGAData/Wall.rda"
)

meta <- rbind(meta, data.frame(
  Title = "TCGA project IDs for 14551 cases",
  Description = paste0("Map case ID to project ID",
                       " (each project ID correspond to a
                       disease type)"),
  BiocVersion = "3.6",
  Genome = "GRCh38",
  SourceType = "JSON",
  SourceUrl = "https://portal.gdc.cancer.gov/repository?searchTableTab=cases",
  SourceVersion = "June 29 2017",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "Harmonized Cancer Datasets Genomic Data Commons Data Portal",
  Maintainer = "Tianle Ma <tianlema@buffalo.edu>",
  RDataClass = "Character",
  DispatchClass = "Rda",
  ResourceName = "project_ids.rda",
  RDataPath = "HarmonizedTCGAData/project_ids.rda"
))

meta <- rbind(meta, data.frame(
  Title = "Survival Data for 12899 TCGA Cases",
  Description = paste0("Data used for overall survival plot ",
                       "of 12899 TCGA cases"),
  BiocVersion = "3.6",
  Genome = "GRCh38",
  SourceType = "JSON",
  SourceUrl = "https://portal.gdc.cancer.gov/exploration?searchTableTab=genes",
  SourceVersion = "June 29 2017",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "Harmonized Cancer Datasets Genomic Data Commons Data Portal",
  Maintainer = "Tianle Ma <tianlema@buffalo.edu>",
  RDataClass = "data.frame",
  DispatchClass = "Rda",
  ResourceName = "surv.plot.rda",
  RDataPath = "HarmonizedTCGAData/surv.plot.rda"
))

## Not run:
## Write the data out and put in the inst/extdata directory.
write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)

## Test the validity of metadata.csv with readMetadataCsv():
library(AnnotationHubData)
readMetadataFromCsv("D:/github/HarmonizedTCGAData")

## End(Not run)
