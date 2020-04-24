# load sleuth
library("sleuth")
# get all data files
base_dir <- "/Users/sophiewalton/rnaseq"
sample_id <- dir(file.path(base_dir, 'WT'))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, 'WT', id))
s2c <- read.table(file.path(base_dir, "sleuth", "infoWT.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
# !!must manually make sure s2c and kal_dirs are consistent before proceeding
print(s2c)
# get gene names
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "celegans_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
"external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# make sleuth object
so <- sleuth_prep(s2c, target_mapping = t2g)
# to group by genes and not transcripts
#, aggregation_column = 'ens_gene')
# set the number of available cores to 4
#options(mc.cores = 4L)


# fit the full model, no interactions
so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')
# compute likelihood ratio
so <- sleuth_lrt(so, 'reduced', 'full')

# compute wald test on each beta
so <- sleuth_wt(so, which_beta = 'conditionyHSS')
#so <- sleuth_wt(so, which_beta = 'genotypeWT')
#so <- sleuth_wt(so, which_beta = 'genotypeWT:conditionNO_HS')

# get results
#results_GenHS <- sleuth_results(so, 'genotypeWT:conditionNO_HS')
results_HS <- sleuth_results(so, 'conditionyHSS')
#results_Gen <- sleuth_results(so, 'genotypeWT')

# and save them
#write.csv(results_GenHS, file=file.path(base_dir, 'sleuth/results_GenHS_good.txt'), quote=FALSE, row.names=FALSE)
write.csv(results_HS, file=file.path(base_dir, 'sleuth/results_WTHS_goodkk.txt'), quote=FALSE, row.names=FALSE)
#write.csv(results_Gen, file=file.path(base_dir, 'sleuth/results_Gen_good.txt'), quote=FALSE, row.names=FALSE)

# launch EDA shiny app
sleuth_live(so)

# save sleuth object
#saveRDS(so, file=file.path(base_dir, 'sleuth/so.rds'))

# to load sleuth object
library("sleuth")
base_dir <- "/Users/sophiewalton/rnaseq"
so <- readRDS(file.path(base_dir, 'sleuth/so.rds'))
