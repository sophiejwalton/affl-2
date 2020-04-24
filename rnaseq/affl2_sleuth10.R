# load sleuth
library("sleuth")
packageVersion('sleuth')
sessionInfo()
# get all data files
base_dir <- "/Users/sophiewalton/rnaseq"
sample_id <- dir(file.path(base_dir, 'kal_quant'))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, 'kal_quant', id))
s2c <- read.table(file.path(base_dir, "sleuth", "infoHS_good_switched.txt"), header = TRUE, stringsAsFactors=FALSE)
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
so <- sleuth_prep(s2c,  ~condition + genotype + replicate + condition*genotype, target_mapping = t2g)
# to group by genes and not transcripts
#, aggregation_column = 'ens_gene')
# set the number of available cores to 4
#options(mc.cores = 4L)

# fit the full model
so <- sleuth_fit(so, ~condition + genotype + replicate + condition*genotype, 'full')
# fit the reduced model, no interactions
so <- sleuth_fit(so, ~replicate + condition + genotype , 'reduced')
# compute likelihood ratio
so <- sleuth_lrt(so, 'reduced', 'full')
lrt_results <- sleuth_results(so, 'reduced:full', 'lrt')
write.csv(lrt_results, file=file.path(base_dir, 'sleuth/lrt_results_rep.txt'), quote=FALSE, row.names=FALSE)
print(models(so))
# compute wald test on each beta
so <- sleuth_wt(so, which_beta = 'conditionyHSS')
so <- sleuth_wt(so, which_beta = 'genotypezaffl2')
so <- sleuth_wt(so, which_beta = 'conditionyHSS:genotypezaffl2')

# get results
results_GenHS <- sleuth_results(so, 'conditionyHSS:genotypezaffl2')
results_HS <- sleuth_results(so, 'conditionyHSS')
results_Gen <- sleuth_results(so, 'genotypezaffl2')

# and save them
write.csv(results_GenHS, file=file.path(base_dir, 'sleuth/results_2GenHS_good_gene_rep2.txt'), quote=FALSE, row.names=FALSE)
write.csv(results_HS, file=file.path(base_dir, 'sleuth/results_2HS_good_gene_rep2.txt'), quote=FALSE, row.names=FALSE)
write.csv(results_Gen, file=file.path(base_dir, 'sleuth/results_2Gen_good_gene_rep2.txt'), quote=FALSE, row.names=FALSE)

# launch EDA shiny app
sleuth_live(so)

units <- 'tpm'
dat_all <- NULL

for (pc_input in (1:2)){
    mat <- sleuth:::spread_abundance_by(so$obs_norm_filt, units)
    pca_calc <- prcomp(mat)
    loadings <- t(pca_calc$x)
    loadings <- pca_calc$x[, pc_input]
    loadings <- abs(loadings)
    loadings <- sort(loadings, decreasing = TRUE)
    label_names <- names(loadings)

    dat <- data.frame(pc = label_names, loadings = loadings)
    dat$pc <- factor(dat$pc, levels = unique(dat$pc))

    if (is.null(dat_all)){
        dat_all <- data.frame(dat)
        print(head(dat_all))
    }
    else {
        dat_all <- merge(dat_all, dat, by='pc')
    }
}
# save
write.csv(dat_all,
        file='pc_loadings_2.csv',
        quote=FALSE, row.names=FALSE)
'


# launch EDA shiny app

