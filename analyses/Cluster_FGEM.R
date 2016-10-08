#Code to (re)generate Betas.
library(FGEM)
library(dplyr)
library(BBmisc)
chunki <- 1
chunki <- as.integer(commandArgs(trailingOnly = F))
sub_all_anno_df<- readRDS("~/Dropbox/BayesianDA/FGEM_Data/sub_annotations.RDS")
datadf <- readRDS("~/Dropbox/BayesianDA/FGEM_Data/TADA.RDS")

featureinds <- unique(sub_all_anno_df$feat_ind)
mchunk <- chunk(featureinds,chunk.size = 1)
mchunk <- mchunk[[chunki]]
miniter <- min(mchunk)
maxiter <- max(mchunk)

s_anno_df <- filter(sub_all_anno_df,between(feat_ind,miniter,maxiter))

results <- group_by(s_anno_df,feature) %>%  do(cfeat_df(.,datadf,impute=F)) %>%  do(sem_df(.)) %>% ungroup()

saveRDS(results,paste0("/scratch/midway/nwknoblauch/FGEM/FGEM_Results/FGEM_",miniter,"_",maxiter,".RDS"))


