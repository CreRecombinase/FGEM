#Code to (re)generate Betas.
library(FGEM)
library(dplyr)
library(BBmisc)
chunki <- 1
# chunki <- as.integer(commandArgs(trailingOnly = F))
sub_all_anno_df<- readRDS("~/Dropbox/BayesianDA/FGEM_Data/sub_annotations.RDS")
sub_all_anno_df <- filter(sub_all_anno_df,grepl("GO",feature))
datadf <- readRDS("~/Dropbox/BayesianDA/FGEM_Data/TADA.RDS")

featureinds <- unique(sub_all_anno_df$feat_ind)
tfi <- filter(sub_all_anno_df, feat_ind==2)
mchunk <- chunk(featureinds,chunk.size = 50)
for(chunki in 90:length(mchunk)){
  mchunkl <- mchunk[[chunki]]
  miniter <- min(mchunkl)
  maxiter <- max(mchunkl)
  s_anno_df <- filter(sub_all_anno_df,between(feat_ind,miniter,maxiter))
  f_results <- group_by(s_anno_df,feature) %>% do(cfeat_df(.,datadf,impute=F))
  %>% do(fisher_comp(.))
  #results <- group_by(s_anno_df,feature) %>%  do(cfeat_df(.,datadf,impute=F)) %>%  do(sem_df(.)) %>% ungroup()
  saveRDS(f_results,paste0("~/Dropbox/BayesianDA/FGEM_Data/FGEM_chisq_",miniter,"_",maxiter,".RDS"))
  gc()
}


