#Code to (re)generate Betas.
library(FGEM)
library(dplyr)
library(BBmisc)
library(methods)
#chunki <- 1
chunki <- as.integer(commandArgs(trailingOnly = T))
cat(chunki,"\n")
sub_all_anno_df<- readRDS("/scratch/midway/nwknoblauch/FGEM/sub_annotations.RDS")
datadf <- readRDS("/scratch/midway/nwknoblauch/FGEM/TADA.RDS")

featureinds <- unique(sub_all_anno_df$feat_ind)
cat("There are ",length(featureinds),"features\n")

mchunk <- chunk(featureinds,chunk.size = 100)
cat("There are ",length(mchunk)," chunks of size 20\n")
mchunk <- mchunk[[chunki]]
miniter <- min(mchunk)
maxiter <- max(mchunk)

s_anno_df <- filter(sub_all_anno_df,between(feat_ind,miniter,maxiter))

results <- group_by(s_anno_df,feature) %>%  do(cfeat_df(.,datadf,impute=F)) %>%  do(sem_df(.)) %>% ungroup()

saveRDS(results,paste0("/scratch/midway/nwknoblauch/FGEM/FGEM_Results/FGEM_",miniter,"_",maxiter,".RDS"))


