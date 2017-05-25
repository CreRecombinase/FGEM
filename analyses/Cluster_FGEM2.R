#Code to (re)generate Betas.
library(FGEM)
library(dplyr)
library(BBmisc)

#chunki <- 1

chunki <- as.integer(commandArgs(trailingOnly = T))
cat(chunki,"\n")
GOdf <- distinct(BPGOdf,feature) %>% mutate(feat_ind=1:n()) %>% inner_join(BPGOdf)
datadf <- BFdata
featureinds <- unique(GOdf$feat_ind)
cat("There are ",length(featureinds),"features\n")

mchunk <- chunk(featureinds,chunk.size = 100)
cat("There are ",length(mchunk)," chunks of size 100\n")
mchunk <- mchunk[[chunki]]
miniter <- min(mchunk)
maxiter <- max(mchunk)

s_anno_df <- filter(GOdf,between(feat_ind,miniter,maxiter))

results <- select(s_anno_df,-feat_ind) %>% group_by(feature) %>%  do(sem_df(cfeat_df(.,datadf,impute=F))) %>% ungroup()

saveRDS(results,paste0("~/Dropbox/BayesianDA/FGEM_Data/Exac_res_2.RDS"))
saveRDS(results,paste0("/scratch/t.cri.nknoblauch/FGEM_Results/FGEM_",miniter,"_",maxiter,".RDS"))



