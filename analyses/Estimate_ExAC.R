library(FGEM)
library(dplyr)


annodf <- exacdf
feature_df <- distinct(annodf,feature) %>% mutate(feat_ind=1:n()) %>% inner_join(annodf)
datadf <- BFdata
featureinds <- unique(feature_df$feat_ind)
cat("There are ",length(featureinds),"features\n")

null_feats <- c("bp","cds_end","cds_start")
test_feats <- c("mu_mis")
feats <- c(null_feats,test_feats)
all_feats <- distinct(feature_df,feature)
tfeature_df <- filter(feature_df,feature %in% feats)

data_df <- cfeat_df(tfeature_df,datadf = datadf)

nresults <- em_df(data_df,null_features=c("bp"),test_features=c("mu_mis"))
test_feats <- filter(feature_df,!feature%in% null_feats) %>% distinct(feature)

sall_results <- group_by(test_feats,feature) %>% do(em_df(cfeat_df(filter(feature_df,feature %in% c(.$feature,null_feats)),datadf),
                                                         prior_mean = 0.02,scale=T,null_feature = null_feats)) %>% ungroup()
all_results


results <- select(tfeature_df,-feat_ind) %>% group_by(feature) %>%  do(sem_df(cfeat_df(.,datadf,impute=F))) %>% ungroup()
results <- mutate(results,qval=p.adjust(pval)) %>% arrange(qval)

saveRDS(results,paste0("/scratch/t.cri.nknoblauch/FGEM_Results/FGEM_",miniter,"_",maxiter,".RDS"))


