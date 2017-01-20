library(dplyr)
require(tidyr)
tmp <- tempfile()
exacurl <- "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"
exacf <- download.file(exacurl,tmp)
exacdf <- read.table(tmp,sep="\t",header=T,stringsAsFactors = F)
exacdf <- dplyr::select(exacdf,-chr,-transcript)  %>% dplyr::rename(Gene=gene)
exacdf <- gather(exacdf,feature,value,-Gene)  %>% mutate(class="ExAC")
save(exacdf,"data/ExAC.rdata")
