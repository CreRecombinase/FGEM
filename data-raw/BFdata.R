library(FGEM)
library(readr)
library(dplyr)
datalink <- "https://www.dropbox.com/s/7k62nly5jzj2j00/TADA_ASC_SSC_results_Dec23.csv?dl=1"
datadf <- read_delim(datalink,delim=",",col_names=T)
datadf <- dplyr::select(datadf,Gene,BF,qvalue,pval.TADA)
BFdata <- datadf
save(datadf,"data/BFdata.rdata")
