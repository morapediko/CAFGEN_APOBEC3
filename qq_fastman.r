# import library
library(data.table) ## allow to read big file using fread function
#devtools::install_github('kaustubhad/fastman')
library(fastman)
data_sumstat<-fread('data/chromosome22_march_2024.hwe') # reading summary statistics
# identify colums contains p-value
head(data_sumstat, 2)
# dimension of data
fastqq (data_sumstat$P)

