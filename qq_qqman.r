# import library
library(data.table) ## allow to read big file using fread function
library(qqman)
data_sumstat<-fread('data/chromosome22_march_2024.hwe') # reading summary statistics
qq(data_sumstat$P) ## implemented same line than before.
