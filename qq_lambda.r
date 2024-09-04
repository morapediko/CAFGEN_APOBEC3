# import library
library(data.table) ## allow to read big file using fread function
data_sumstat<-fread('data/chromosome22_march_2024.hwe') # reading summary statistics

chisq <- qchisq(1-data_sumstat$P[!is.na(data_sumstat$P)],1) # transform p in chisq using density
lambda<-median(chisq)/qchisq(0.5,1) # computed lambda

print(lambda)

