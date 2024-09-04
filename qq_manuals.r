# import library
library(data.table) ## allow to read big file using fread function
data_sumstat<-fread('data/chromosome22_march_2024.hwe') # reading summary statistics
# identify colums contains p-value
head(data_sumstat, 2)
# dimension of data
dim(data_sumstat)
# extraction of column contains p_value (P for plink or p_wald for gemma)
pvector<-data_sumstat$P
# clean na value, null value and allow finite vale
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0]
# sort and transform p -value 
o = -log10(sort(pvector, decreasing = FALSE))
# computed value expected, pppoints probability points for the continuous sample quantile 
e = -log10(ppoints(length(pvector)))
# plot of expexted points and observed in y-axis
plot(e,o, xlab='Expected', ylab='Observed')
# add red line
abline(0, 1, col = "red")

