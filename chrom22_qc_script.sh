#!/bin/bash
#this script executes quality checks on the WES data, following best practices as captured from articles by:
# 1. Marees AT, de Kluiver H, Stringer S, Vorspan F, Curis E, Marie-Claire C, Derks EM.
	# A tutorial on conducting genome-wide association studies: Quality control and statistical analysis.
	# Int J Methods Psychiatr Res. 2018 Jun;27(2):e1608. doi: 10.1002/mpr.1608. Epub 2018 Feb 27. PMID: 29484742; PMCID: PMC6001694
# 2. Anderson, C., Pettersson, F., Clarke, G. et al. 
	#Data quality control in genetic case-control association studies.
	#Nat Protoc 5, 1564–1573 (2010). https://doi.org/10.1038/nprot.2010.116
# 3. Truong, V. Q., Woerner, J. A., Cherlin, T. A., Bradford, Y., Lucas, A. M., Okeh, C. C., 
	#Shivakumar, M. K., Hui, D. H., Kumar, R., Pividori, M., Jones, S. C., Bossa, A. C., Turner, S. D., Ritchie, M. D., & Verma, S. S.
	#(2022). Quality Control Procedures for Genome-Wide Association Studies. Current protocols, 2(11), e603.
	#https://doi.org/10.1002/cpz1.603
# 4. Turner S, Armstrong LL, Bradford Y, Carlson CS, Crawford DC, Crenshaw AT, de Andrade M, Doheny KF, Haines JL, Hayes G, Jarvik G,
	# Jiang L, Kullo IJ, Li R, Ling H, Manolio TA, Matsumoto M, McCarty CA, McDavid AN, Mirel DB, Paschall JE, Pugh EW,
	# Rasmussen LV, Wilke RA, Zuvich RL, Ritchie MD. Quality control procedures for genome-wide association studies. 
	#Curr Protoc Hum Genet. 2011 Jan;Chapter 1:Unit1.19.
	# doi: 10.1002/0471142905.hg0119s68. PMID: 21234875; PMCID: PMC3066182.

#### THERE ARE SEVERAL STEPS INVOLVED #### 

#### OVERALL INDIVIDUAL AND SNP MISSINGNESS
	### --missing produces a report of missing  genotype data per individual and marker
plink --bfile cleaned_out_feb_2024 --missing

###### NEED TO RUN AN R-SCRIPT THAT WILL PLOT HISTOGRAMS OF SNP AND INDIVIDUAL MISSINGNESS #####
	###THIS SCRIPT IS CALLED hist_miss.R  and will be executed by command Rscript --no-save hist_miss.R #####

####################################################################################################################################
###### SEX CHECK AND CLEAN UP 
		### a. This command uses Plink to check the sex chromosomes of each individual in the data set cleaned_out_feb_2024.
		# If any individuals have discrepancies between their reported sex and their sex chromosomes,
		# they will be flagged with the label "PROBLEM" in the cleaned_out_feb_2024.sexcheck file.
plink --bfile cleaned_out_feb_2024 --check-sex --make-bed --out cleaned_out_feb_2024_sex_check

		#### NEED TO RUN AN R-SCRIPT FOR PLOTTING HISTOGRAM PROPORTIONS OF SEX ####
			#### THIS SCRIPT IS CALLED gender_check.R and will be executed by command  Rscript --no-save gender_check.R ####
		
		### b. This command uses the grep utility to extract any lines in the cleaned_out_feb_2024.sexcheck file 
		#that contain the label "PROBLEM", and saves them to a new file called cleaned_out_feb_2024_sexprobs.txt.
grep PROBLEM cleaned_out_feb_2024_sex_check.sexcheck > cleaned_out_feb_2024_sexprobs.txt

		### c. This command uses Plink to exclude any individuals with sex chromosome issues from the data set cleaned_out_feb_2024, 
		#by using the --exclude option and specifying the file cleaned_out_feb_2024_sexprobs.txt.
		# Plink creates a new data set called cleaned_out_feb_2024 without the individuals with sex chromosome issues.
plink --bfile cleaned_out_feb_2024_sex_check --remove cleaned_out_feb_2024_sexprobs.txt --make-bed --out cleaned_out_feb_2024_sex_pass

##########################################################################################################################################


######## 2. AWK command for filtering-out sex chromosomes. The study is not looking for sex related disorders, therefore, sex chromosomes 
	# are not neccessary from this point onwards.
	#In this case, the command is checking if the first column of the HapMap_3_r3_6.bim file is between 1 and 22 (inclusive). 
	#If the first column meets this condition,
	#the corresponding value in the second column (SNP ID) is printed and saved to a new file called snp_1_22.txt.
 awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' cleaned_out_feb_2024_sex_pass.bim > snp_1_22.txt
plink  --bfile cleaned_out_feb_2024_sex_pass --extract snp_1_22.txt --make-bed --out cleaned_out_feb_2024_autoSNPs


#############################################################################################################################################
######### 3. This command filters the genotype data file cleaned_out_feb_2024_autoSNPs based on genotyping call rate (--geno 0.10 option),
	# and creates a new file cleaned_out_feb_2024_clean_geno  that only includes SNPs with a genotyping call rate greater than
	# or equal to 90%.
plink --bfile cleaned_out_feb_2024_autoSNPs --geno 0.10 --make-bed --out cleaned_out_feb_2024_clean_geno

#############################################################################################################################################

######## 4. This command filters the genotype data file cleaned_out_feb_2024_clean_geno based on individual call rate (--mind 0.10 option),
	# and creates a new file cleaned_out_feb_2024_clean_mind_10 that only includes individuals with a call rate greater than 
	#or equal to 90%.
plink --bfile cleaned_out_feb_2024_clean_geno --mind 0.10 --make-bed  --out cleaned_out_feb_2024_clean_mind_10

####### 5.#A# This command will calculate the minor allele frequency distribution and also plot a histogram of the discribution #####
plink --bfile cleaned_out_feb_2024_clean_mind_10 --freq --out MAF_check
	######  #B# The R script called MAF_check.R will plot histogram of MAF distribution using command: Rscript --no-save MAF_check.R


############################################################################################################################################
######## 6. This command filters the genotype data file cleaned_out_feb_2024_clean_mind_10 based on minor allele frequency (--maf 0.05 option),
	# and creates a new file MAF_greater_5 that only includes SNPs with a minor allele frequency greater than or equal to 5%.
plink --bfile cleaned_out_feb_2024_clean_mind_10 --maf 0.05 --make-bed  --out MAF_greater_5

#############################################################################################################################################
####### 7. This command extracts SNPs with a minor allele frequency greater than or equal to 5% from the file named MAF_greater_5.bim
	 # and directs the output to a text file named MAF_greater_5.txt
cat MAF_greater_5.bim > MAF_greater_5.txt

####### 8. This command filters the genotype data file cleaned_out_feb_2024_clean_mind_10 based on the list of 
	#SNPs excluded from the previous command (--exclude MAF_greater_5.txt option),
	# and creates a new file MAF_less_5 that only includes the remaining SNPs of MAF less than 5%.
plink --bfile cleaned_out_feb_2024_clean_mind_10 --exclude MAF_greater_5.txt --make-bed --out MAF_less_5


############################################################################################################################################
###### 9. This command filters the genotype data file MAF_greater_5 based on genotyping call rate (--geno 0.05 option), 
	#and creates a new file MAF_greater_5_clean that only includes SNPs with a genotyping call rate greater than or equal to 95%.
plink --bfile MAF_greater_5 --geno 0.05 --make-bed --out MAF_greater_5_clean


#########################################################################################################################################
####### 10.This command filters the genotype data file MAF_less_5 based on genotyping call rate (--geno 0.01 option),
	# and creates a new file MAF_less_5_clean that only includes SNPs with a genotyping call rate greater than or equal to 99%. 
plink --bfile MAF_less_5 --geno 0.01 --make-bed --out MAF_less_5_clean


##########################################################################################################################################
####### 11. This command merges the genotype data files MAF_greater_5_clean and MAF_less_5_clean into a single file cafgen_WES_MAF_clean.
plink --bfile MAF_greater_5_clean --bmerge MAF_less_5_clean --make-bed --out cafgen_WES_MAF_clean


##########################################################################################################################################
####### 12.This script filters the cafgen_WES_MAF_clean data set based on the minimum individual call rate. 
	#The --mind option specifies the minimum call rate required for an individual to be included in the data set,
	# and the --make-bed option creates a new Plink binary format file with the filtered data. 
	#The minimum individual call rate is set to 2%,meaning that any individual with more than 98% missing data will be excluded.
plink --bfile cafgen_WES_MAF_clean --mind 0.02 --make-bed --out cafgen_WES_MAF_clean2


#########################################################################################################################################
###### 13. This script will calculate the MAF distribution of cafgen_WES_MAF_clean2 and plot out a histogram of the distribution
	### this is to show the variance between distributions before and after MAF prunning####
plink --bfile cafgen_WES_MAF_clean2 --freq --out MAF_check_clean2
  #### R-script is called  MAF_check_clean2.R and will be run using command   Rscript --no-save MAF_check_clean2.R

#########################################################################################################################################
###### 14.This script checks for distribution of Hardy-Weinberg equilibrium (HWE) in SNPs
plink --bfile cafgen_WES_MAF_clean2 --hardy

###### 15. This script selects SNPs with HWE p-value below 0.0001

awk '{ if ($9 <0.0001) print $0 }' plink.hwe>plinkzoomhwe.hwe
			#Rscript --no-save hwe.R

#############################################################################################################################################

# Filtering HWE for cases and controls using the same threshold (1e-4)
plink --bfile cafgen_WES_MAF_clean2 --hwe 1e-4 --make-bed --out cafgen_WES_hwe_filter_step1

# Use the same HWE threshold for cases
plink --bfile cafgen_WES_hwe_filter_step1 --hwe 1e-4 include-nonctrl --make-bed --out cafgen_WES_hwe_clean

#############################################################################################################################################
#STEPS TO REMOVE HETEROZYGOSITY RATE OUTLIERS
####The command --indep-pairwise will prune the SNPs
#The parameters ‘1000 100 0.2' refers to the following:
#1000: window size, is the number of SNPs in a window
#100: number of SNPs to shift the window at each step
#0.2: multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously
#The window size and the number of SNPs to shift the window at each step are used to divide the genome into non-overlapping windows
#Multiple correlation coefficient is used to determine the strength of the association between a SNP and the other SNPs in the window

plink --bfile cafgen_WES_hwe_clean --maf 0.05 --indep-pairwise 1000 100 0.2 --out pruned_data
## Prunned dataset
plink --bfile cafgen_WES_hwe_clean --extract pruned_data.prune.in --het --out R_check 
### visualization of heterozygosity rate distribution
Rscript --no-save check_heterozygosity_rate.R
#list of individuals who deviate >3 standard deviations from the heterozygosity rate mean
Rscript --no-save heterozygosity_outliers_list.R
#Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
#Remove heterozygosity rate outliers
plink --bfile cafgen_WES_hwe_clean --remove het_fail_ind.txt --make-bed --out pruned_data_final

#############################################################################################################################################
#############################################################################################################################################
### CRYPTIC RELATEDNESS ####
	# Execute the command to use pruned_data_final  dataset and extract the variants specified in pruned_data.prune.in. 
	#Compute the IBD estimates using --genome’ flag by filtering the estimates above 0.2 and then output the results to files with 
	#prefix "pihat_min0.2"
# Check for relationships between individuals with a pihat > 0.2.
plink  --bfile pruned_data_final --extract pruned_data.prune.in --genome  --out clean_relation 

#####################################################################################################################
 
#  we are looking to remove individuals with pihat >0.2 with this R script 
	#In RStudio:
#ibd <- read.table("clean_relation.genome", header=T)
#exclusions = ibd[ ibd$PI_HAT > 0.1, c('FID2','IID2')]
#write.table( exclusions, file="related_samples.txt", col.names = F, row.names = F, quote = F )

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.1
plink --bfile pruned_data_final --remove related_samples.txt --make-bed --out cafgen_WES_qc_passed

echo "quality checks performed according to best gwas cleaning practices"
