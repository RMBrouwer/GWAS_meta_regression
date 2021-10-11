library(data.table)
library(foreach)
library(metafor)
library(iterators)

#s("~/Documents/Meta_regression")

cmdargs = commandArgs(trailingOnly=TRUE)

VAR=cmdargs[1] # phenotype
CHR=cmdargs[2] # chromosome
TYPE=cmdargs[3] # type of meta-age-analysis - n_cm with n is the number of beta's in the analysis (2 linear or 3 quadratic) and m is the number of variables that are tested (1/2 for the age-beta's versus 2/3 including the intercept), see description below.

# this file should contain columns "cohort", "mean_age" and "N"
metadata <- read.table(sprintf("summary_measures_%s.txt",VAR),header=T)
inlist <- metadata[,"cohort"]

# we do this per chromosome to limit usage of memory

# assembling the datafile
inputfile=sprintf("gunzip -cq cohort_input_meta_regression/aligned_%s/%s/%s_meta_regression_%s_chr%s.gz",VAR,inlist[1],inlist[1],VAR,CHR)
data <- fread(cmd=inputfile,header=F,colClasses=c(NA,NA,NA,NA,NA,"NULL","NULL","NULL"))
names(data) <- c("markername","allele1","allele2",sprintf("effect_%s",inlist[1]),sprintf("se_%s",inlist[1]))
print(sprintf("datafile %s read, number of snps is now %i",sprintf("%s_meta_regression_%s_chr%s",inlist[1],VAR,CHR),dim(data)[1]))
for (C in inlist[2:length(inlist)])
{
	datatmp <- fread(cmd=sprintf("gunzip -cq cohort_input_meta_regression/aligned_%s/%s/%s_meta_regression_%s_chr%s.gz",VAR,C,C,VAR,CHR),header=F,colClasses=c(NA,NA,NA,NA,NA,"NULL","NULL","NULL"))
	names(datatmp) <- c("markername","allele1","allele2",sprintf("effect_%s",C),sprintf("se_%s",C))
	data <- merge(datatmp,data,by=c("markername","allele1","allele2"),all.x=T,all.y=T)
	print(sprintf("datafile %s read, number of SNPs is now %i",sprintf("%s_meta_regression_%s_chr%s",C,VAR,CHR),dim(data)[1]))
}

# storing names of effects, standard error of effects, and independent variables in the meta-regression
effects <- sprintf("effect_%s",inlist)
ses <- sprintf("se_%s",inlist)
ages <- metadata[,"mean_age"]
ages2 <- ages^2

# filter SNPS on the number of cohorts (>50%) and the number of subjects (>50%)
N <- metadata[,"N"]
maxN <- sum(metadata$N,na.rm=T) 
maxC <- length(inlist)

# this function works on a vector X where X is a row of the data frame
check <- function (X){
return((length(which(is.finite(t(X[,..effects])))) > maxC/2) & (sum(N[which(is.finite(t(X[,..effects])))]) > maxN/2))
}

print("--- selecting snps ---")
# Use this to get estimates only in the SNPS that survive the filtering on number of cohorts/subjects 
inclusion_filter <- foreach(snp=iter(data, by='row')) %dopar% {check(snp)}
data <- data[unlist(inclusion_filter)]
rm(inclusion_filter)

# Define functions to compute various versions of the meta-regression
# These are named after the number of estimates in the meta-regression and the number of covariates to be tested.
# X is a row of the data containing all information for a certain SNP 

# Testing a quadratic function of age, testing for an age effect. Null hypothesis: there might be an effect of the SNP, but it is not dependent on age.
metareg3_c2 <- function(X) {
	result <- try(rma.uni(t(X[,..effects]),sei=t(X[,..ses]),mods=cbind(ages,ages2),method="FE",btt=c(2:3)))
}

# Testing a linear function of age, testing for an age effect. Null hypothesis: there might be an effect of the SNP, but it is not dependent on age.
metareg2_c1 <- function(X) {
	result <- try(rma.uni(t(X[,..effects]),sei=t(X[,..ses]),mods=cbind(ages),method="FE",btt=c(2)))
}

# functions to extract the necessary output from the meta-regressions

getQMp <- function(lst){try(lst$QMp)}
getQEp <- function(lst){try(lst$QEp)}
getQM <- function(lst){try(lst$QM)}
getest <- function(lst){try(t(lst$beta))}
getse <- function(lst){try(t(lst$se))}

print("starting analysis")

if (TYPE == "3_c2")
{
results2 <- foreach(snp=iter(data, by='row')) %dopar% { metareg3_c2(snp) }
print(sprintf("results done for chromosome %s",CHR))
data[,"QM2"] <- unlist(lapply(results2,getQM))
data[,"QMp2"] <- unlist(lapply(results2,getQMp))
data[,"QEp2"] <- unlist(lapply(results2,getQEp))
temp <- do.call(rbind,lapply(results2,getest))
colnames(temp) <- c("int2","beta_age2","beta_agesq2")
data <- cbind(data,temp)
temp <- do.call(rbind,lapply(results2,getse))
colnames(temp) <- c("se_int2","se_beta_age2","se_beta_agesq2")
data <- cbind(data,temp)

rm(temp,results2)
}


if (TYPE == "2_c1")
{
results1 <- foreach(snp=iter(data, by='row')) %dopar% {metareg2_c1(snp) }
print(sprintf("results 1 done for chromosome %s",CHR))
data[,"QMp1"] <- unlist(lapply(results1,getQMp))
data[,"QEp1"] <- unlist(lapply(results1,getQEp))
data[,"QM1"] <- unlist(lapply(results1,getQM))
temp <- do.call(rbind,lapply(results1,getest))
colnames(temp) <- c("int2","beta_age2")
data <- cbind(data,temp)
temp <- do.call(rbind,lapply(results1,getse))
colnames(temp) <- c("se_int2","se_beta_age2")
data <- cbind(data,temp)
rm (temp,results1)
} 


save.image(file=sprintf("Meta_regression_%s_chr%s_%s.Rdata",VAR,CHR,TYPE))
rm(data)

quit()

