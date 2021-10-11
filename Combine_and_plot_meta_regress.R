# script to create Manhattan and QQ plots, and to create files for LD score regression, adding rsIDs.
library(stringr)	
library(qqman)

doplot=1

cmdargs = commandArgs(trailingOnly=T);	
var2=cmdargs[1]	
analysis=cmdargs[2]

#setwd("~/Documents/Meta_regression")
if (analysis == "2_c1")
{
Qp="QMp1"
effect_names <- c("int2","se_int2","beta_age2","se_beta_age2")
nice_effect_names <- c("intercept","se_intercept","beta_age","se_beta_age")
Qe="QEp1"
}
if (analysis == "3_c2")
{
Qp="QMp2"
effect_names <- c("int2","se_int2","beta_age2","se_beta_age2", "beta_agesq2","se_beta_agesq2")
nice_effect_names <- c("intercept","se_intercept","beta_age","se_beta_age","beta_age_squared","se_beta_age_squared")
Qe="QEp2"
}



data_combine <- NULL
for (chr in 1:22)
{
print(sprintf("loading chromosome %i",chr))
load(sprintf("Meta_regression_%s_chr%s_%s.Rdata",var2,chr,analysis))
data_combine <- rbind(data_combine,data)
print(sprintf("after combining, the number of SNPs is %i",dim(data_combine)[1]))
}

data_combine[,Qp] <- as.numeric(data_combine[,Qp])
print(sprintf("minimal p-value for %s = %e",var2,min(data_combine[,Qp],na.rm=T)))

data_combine[,c("chr","bp")] <- str_split_fixed(str_split_fixed(data_combine[,1],"# ",2)[,2],":",3)[,1:2]
data_combine[,"chr2"] <- as.numeric(data_combine$chr)
data_combine[,"bp2"] <- as.numeric(data_combine$bp)

compute_samplesize <- function(snp) {
	return(sum(N[which(!is.na(data_combine[snp, cbind(sprintf("effect_%s", inlist))]))]))
}
data_combine[,"TotalSampleSize"] <- unlist(lapply(1:dim(data_combine)[1],compute_samplesize))
data_combine <- data_combine[,c('markername','chr2','bp2','allele1','allele2',Qp,'TotalSampleSize',effect_names,Qe)]
names(data_combine) <- c('markername','chr','bp','allele1','allele2','P_age','TotalSampleSize',nice_effect_names,"P_resid_heterogeneity")

write.table(data_combine,sprintf("Meta_%s_%s.txt",var2,analysis),col.names=T,row.names=F,quote=F)

if (doplot == 1)
{
png(sprintf("Meta_regression_Manhattan_%s_%s.png",var2,analysis))
manhattan(data_combine[which(is.finite(data_combine[,'P_age'])),],p="P_age",chr="chr",bp="bp",snp="markername",suggestiveline=F,xlab=sprintf("based on %i SNPs",dim(data_combine[which(is.finite(data_combine[,'P_age'])),])[1]))
dev.off()

# create QQ plot 
chisq <- qchisq(1-data_combine[which(is.finite(data_combine[,'P_age'])),'P_age'],1)
GCl=median(chisq)/qchisq(0.5,1)
png(sprintf("Meta_regression_qq_%s_%s.png",var2,analysis))
qq(data_combine[,"P_age"],xlab=sprintf("genomic inflation factor = %f",GCl))
dev.off()
}

