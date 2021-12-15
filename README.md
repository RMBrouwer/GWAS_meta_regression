# GWAS_meta_regression
This repository contains the script that was used to perform a genome-wide associaton meta-regression as was used in https://www.biorxiv.org/content/10.1101/2020.04.24.031138v2. Please make sure that there are enough cohorts in the study to make a regression feasible. 

## Prerequisites
- Cohort specific summary statistics for the phenotype that are ready for meta-analysis (QCed and filtered). 
- A summary file (see below for an example) that contains the necessary information per cohort.

## Step 1: align SNPs and store per cohort, per chromosome. 
Since one probably also wants to run a standard meta-analysis on the data, we will make use of METAL<sup>1</sup> - a program for meta-analysis that performs alignment - using the VERBOSE option to write all the aligned SNPs to a temporary file. See https://genome.sph.umich.edu/wiki/METAL_Documentation for documentation. Separating the chromosomes is done in this step to limit the number of SNPs that need to be kept in memory in subsequent steps. Input for this script are tempfile that contains the verbose output of METAL, and a textfile containing a list of the cohorts that went into that analysis. If you rerun the analysis (e.g. adding more cohorts) you should also rerun this alignment and separation of the chromosomes.


Call:
```
./separate_chr.do TRAIT Cohortlist.txt METAL_verbose_tempfile
```

After running this step, directory structure should be like this:

cohort_input_meta_regression/aligned_TRAIT <br>
cohort_input_meta_regression/aligned_TRAIT/Cohort1 <br>
cohort_input_meta_regression/aligned_TRAIT/Cohort2 <br>
...

within each of these directories, there should be files per chromosome: Cohort1_meta_regression_TRAIT_chr1.gz to Cohort1_meta_regression_TRAIT_chr22.gz.

Files should look like this: 

```
# 1:100000012	t	  g	    -72.943	148.537	0.6234	0.248	Cohort1_TRAIT_filtered_input.txt
# 1:100000827	t	  c	    -132.176	140.684	0.3475	0.295	Cohort1_TRAIT_filtered_input.txt
# 1:100001201	t	  g	    -193.362	225.048	0.3902	0.093	Cohort1_TRAIT_filtered_input.txt
# 1:100001671	ct        c	    -121.134	155.883	0.4371	0.263	Cohort1_TRAIT_filtered_input.txt
# 1:100002155	g	  gttagt    -192.873	225.047	0.3914	0.093	Cohort1_TRAIT_filtered_input.txt
# 1:100002713	t	  c	    -64.611	213.359	0.762	0.104	Cohort1_TRAIT_filtered_input.txt
...
```

## Step 2: preprare the summary file
The summary file named "summary_measures_TRAIT.txt" contains the per-cohort information that is used in the meta-regression, age in this example. Make sure there are no missing values. All - and only - cohorts that are listed in this file are included in the meta-regression. For these cohorts, it is assumed that the aligned SNP data as described above is present. The (space separated) file should contain columns "cohort", "N" (indicating sample size) and "mean_age". Other columns are allowed but not used (unless the Meta_regression_GWAS script is adapted). 

Example file: 

```
cohort N mean_age
Cohort1 840 15.9
Cohort2 130 46.3
Cohort3 1095 87.0
...
```
## Step 3: run the meta-regression
The meta-regression is run per chromosome, using the metafor<sup>2</sup> package in R<sup>3</sup>. Additional R packages needed are data.table, iterators and foreach. There are currently two types of meta-regression that can be run: 
1) Effect_SNPC ~ b0 + b1\*age + ε under the null hypothesis that b1=0 (1 degree of freedom), and
2) Effect_SNPC ~ b0 + b1\*age + b2\*(age)^2 + ε under the null hypothesis that (b1=b2=0, 2 degrees of freedom). 

Call:
```
R --no-save < Meta_regression_GWAS.R --args TRAIT CHR TYPE
```
where TYPE is either "2_c1" or "3_c2". 
These numbers refer to the number of parameters in the model and the number of covariates that are hypothesised zero - "2_c1" refers to 1) above, testing for linear age effect versus no age effect, "3_c2" refers to 2) above, testing for quadratic age effects versus no age effect. 

Notes: 
1) Set your working directory in line 6
2) this script automatically filters out all the SNPs that are present in less than 50% of cohorts and in less than 50% of total number of subjects (lines 39-52). Adjust if necessary. 
3) Depending on the number of cohorts and number of SNPs included, this script may take a couple of hours to complete. Ideally one would run each chromosome separately on a cluster. 

Meta-regression output is stored in an .Rdata file per chromosome. 

## Step 4: combine the information per chromosome
Combine the information per chromosome for followup analyses. R packages needed are stringr and qqman. 
Call:
```
R --no-save < Combine_and_plot_meta_regression --args TRAIT TYPE

```

Notes: 
1) Set your working directory in line 11
2) Set whether you want Manhattan and qq plots in line 5

Output of this script is similar to a GWAS summary stats file; but instead of effect sizes there are regression estimates and standard errors for the meta-regression (intercept and linear/quadratic components), a p-value for the age-effect (see step 3) and a p-value for residual heterogeneity. PLEASE NOTE: one cannot directly use these files in followup analysis tools that expect a (direction) of effect as these depend on age. 

Example output:
```
markername chr bp allele1 allele2 P_age TotalSampleSize intercept se_intercept beta_age se_beta_age P_resid_heterogeneity
# 21:14601415 21 14601415 a g 0.90778 10626 -83.5286 100.9246 1.1524 1.9022 0.1651
# 21:14604190 21 14604190 t c 0.69777 10626 68.0641 110.7797 -0.8302 2.1392 0.8478
# 21:14608578 21 14608578 a t 0.8965 10626 66.2267 112.3745 -0.8694 2.1643 0.6244
```

### References
<sup>1</sup> Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).
<sup>2</sup> Viechtbauer, W. Journal of statistical software. J. Stat. Softw. 36, 1–48 (2010).
<sup>3</sup> The R Core Team. R: A language and environment for statistical computing. (2018).






