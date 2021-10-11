#!/bin/bash

[ $# -eq 0 ] && { echo "Usage: $0 phenotype include_cohorts.txt METAL_tempfile"; exit 1; }

k=$1 # phenotype
inlist=$2 # text file containing names of all the cohorts - matching the names of the inputfiles that went into METAL
tempFILE=$3 # the verbose output of METAL

mkdir cohort_input_meta_regression
mkdir cohort_input_meta_regression/aligned_${k}

for cohort in `cat $inlist` 
do
# creating tempfiles separating out the cohorts 
cat $tempFILE | grep $cohort > cohort_input_meta_regression/aligned_${k}/${cohort}_${k}

cd cohort_input_meta_regression/aligned_${k}
mkdir ${cohort}

for j in `seq 1 22` # looping over the chromosomes
do 
cat ${cohort}_${k} | grep "# $j:" > ${cohort}/${cohort}_meta_regression_${k}_chr${j}
done

mv ${cohort}_${k} ${cohort}
gzip -f ${cohort}/*

cd -
done

