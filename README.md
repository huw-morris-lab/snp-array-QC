# snp-array-QC
Quality control steps for SNP array data, e.g. NeuroChip

Date: April 2020

Last updated: 22/04/2020

Authors: Manuela Tan

## General description and purpose

QC and data cleaning for SNP array data (e.g. NeuroChip)

This covers:
1. Extracting a particular set of individuals from a wider dataset
2. Variant filtering
3. Sample/individual filtering
4. Sex checking
5. HWE filtering (optional)
6. Sample relatedness - IBD
7. PCA for ancestry

Recommend looking at these websites if you are not familiar with plink.

http://zzz.bwh.harvard.edu/plink/

https://www.cog-genomics.org/plink/


When you see code with backslash this just means continue the same command but on a new line. 

Also I recommend always looking at the plink log/output from your commands, check the total genotyping rate, number of variants, and number of individuals to make sure everything looks sensible.


# 1. Extract set of individuals

You need to make a list of individuals that you want to keep from your wider dataset. E.g. if I just want the QSBB PD cases from the whole NeuroChip data.

You should do you QC steps (variant filtering etc.) just in the individuals that you are using for your analysis, not all the individuals you have in your NeuroChip data. This is because different groups may come from slightly different populations, e.g. if I am just doing a GWAS in PD cases, I don't want to include any controls (even if my study recruited controls) because the cases and controls may be sampled from slightly different populations. This will change the variants and individuals that you remove during your QC, as the distributions and means etc. will be different depending on who you incude.

Make a tab or separated .txt file with the FID and IID of the individuals you want to keep (e.g. mine is called QSBB_PD_extract.txt)

Extract these cases

```
plink --bfile NeuroChip \
--keep QSBB_PD_extract.txt \
--make-bed \
--out NeuroChip.QSBB_PD
```

# 2. QC - sample and SNP filtering

Variant filtering - remove variants with MAF < 1% and genotyping rate < 95%

```
plink --bfile NeuroChip.QSBB_PD \
--geno 0.05 \
--maf 0.01 \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01
```


Sample filtering - generate statistics for missingness and heterozygosity

```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01 \
--missing \
--het \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc
```


Plot sample heterozygosity and genotyping rates. There are some QC R scripts in the /kronos/hmorris/QC_Rscripts folder.

This command uses the sampleqc.R script to make a PDF scatterplot of your sample missingness vs. heterozygosity. Helpful to visualise and see if there are outliers. This is optional as you can also use the R code below to make these plots and format however you want.

```
R --vanilla --slave \
--args NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.imiss NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.het \
imiss_het_plot.pdf < /data/kronos/hmorris/QC_Rscripts/sampleqc.R
```
This makes a plot called imiss_het_plot.pdf (you can change the name).


Plot sample heterozygosity and genotyping rates in R. This code makes a text file with samples to remove with FID and IID. Remove samples who do not meet call rate (>98%) or heterozygosity (2SDs away from mean) cutoffs. You can edit cutoffs depending on your data.

```
library(tidyverse)

#Read in sample heterozygosity table
sample_het <- read.table("NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.het", header = TRUE)

#Read in sample call rate table
sample_callrate <- read.table("NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.imiss", header = TRUE)

#Calculate proportion of heterozygosity
sample_het <- sample_het %>% 
  mutate(het = (N.NM. - O.HOM.)/N.NM.)

#Calculate mean and SD of heterozygosity
summary <- sample_het %>% 
  summarise(mean_het = mean(het),
            sd_het = sd(het))

mean_het <- summary[1,1]
sd_het <- summary[1,2]

#Write list of samples who are > 2 SDs from mean of heterozygosity
sample_het <- sample_het %>% 
  mutate(remove_het = ifelse(het > 2*sd_het + mean_het, "remove",
                             ifelse(het < mean_het - 2*sd_het, "remove", "keep")))


#Merge with callrate table
sample_stats <- sample_het %>% 
  left_join(sample_callrate, by = c("FID", "IID"))

#Calculate genotyping rate (1 minus missing rate)
sample_stats <- sample_stats %>% 
  mutate(callrate = 1 - F_MISS)

#Plot scatterplot
ggplot(data = sample_stats, mapping = aes(x = het, y = callrate, color = remove_het)) +
  geom_point() +
  theme_bw()

#Write list of samples to remove - if heterozygosity outliers or if callrate is <98%
samples_to_remove <- sample_stats %>% 
  filter(remove_het == "remove" | callrate < 0.98) %>% 
  select(FID, IID)

#Export as text file with just FID and IID
write.table(samples_to_remove, "samples_to_remove.txt",
            quote=F, col.names = F, row.names = F)

```

Now you should have a .txt file with the FIDs and IIDs of individuals you want to exclude, who do not meet call rate (>98%) or heterozygosity (2SDs away from mean) cutoffs.

Exclude these individuals using plink.

```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01 \
--remove samples_to_remove.txt \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD
```

# 3. Sex checking 

First update clincal genders in fam file (if this has not been added in GenomeStudio).

Need to create a .txt file with FID, IID and clinical sex (1 or M = male, 2 or F = female, 0 = missing)

```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD \
--update-sex QSBB_clinicalGenders.txt \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex
```

If you look at the fam file now you should see the sex column has been updated.

Sex checking in plink
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex \
--check-sex 0.2 0.7 \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexcheck_results
```

Can open the .sexcheck file in text editor. Check the mismatches. If there are lot, espcially on a plate, then there may have been a problem with the genotyping (plate flip, sample mixup)


From your sexcheck results, write list of samples that pass sexcheck (FID and IID). If the clinical gender is missing (column 3 is 0) then these samples will still be included in this step. But you should definitely try to get the clinical gender for these samples
```
cat NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexcheck_results.sexcheck | awk '($3=="0" || $5=="OK") {print $1 "\t"$2}' > sex_samples_to_keep.txt
```

Remove sex discordant samples using plink
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex \
--keep sex_samples_to_keep.txt \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass
```

# 4. Hardy Weinberg Equilibrium filtering

You could do this filtering step in controls only, not including cases. 

I have applied this filter when working with PD cases only. It is up to you.

Also if you are working with rare variants (e.g. not GWAS), you may not want to apply HWE filters for your rare variants. You could filter for common variants first, then write the list of SNPs that do not meet HWE, and extract these from the full dataset of rare + common variants.


Generate stats file which shows HWE stats for each SNP (this doesn't filter any variants, just makes some new files with the stats.)
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass \
--hardy \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe
```

Filter out variants with HWE p value < 0.00001 (you decide what cutoff you want to use).
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass \
--hwe 0.00001 \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe
```

# 5. IBD

You don't need to do this step if you are doing family-based studies etc.

First create pruned list of variants (independent SNPs, removed variants that are in linkage) (without restriction on MAF)
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe \
--indep-pairwise 50 5 0.05 \
--out NeuroChip.pruned
```


Run IBD only on the pruned SNP list - called prune.in

The min 0.1 means that plink will only output pairs of samples that have PI-HAT > 0.1. You can adjust this if you want to look at samples that are more distantly related
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe \
--extract NeuroChip.pruned.prune.in \
--genome \
--min 0.1 \
--out NeuroChip.IBD
```


Look at your related samples (can open the IBD.genome file in text editor or Excel). PI-HAT of 1 indicates that the samples are the same individual or identical twins. PI-HAT of 0.5 indicates parent/child relationship.

You need to use some judgement here to decide which samples to remove. If you see one sample is related to lots of other people, this may indicate sample contamination. You can remove one individual from each pair

However if one related pair with PI-HAT close to 1 also has very similar sample IDs, I have removed both because this suggests there has been sample mixup.


Write a list of individuals to remove (FID and IID in .txt file).

Remove related individuals (pi-hat > 0.1)
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe \
--remove IBD_remove.txt \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe.IBD_0.1
```

# 6. PCA

This does PCA using pruned SNP list. Only merging with CEU individuals

HapMap SNP list and binary files are available in /data/kronos/NGS_Reference/HapMap_Reference/.

I made a copy of these files in my local directory

Extract HapMap SNPs from your dataset.
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe.IBD_0.1 \
--extract /data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt \
--make-bed \
--out NeuroChip.hapmap_SNPs
```

Extract CEU individuals from HapMap dataset.

Downloaded from https://ftp.ncbi.nlm.nih.gov/hapmap/samples_individuals/relationships_w_pops_121708.txt

I filtered in R to just get the IDs of the CEU individuals from this list (included in the PCA R code later on).

```
plink --bfile /data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps \
--keep /data/kronos/mtan/reference/hapmap/HapMap_CEU.txt \
--make-bed \
--out /data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only
```

Merge your data with HapMap CEU data and extract pruned SNPs
```
plink --bfile NeuroChip.hapmap_SNPs \
--bmerge /data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.bed \
/data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.bim \
/data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.fam \
--extract NeuroChip.pruned.prune.in \
--make-bed \
--out NeuroChip.hapmap_SNPs.CEU_only.merged-pruned
```

At first pass this will come up witth errors for variants with 3+ alleles present. This may be because the alleles are swapped so don't match in your data to the HapMap data. Need to flip missnps.
```
plink --bfile NeuroChip.hapmap_SNPs \
--flip NeuroChip.hapmap_SNPs.CEU_only.merged-pruned-merge.missnp \
--make-bed \
--out NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps
```

Remerge and extract pruned SNPs. Only CEU individuals from HapMap.
```
plink --bfile NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps \
--bmerge /data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.bed \
/data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.bim \
/data/kronos/mtan/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.CEU_only.fam \
--extract NeuroChip.pruned.prune.in \
--make-bed \
--out NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned
```

Check how many variants you have here and what the genotyping rate is. Cornelis recommends you need at least 10K independent SNPs (linkage pruned) to do your PCA, so if you have less you should check why.

Create genetic principal components from the merged (your data + HapMap CEU) and pruned dataset

```
gcta --bfile NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned \
--make-grm \
--autosome \
--out NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned.matrix

#Run PCA to generate 10 principal components
gcta --grm NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned.matrix --pca 10
```

Alternatively you can run the PCA in plink, rather than gcta. These come up with similar results
```
plink --bfile NeuroChip.hapmap_SNPs.CEU_only.flipped_missnps.merged-pruned \
--pca \
--out PCA.plink
```


In R, I run this code to make PCA plots, and writes list of PCA outliers to remove who are >6SD away from the mean of any of the first 10 PCs. Adjust cutoffs if necessary.

```
library(tidyverse)

#---Load principal components generated from GCTA---####
#This is the PCA from your data merged with HapMap
#Using linkage-pruned SNPs, MAF > 5%, excluding palindromic SNPs and flipping missnps

#Read in eigenvectors
geneticPCA_hapmap.eigenvec <- as_tibble(read.table("gcta.eigenvec", sep = ""))

#Change column names
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>% 
  dplyr::rename(FID = V1,
                IID = V2,
                PC1 = V3,
                PC2 = V4,
                PC3 = V5,
                PC4 = V6,
                PC5 = V7,
                PC6 = V8,
                PC7 = V9,
                PC8 = V10,
                PC9 = V11,
                PC10 = V12)


#---Read in data about HapMap samples (ethnicity)---####

#Read in population data for HapMap samples (adjust your directory to wherever you have kept this file).
HapMap_pops <- as_tibble(read.table("../../../../reference/hapmap/relationships_w_pops_121708.txt",
                                    header = TRUE))

#---Write list of CEU HapMap samples only---####
#To extract from HapMap dataset

#Write list of CEU HapMap samples only
HapMap_pops_CEU <- HapMap_pops %>% 
  filter(population == "CEU") %>% 
  select(FID, IID)

write.table(HapMap_pops_CEU, "../../../../reference/hapmap/HapMap_CEU.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#---Merge PROBAND PCA data with HapMap population data---####

geneticPCA_hapmap.eigenvec$FID <- as.factor(geneticPCA_hapmap.eigenvec$FID)

#Join with PCA data table
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>% 
  left_join(HapMap_pops)

#If the population is missing, these should be our samples (check this)
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>% 
  dplyr::mutate(group = ifelse(!is.na(population), population, "QSBB"))

#Make new variable for group (population)
geneticPCA_hapmap.eigenvec <- geneticPCA_hapmap.eigenvec %>% 
  mutate(group = ifelse(group == 2, "CEU",
                        ifelse(group == "QSBB", "QSBB", NA)))

#---Plot first two Principal Components with HapMap samples---####

ggplot(data = geneticPCA_hapmap.eigenvec, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_color_discrete(breaks=c("QSBB","CEU")) +
  theme_bw() +
  ggsave("PCA_HapMap_CEUsamples.png")

#---Classify outliers who are > 6 SDs away from mean for any of the first 10 PCs---####

#Individuals who were more than 6 standard deviations away from the mean of the any of the first 10 principal components were removed.
#Each PC mean is calculated, and any individual who is away from the mean on any of the PCs is removed

#Look at the mean for each PC
geneticPCA.eigenvec.PCmeans <- geneticPCA_hapmap.eigenvec %>% 
  dplyr::mutate(mean_PC1 = mean(PC1),
                mean_PC2 = mean(PC2),
                mean_PC3 = mean(PC3),
                mean_PC4 = mean(PC4),
                mean_PC5 = mean(PC5),
                mean_PC6 = mean(PC6),
                mean_PC7 = mean(PC7),
                mean_PC8 = mean(PC8),
                mean_PC9 = mean(PC9),
                mean_PC10 = mean(PC10))

## Remove individuals who are outliers on any PC

#Make results table, with one row for each individual
PC.outlierResults <- as.data.frame(matrix(ncol = 11, nrow = nrow(geneticPCA_hapmap.eigenvec)))
PC.outlierResults[, 1] <- geneticPCA_hapmap.eigenvec$IID #Put your IIDs in the first column


#For loop to calculate the SDs of each Principal Component (first 10 PCs only)
#This outputs into a results table
for (i in 3:12) {
  mean <- mean(geneticPCA.eigenvec.PCmeans[[i]])
  sd <- sd(geneticPCA.eigenvec.PCmeans[[i]])
  
  PC.outlierResults[, i-1] <- geneticPCA.eigenvec.PCmeans %>% 
    mutate(outlier = ifelse(geneticPCA.eigenvec.PCmeans[[i]] > mean + 6*sd, "outlier",
                            ifelse(geneticPCA.eigenvec.PCmeans[[i]] < mean - 6*sd, "outlier", "keep"))) %>% 
    dplyr::select(outlier)
}

#Rename column names
PC.outlierResults <- PC.outlierResults %>% 
  dplyr::rename(ID = V1,
                PC1_result = V2,
                PC2_result = V3,
                PC3_result = V4,
                PC4_result = V5,
                PC5_result = V6,
                PC6_result = V7,
                PC7_result = V8,
                PC8_result = V9,
                PC9_result = V10,
                PC10_result = V11)

#Now merge the outlier results with the main dataset
geneticPCA.eigenvec.PCmeans <- geneticPCA.eigenvec.PCmeans %>% 
  left_join(PC.outlierResults, by = c("IID" = "ID"))

#If any of the PC results are outliers, flag as outlier
geneticPCA.eigenvec.PCmeans <- geneticPCA.eigenvec.PCmeans %>% 
  mutate(PCA_outlier = ifelse(PC1_result == "outlier" |
                                PC2_result == "outlier" |
                                PC3_result == "outlier" |
                                PC4_result == "outlier" |
                                PC5_result == "outlier" |
                                PC6_result == "outlier" |
                                PC7_result == "outlier" |
                                PC8_result == "outlier" |
                                PC9_result == "outlier" |
                                PC10_result == "outlier", "outlier final", "keep final"))

#Plot first 2 PCs by outlier status (this includes both the HapMap samples and your samples).
ggplot(data = geneticPCA.eigenvec.PCmeans, mapping = aes(x = PC1, y = PC2, color = PCA_outlier)) +
  geom_point(alpha = 0.2) +
  theme_bw()

#Plot third and fourth PC by outlier status
ggplot(data = geneticPCA.eigenvec.PCmeans, mapping = aes(x = PC3, y = PC4, color = PCA_outlier)) +
  geom_point(alpha = 0.2) +
  theme_bw()

#Count how many outliers
geneticPCA.eigenvec.PCmeans %>% 
  group_by(PCA_outlier) %>% 
  filter(group == "QSBB") %>% 
  dplyr::summarise(count = n())

#---Write list of samples that are population outliers---####
#To remove in PLINK

PCA_outliers <- geneticPCA.eigenvec.PCmeans %>% 
  filter(PCA_outlier == "outlier final") %>% 
  filter(group == "QSBB") %>% 
  select(FID, IID)

#Write text file of FID and IID
write.table(PCA_outliers, "PCA_outliers.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
```

Remove PCA outliers using plink, using the list of individuals that you just made in R.
```
plink --bfile NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe.IBD_0.1 \
--remove PCA_outliers.txt \
--make-bed \
--out NeuroChip.QSBB_PD.geno_0.95.maf_0.01.sampleqc.sample_0.98.het_2SD.updatedsex.sexpass.hwe.IBD_0.1.PCA_keep
```


Good luck!!
