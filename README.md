# MQMVtest
An R package of X chromosome-wide association studies for quantitative trait loci based on the mixture of general pedigrees and additional unrelated individuals. In this package, the nine methods (MQXcat, MQZmax, MTplink, MTchen, MwM3VNA, MQMVXcat, MQMVZmax, MpMV and McMV) can not only handle the mixed data but also be applied to only general pedigrees, where the latter case simply requires reducing the block matrix to a kinship matrix. Specifically, for the mixed data, mean_test() performs MQXcat and MQZmax; MTplinkw_test() and MTchenw_test() perform MTplinkw and MTchenw; variance_test() performs MwM3VNA; MQMV_test() performs MQMVXcat and MQMVZmax; MpMV_test() and McMV_test() perfrom MpMV and McMV. For general pedigrees, we also can use these functions to perfrom corresponding methods, these methods respectively denoted as PQXcat, PQZmax, PTplinkw, PTchenw, PwM3VNA, PQMVXcat, PQMVZmax, PpMV and PcMV.

# Installation
It is easy to install the development version of MQMVtest package using the 'devtools' package.
```r
# Install devtools if you haven't already
if (!require("devtools")) {
    install.packages("devtools")
}

# Load devtools
library(devtools)

# Install MQMVtest
install_github("yifangw/MQMVtest")
```

# Dependencies
```r
#install package dependencies

install.packages(c('GMMAT','lmtest','quantreg', 'gJLS2','mvtnorm', 'expm', 'bigsnpr', 'xgboost'))

```

# Usage
```r
# Required R packages
#============================================================#
#--------------------------------------------------------------#
# Usage Examples
#--------------------------------------------------------------#
#============================================================#



#============================================================#
# Loading packages
#============================================================#
library(GMMAT)
library(MQMVtest)
library(data.table)
library(kinship2)
library(RNOmni)

#============================================================#
# Loading required data
#============================================================#


#============================================================#
# step 1 for mixed data
#============================================================#

#==============================
# Calculate the kinship coefficient
#==============================
ped.ID <- pedigree(id=pedigrees$IID, dadid=pedigrees$PID, momid=pedigrees$MID,
                   sex=pedigrees$SEX, famid=pedigrees$FID, missid=0)
kin <- kinship(ped.ID,"x")
dim(kin)

kin_order = colnames(kin)
match_idx <- match(kin_order, pedigrees$IID)
class(kin_order)
ped_final_sort = pedigrees[match_idx, ]

## Convert to a genetic correlation matrix
kin[ped_final_sort$SEX==1,ped_final_sort$SEX==2]=sqrt(2)*kin[ped_final_sort$SEX==1,ped_final_sort$SEX==2]
kin[ped_final_sort$SEX==2,ped_final_sort$SEX==1]=sqrt(2)*kin[ped_final_sort$SEX==2,ped_final_sort$SEX==1]
kin[ped_final_sort$SEX==2,ped_final_sort$SEX==2]=2*kin[ped_final_sort$SEX==2,ped_final_sort$SEX==2]
kin<-as.matrix(kin)
kin1<-kin[row.names(kin) %in% mixed$IID,row.names(kin) %in% mixed$IID]

GRM <- matrix(0,nrow=nrow(mixed),ncol=nrow(mixed),dimnames=list(mixed$IID,mixed$IID))
diag(GRM) <- rep(1,nrow(mixed))
common_ids <- intersect(rownames(kin), mixed$IID)
GRM[common_ids, common_ids] <- kin[common_ids, common_ids]
GRM <- as.matrix(GRM)


GRM <- as(GRM, "sparseMatrix")
rownames(GRM) <- c(1:nrow(GRM))
colnames(GRM) <- c(1:nrow(GRM))

#-----------------------------------------------------------------------------#
##############  Run association studies for the mixed data  #############
#-----------------------------------------------------------------------------#

mixed$HDL = RankNorm(mixed$HDL)
colnames(mixed)[1:30]

mixed_results_1 = MQMV_test(Genotype = mixed[,19],
                            Y = mixed[,'HDL'],
                            Sex = mixed[,'SEX'],
                            Covariate = mixed[,c("Age_at_recruitment", "PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")],
                            missing_cutoff = 0.15,
                            MGC_Cutoff = 2,
                            kins = GRM,
                            method = 'all')
                            
# mixed_results_1
# $Tstat
#    SNP   MQXcat   MQZmax   MwM3VNA MQMVXcat MQMVZmax
# 1 snp1 6.678427 1.397362 0.4122579  2.94761 3.994008
# 
# $pvalue
#    SNP   MQXcat    MQZmax   MwM3VNA  MQMVXcat  MQMVZmax
# 1 snp1 0.307779 0.1823965 0.7442103 0.5666307 0.4068174

mixed_results_2 = MpMV_test(Genotype = mixed[,19],
                            Y = mixed[,'HDL'],
                            Sex = mixed[,'SEX'],
                            Covariate = mixed[,c("Age_at_recruitment", "PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")],
                            missing_cutoff = 0.15,
                            MGC_Cutoff = 2,
                            kins = GRM,
                            method = 'all')

mixed_results_3 = McMV_test(Genotype = mixed[,19],
                            Y = mixed[,'HDL'],
                            Sex = mixed[,'SEX'],
                            Covariate = mixed[,c("Age_at_recruitment", "PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")],
                            missing_cutoff = 0.15,
                            MGC_Cutoff = 2,
                            kins = GRM,
                            method = 'all')


#============================================================#
# step 2  for general pedigrees
#============================================================#

#==============================
# Calculate the kinship coefficient
#==============================
ped.ID <- pedigree(id=pedigrees$IID, dadid=pedigrees$PID, momid=pedigrees$MID,
                   sex=pedigrees$SEX, famid=pedigrees$FID, missid=0)
kin <- kinship(ped.ID,"x")
dim(kin)

#=================================
#Ensure that the order of the samples matches that of the kin when converting to the genetic correlation matrix
#=================================

kin_order = colnames(kin)
match_idx <- match(kin_order, pedigrees$IID)
ped_final_sort = pedigrees[match_idx, ]
ped_final_sort$IID[1:10]

## genetic correlation matrix
kin[ped_final_sort$SEX==1,ped_final_sort$SEX==2]=sqrt(2)*kin[ped_final_sort$SEX==1,ped_final_sort$SEX==2]
kin[ped_final_sort$SEX==2,ped_final_sort$SEX==1]=sqrt(2)*kin[ped_final_sort$SEX==2,ped_final_sort$SEX==1]
kin[ped_final_sort$SEX==2,ped_final_sort$SEX==2]=2*kin[ped_final_sort$SEX==2,ped_final_sort$SEX==2]
kin<-as.matrix(kin)


# for general peigrees
GRM = as(kin, "sparseMatrix")
rownames(GRM)=c(1:nrow(GRM))
colnames(GRM)=c(1:nrow(GRM))


#-----------------------------------------------------------------------------#
##############  Run association studies for general pedigrees  #############
#-----------------------------------------------------------------------------#
library(RNOmni)
pedigrees$HDL = RankNorm(pedigrees$HDL)
colnames(pedigrees)[1:30]

ped_results_1 = MQMV_test(Genotype = pedigrees[,19],
                          Y = pedigrees[,'HDL'],
                          Sex = pedigrees[,'SEX'],
                          Covariate = pedigrees[,"Age_at_recruitment"],
                          missing_cutoff = 0.15,
                          MGC_Cutoff = 2,
                          kins = GRM,
                          method = 'all')
ped_results_1
$Tstat
   SNP   MQXcat   MQZmax   MwM3VNA MQMVXcat MQMVZmax
1 snp1 8.874293 1.786261 0.9477267 5.852146  6.68318

$pvalue
   SNP    MQXcat     MQZmax   MwM3VNA  MQMVXcat  MQMVZmax
1 snp1 0.1286391 0.08490163 0.4167252 0.2104655 0.1536083                          

ped_results_2 = MpMV_test(Genotype = pedigrees[,19],
                          Y = pedigrees[,'HDL'],
                          Sex = pedigrees[,'SEX'],
                          Covariate = pedigrees[,"Age_at_recruitment"],
                          missing_cutoff = 0.15,
                          MGC_Cutoff = 2,
                          kins = GRM,
                          method = 'all')
# ped_results_2
# $Tstat
#    SNP MTplinkw   MwM3VNA    MpMV
# 1 snp1 2.138653 0.9477267 3.88931
# 
# $pvalue
#    SNP  MTplinkw   MwM3VNA      MpMV
# 1 snp1 0.3432396 0.4167252 0.4211934

ped_results_3 = McMV_test(Genotype = pedigrees[,19],
                          Y = pedigrees[,'HDL'],
                          Sex = pedigrees[,'SEX'],
                          Covariate = pedigrees[,"Age_at_recruitment"],
                          missing_cutoff = 0.15,
                          MGC_Cutoff = 2,
                          kins = GRM,
                          method = 'all')
# ped_results_3
# $Tstat
#    SNP  MTchenw   MwM3VNA     McMV
# 1 snp1 4.406631 0.9477267 4.771912
# 
# $pvalue
#    SNP   MTchenw   MwM3VNA      McMV
# 1 snp1 0.2207713 0.4167252 0.3115113
#============================================================#
# step 3  for unrelated individuals
#============================================================#
# library(MQMVtest)
library(QMVtest)
library(RNOmni)
library(lmtest)


unrelated = inddata_filtered 

unre_results_1 = TpMV_test(Genotype = unrelated[,19],
                           Y = unrelated[,'HDL'],
                           Sex = unrelated[,'SEX'],
                           Covariate = unrelated[,c("Age_at_recruitment", "PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")],
                           missing_cutoff = 0.15,
                           MGC_Cutoff = 2,
                           method = 'all')

# unre_results_1
# $Tstat
#    SNP Tplinkw    wM3VNA     TpMV
# 1 snp1      NA 0.4310167 2.506031
# 
# $pvalue
#    SNP   Tplinkw    wM3VNA      TpMV
# 1 snp1 0.3908545 0.7308142 0.6435561

unre_results_2 = TcMV_test(Genotype = unrelated[,19],
                           Y = unrelated[,'HDL'],
                           Sex = unrelated[,'SEX'],
                           Covariate = unrelated[,c("Age_at_recruitment", "PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")],
                           missing_cutoff = 0.15,
                           MGC_Cutoff = 2,
                           method = 'all')
                    
# unre_results_2
# $Tstat
#    SNP Tchenw    wM3VNA     TcMV
# 1 snp1     NA 0.4310167 1.701122
# 
# $pvalue
#    SNP    Tchenw    wM3VNA      TcMV
# 1 snp1 0.5845196 0.7308142 0.7905138


```
# MQMVtest
