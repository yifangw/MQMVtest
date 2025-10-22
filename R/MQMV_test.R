###################################################################################
#mean-variance-based methods: MQMVXcat and MQMVZmax
###################################################################################

#' @title The mean-variance-based methods: MQMVXcat and MQMVZmax, where MQMVXcat is a method that accounts for different XCI patterns and
#' simultaneously test for differences in both the means and variances of a quantitative trait,
#' while MQMVZmax is a mean-variance-based method that considers different DC patterns.
#'
#'
#'@description
#'A function to obtain the p values and the test statistics of the mean_test (i.e., MQXcat and MQZmax),
#'the variance_test (i.e.,MwM3VNA), the MQMV_test (i.e., MQMVXcat and MQMVZmax) or all.
#'MQMVXcat and MQMVZmax are designed to test for both the mean differences and the variance heterogeneity of the trait value across genotypes.
#'MQXcat and MQZmax are used for testing the mean differences of the trait value only.
#'MwM3VNA is for testing the variance heterogeneity only.
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate  Optional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' When it is NULL, 1/ sqrt(2 SampleSize) will be used (Ionita-Laza et al. 2013). Only common variants are included in the analysis.
#' @param MGC_Cutoff  Cutoff for the minimum genotype count in either females or males.
#' @param kins For the mixed data, kins is a block matrix, interpreted as the genetic relatedness among these individuals,
#' including general pedigrees and unrelated individuals.
#' For general pedigrees, it is a kinship matrix, indicating the genetic relatedness among individuals within general pedigrees,
#' calculated using the dedicated method for kinship coefficients of X chromosome provided by the “kinship2” package in R software.
#' @param method  Optional: A character string indicating which kind of association tests is to be conducted.
#' There are four options: "mean", "variance", "joint" (default) and "all".
#' method="mean": MQXcat and MQZmax; method="variance": MwM3VNA; method="joint": MQMVXcat and MQMVZmax;
#' method="all": all of the above association tests.
#'
#' @return The p values and the test statistics of association tests selected by the method option for each SNP.
#' @export
#'
#' @examples   MQMV_test(Genotype,Y,Sex,
#'                       Covariate=mixed[,"Age"],
#'                       missing_cutoff=0.15,
#'                       MAF_Cutoff=NULL,
#'                       MGC_Cutoff=20,
#'                       kins=GRM,
#'                       method='joint')
MQMV_test <- function(Genotype,Y,Sex,
                      Covariate=NULL,
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      kins=NULL,
                      method='joint'){

  if(method=='mean') {
    mean <- mean_test(Genotype,Y,Sex,
                              Covariate,
                              missing_cutoff,
                              MAF_Cutoff,
                              MGC_Cutoff,
                              kins)
    return(mean)
  } else if(method=='variance'){
    variance <- variance_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)
    return(variance)
  } else if(method=='joint'){
    mean <- mean_test(Genotype,Y,Sex,Covariate,
                              missing_cutoff,
                              MAF_Cutoff,
                              MGC_Cutoff,
                              kins)
    variance <- variance_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)
    MQMV <- matrix(nrow = nrow(mean$pvalue),ncol=4)
    log_mean <- log(mean$pvalue[,-1,drop=F])
    log_variance <- log(variance$pvalue[,-1,drop=F])
    for (i in 1:2) {
      MQMV[,i] <- -2*(log_mean[,i]+log_variance[,1]) #i=1,MQMVXcat;i=2,MQMVZmax
      MQMV[,i+2] <- pchisq(MQMV[,i],df = 4, lower.tail = FALSE)
    }
    colnames(MQMV) <- c('MQMVXcat.test','MQMVZmax.test','MQMVXcat.p','MQMVZmax.p')
    Tstat <- data.frame(SNP=mean$pvalue[,1],MQMV[,1:2,drop=F],row.names = NULL)
    colnames(Tstat) <- c('SNP','MQMVXcat','MQMVZmax')
    pvalue <- data.frame(SNP=mean$pvalue[,1],MQMV[,3:4,drop=F],row.names = NULL)
    colnames(pvalue) <- c('SNP','MQMVXcat','MQMVZmax')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }else if(method=='all'){
    mean <- mean_test(Genotype,Y,Sex,
                              Covariate,
                              missing_cutoff,
                              MAF_Cutoff,
                              MGC_Cutoff,
                              kins)
    variance <- variance_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)
    MQMV <- matrix(nrow = nrow(mean$pvalue),ncol=4)
    log_mean <- log(mean$pvalue[,-1,drop=F])
    log_variance <- log(variance$pvalue[,-1,drop=F])
    for (i in 1:2) {
      MQMV[,i] <- -2*(log_mean[,i]+log_variance[,1])
      MQMV[,i+2] <- pchisq(MQMV[,i],df = 4, lower.tail = FALSE)
    }
    colnames(MQMV) <- c('MQMVXcat.test','MQMVZmax.test','MQMVXcat.p','MQMVZmax.p')
    Tstat <- data.frame(mean$Tstat,
                        variance$Tstat[,-1,drop=F],
                        MQMV[,1:2,drop=F],
                        row.names = NULL)
    colnames(Tstat) <- c('SNP','MQXcat','MQZmax','MwM3VNA','MQMVXcat','MQMVZmax')
    pvalue <- data.frame(mean$pvalue,
                         variance$pvalue[,-1,drop=F],
                         MQMV[,3:4,drop=F],
                         row.names = NULL)
    colnames(pvalue) <- c('SNP','MQXcat','MQZmax','MwM3VNA','MQMVXcat','MQMVZmax')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }
}
