###################################################################################
#mean-variance-based method:MpMV
###################################################################################


#' @title mean-variance-based method:MpMV, simultaneously test for differences in both the means and variances of a quantitative trait,
#'
#'@description
#'A function to obtain the p value and the test statistics of the MTplinkw_test (i.e., MTplinkw), the variance test (i.e.,MwM3VNA), the MpMV_test(i.e., MpMV) or all.
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate Optional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males.
#' @param kins  For the mixed data, kins is a block matrix, interpreted as the genetic relatedness among these individuals,
#' including general pedigrees and unrelated individuals.
#' For general pedigrees, it is a kinship matrix, indicating the genetic relatedness among individuals within general pedigrees,
#' calculated using the dedicated method for kinship coefficients of X chromosome provided by the “kinship2” package in R software.
#' @param method Optional: A character string indicating which kind of association tests is to be conducted.
#' There are four options: "MTplinkw", "variance", "joint" (default) and "all".
#' method="MTplinkw": MTplinkw; method="variance": MwM3VNA; method="joint": MpMV; method="all": all of the above association tests.
#'
#' @return The p values and the test statistics of association tests selected by the method option for each SNP.
#' @export
#'
#' @examples MpMV_test(Genotype,Y,Sex,
#'                     Covariate=mixed[,"Age"],
#'                     missing_cutoff=0.15,
#'                     MAF_Cutoff=NULL,
#'                     MGC_Cutoff=20,
#'                     kins=GRM,
#'                     method='joint')
MpMV_test <- function(Genotype,Y,Sex,
                      Covariate=NULL,
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      kins=NULL,
                      method='joint'){

  if(method=='MTplinkw') {
    MTplinkw <- MTplinkw_test(Genotype,Y,Sex,Covariate,
                            missing_cutoff,
                            MAF_Cutoff,
                            MGC_Cutoff,
                            kins=kins)
    return(MTplinkw)
  } else if(method=='variance'){

    variance <- variance_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)

    return(variance)
  } else if(method=='joint'){
    MTplinkw <- MTplinkw_test(Genotype,Y,Sex,
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
    MpMV <- matrix(nrow = nrow(MTplinkw$pvalue),ncol=2)
    log_MTplinkw <- log(MTplinkw$pvalue[,-1,drop=F])
    log_variance <- log(variance$pvalue[,-1,drop=F])

    MpMV[,1] <- -2*(log_MTplinkw[,1]+log_variance[,1])
    MpMV[,2] <- pchisq(MpMV[,1],df = 4, lower.tail = FALSE)

    colnames(MpMV) <- c('MpMV.test','MpMV.p')
    Tstat <- data.frame(SNP=MTplinkw$pvalue[,1],MpMV[,1,drop=F],row.names = NULL)
    colnames(Tstat) <- c('SNP','MpMV')
    pvalue <- data.frame(SNP=MTplinkw$pvalue[,1],MpMV[,2,drop=F],row.names = NULL)
    colnames(pvalue) <- c('SNP','MpMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }else if(method=='all'){
    MTplinkw <- MTplinkw_test(Genotype,Y,Sex,
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
    MpMV <- matrix(nrow = nrow(MTplinkw$pvalue),ncol=2)
    log_MTplinkw <- log(MTplinkw$pvalue[,-1,drop=F])
    log_variance <- log(variance$pvalue[,-1,drop=F])
    MpMV[,1] <- -2*(log_MTplinkw[,1]+log_variance[,1])
    MpMV[,2] <- pchisq(MpMV[,1],df = 4, lower.tail = FALSE)
    colnames(MpMV) <- c('MpMV.test','MpMV.p')
    Tstat <- data.frame(MTplinkw$Tstat,
                        variance$Tstat[,-1,drop=F],
                        MpMV[,1,drop=F],
                        row.names = NULL)
    colnames(Tstat) <- c('SNP','MTplinkw','MwM3VNA','MpMV')
    pvalue <- data.frame(MTplinkw$pvalue,
                         variance$pvalue[,-1,drop=F],
                         MpMV[,2,drop=F],
                         row.names = NULL)
    colnames(pvalue) <- c('SNP','MTplinkw','MwM3VNA','MpMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }
}
