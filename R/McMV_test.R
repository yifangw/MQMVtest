###################################################################################
# mean-variance-based method:McMV
###################################################################################

#' @title mean-variance-based method:McMV,simultaneously test for differences in both the means and variances of a quantitative trait.
#'
#'@description
#'A function to obtain the p value and the test statistics of the MTchen_test (i.e., MTchen), the scale test (i.e.,MwM3VNA), the McMV_test(i.e., McMV) or all.
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate ptional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males.
#' @param kins For the mixed data, kins is a block matrix, interpreted as the genetic relatedness among these individuals,
#' including general pedigrees and unrelated individuals.
#' For general pedigrees, it is a kinship matrix, indicating the genetic relatedness among individuals within general pedigrees,
#' calculated using the dedicated method for kinship coefficients of X chromosome provided by the “kinship2” package in R software.
#' @param method Optional: A character string indicating which kind of association tests is to be conducted.
#' There are four options: "MTchen", "scale", "joint" (default) and "all".
#' method="MTchen": MTchen; method="scale": MwM3VNA; method="joint": McMV; method="all": all of the above association tests.
#'
#' @return The p values and the test statistics of association tests selected by the method option for each SNP.
#' @export
#'
#' @examples McMV_test(Genotype,Y,Sex,
#' Covariate=Data[,"Age"],
#' missing_cutoff=0.15,
#' MAF_Cutoff=NULL,
#' MGC_Cutoff=20,
#' kins=GRM,
#' method='joint')
McMV_test<- function(Genotype,Y,Sex,
                     Covariate=NULL,
                     missing_cutoff=0.15,
                     MAF_Cutoff=NULL,
                     MGC_Cutoff=20,
                     kins=NULL,
                     method='joint'){

  if(method=='MTchen') {
    MTchen <- MTchen_test(Genotype,Y,Sex,Covariate,
                          missing_cutoff,
                          MAF_Cutoff,MGC_Cutoff,
                          kins=kins)
    return(MTchen)
  } else if(method=='scale'){

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)

    return(scale)
  } else if(method=='joint'){
    MTchen <- MTchen_test(Genotype,Y,Sex,Covariate,
                          missing_cutoff,
                          MAF_Cutoff,
                          MGC_Cutoff,
                          kins)
    scale <- scale_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)
    McMV <- matrix(nrow = nrow(MTchen$pvalue),ncol=2)
    log_MTchen <- log(MTchen$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])

    McMV[,1] <- -2*(log_MTchen[,1]+log_scale[,1])
    McMV[,2] <- pchisq(McMV[,1],df = 4, lower.tail = FALSE)

    colnames(McMV) <- c('McMV.test','McMV.p')
    Tstat <- data.frame(SNP=MTchen$pvalue[,1],McMV[,1,drop=F],row.names = NULL)
    colnames(Tstat) <- c('SNP','McMV')
    pvalue <- data.frame(SNP=MTchen$pvalue[,1],McMV[,2,drop=F],row.names = NULL)
    colnames(pvalue) <- c('SNP','McMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }else if(method=='all'){
    MTchen <- MTchen_test(Genotype,Y,Sex,Covariate,
                          missing_cutoff,
                          MAF_Cutoff,
                          MGC_Cutoff,
                          kins)
    scale <- scale_test(Genotype,Y,Sex,
                        Covariate,
                        missing_cutoff,
                        MAF_Cutoff,
                        MGC_Cutoff,
                        kins)
    McMV <- matrix(nrow = nrow(MTchen$pvalue),ncol=2)
    log_MTchen <- log(MTchen$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])
    McMV[,1] <- -2*(log_MTchen[,1]+log_scale[,1])
    McMV[,2] <- pchisq(McMV[,1],df = 4, lower.tail = FALSE)
    colnames(McMV) <- c('McMV.test','McMV.p')
    Tstat <- data.frame(MTchen$Tstat,
                        scale$Tstat[,-1,drop=F],
                        McMV[,1,drop=F],
                        row.names = NULL)
    colnames(Tstat) <- c('SNP','MTchen','MwM3VNA3.3','McMV')
    pvalue <- data.frame(MTchen$pvalue,
                         scale$pvalue[,-1,drop=F],
                         McMV[,2,drop=F],
                         row.names = NULL)
    colnames(pvalue) <- c('SNP','MTchen','MwM3VNA3.3','McMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }
}
