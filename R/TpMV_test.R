###################################################################################
#Unrelated individuals
# mean-variance-based method:TpMV
###################################################################################

#' @title mean-variance-based method for unrelated individuals:TpMV, simultaneously test for differences in both the means and variances of a quantitative trait.
#'
#'@description
#'A function to obtain the p value and the test statistics of the Tplinkw_test (i.e., Tplinkw), the scale test (i.e.,wM3VNA) from the "QMVtest" package, the TpMV_test(i.e., TpMV) or all.
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate ptional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males.
#' @param method Optional: A character string indicating which kind of association tests is to be conducted.
#' There are four options: "Tplinkw", "scale", "joint" (default) and "all".
#' method="Tplinkw": Tplinkw; method="scale": wM3VNA; method="joint": TpMV; method="all": all of the above association tests.
#'
#'
#' @return The p values and the test statistics of association tests selected by the method option for each SNP.
#' @export
#'
#' @examples
TpMV_test <- function(Genotype,Y,Sex,
                      Covariate=NULL,
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      method='joint'){

  if(method=='Tplinkw') {
    Tplinkw <- Tplinkw_test(Genotype,Y,Sex,Covariate)
    return(Tplinkw)
  } else if(method=='scale'){

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    return(scale)
  } else if(method=='joint'){
    Tplinkw <- Tplinkw_test(Genotype,Y,Sex,Covariate)

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    TpMV_ind <- matrix(nrow = nrow(Tplinkw$pvalue),ncol=2)
    log_Tplinkw <- log(Tplinkw$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])

    TpMV_ind[,1] <- -2*(log_Tplinkw[,1]+log_scale[,1])
    TpMV_ind[,2] <- pchisq(TpMV_ind[,1],df = 4, lower.tail = FALSE)

    colnames(TpMV_ind) <- c('TpMV.test','TpMV.p')
    Tstat <- data.frame(SNP=Tplinkw$pvalue[,1],TpMV_ind[,1,drop=F],row.names = NULL)
    colnames(Tstat) <- c('SNP','TpMV')
    pvalue <- data.frame(SNP=Tplinkw$pvalue[,1],TpMV_ind[,2,drop=F],row.names = NULL)
    colnames(pvalue) <- c('SNP','TpMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }else if(method=='all'){
    Tplinkw <- Tplinkw_test(Genotype,Y,Sex,Covariate)

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    TpMV_ind <- matrix(nrow = nrow(Tplinkw$pvalue),ncol=2)
    log_Tplinkw <- log(Tplinkw$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])
    TpMV_ind[,1] <- -2*(log_Tplinkw[,1]+log_scale[,1])
    TpMV_ind[,2] <- pchisq(TpMV_ind[,1],df = 4, lower.tail = FALSE)
    colnames(TpMV_ind) <- c('TpMV.test','TpMV.p')
    Tstat <- data.frame(Tplinkw$Tstat,
                        scale$Tstat[,-1,drop=F],
                        TpMV_ind[,1,drop=F],
                        row.names = NULL)
    colnames(Tstat) <- c('SNP','Tplinkw','wM3VNA3.3','TpMV')
    pvalue <- data.frame(Tplinkw$pvalue,
                         scale$pvalue[,-1,drop=F],
                         TpMV_ind[,2,drop=F],
                         row.names = NULL)
    colnames(pvalue) <- c('SNP','Tplinkw','wM3VNA3.3','TpMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }
}
