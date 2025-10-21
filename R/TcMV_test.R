###################################################################################
#Unrelated individuals
# mean-variance-based method:TcMV
###################################################################################

#' @title mean-variance-based method for unrelated individuals:TcMV, simultaneously test for differences in both the means and variances of a quantitative trait.
#'@description
#'A function to obtain the p value and the test statistics of the Tchenw_test (i.e., Tchenw), the scale test (i.e.,wM3VNA) from the "QMVtest" package, the TcMV_test(i.e., TcMV) or all.
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
#' There are four options: "Tchenw", "scale", "joint" (default) and "all".
#' method="Tchenw": Tchenw; method="scale": wM3VNA; method="joint": TcMV; method="all": all of the above association tests.
#'
#' @return The p values and the test statistics of association tests selected by the method option for each SNP.
#' @export
#'
#' @examples
TcMV_test <- function(Genotype,Y,Sex,
                      Covariate=NULL,
                      missing_cutoff=0.15,
                      MAF_Cutoff=NULL,
                      MGC_Cutoff=20,
                      method='joint'){

  if(method=='Tchenw') {

    Tchenw <- Tchenw_test(Genotype,Y,Sex,Covariate)
    return(Tchen_ind)

  } else if(method=='scale'){

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    return(scale)
  } else if(method=='joint'){

    Tchenw <- Tchenw_test(Genotype,Y,Sex,Covariate)

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    TcMV_ind <- matrix(nrow = nrow(Tchenw$pvalue),ncol=2)
    log_Tchenw <- log(Tchenw$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])

    TcMV_ind[,1] <- -2*(log_Tchenw[,1]+log_scale[,1])
    TcMV_ind[,2] <- pchisq(TcMV_ind[,1],df = 4, lower.tail = FALSE)

    colnames(TcMV_ind) <- c('TcMV.test','TcMV.p')
    Tstat <- data.frame(SNP=Tchenw$pvalue[,1],TcMV_ind[,1,drop=F],row.names = NULL)
    colnames(Tstat) <- c('SNP','TcMV')
    pvalue <- data.frame(SNP=Tchenw$pvalue[,1],TcMV_ind[,2,drop=F],row.names = NULL)
    colnames(pvalue) <- c('SNP','TcMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }else if(method=='all'){
    Tchen_ind <- Tchenw_test(Genotype,Y,Sex,Covariate)

    scale <- scale_test(Genotype,Y,Sex,Covariate,
                        missing_cutoff,
                        MAF_Cutoff,MGC_Cutoff)

    TcMV_ind <- matrix(nrow = nrow(Tchenw$pvalue),ncol=2)
    log_Tchenw <- log(Tchenw$pvalue[,-1,drop=F])
    log_scale <- log(scale$pvalue[,-1,drop=F])
    TcMV_ind[,1] <- -2*(log_Tchenw[,1]+log_scale[,1])
    TcMV_ind[,2] <- pchisq(TcMV_ind[,1],df = 4, lower.tail = FALSE)
    colnames(TcMV_ind) <- c('TcMV.test','TcMV.p')
    Tstat <- data.frame(Tchenw$Tstat,
                        scale$Tstat[,-1,drop=F],
                        TcMV_ind[,1,drop=F],
                        row.names = NULL)
    colnames(Tstat) <- c('SNP','Tchenw','wM3VNA3.3','TcMV')
    pvalue <- data.frame(Tchenw$pvalue,
                         scale$pvalue[,-1,drop=F],
                         TcMV_ind[,2,drop=F],
                         row.names = NULL)
    colnames(pvalue) <- c('SNP','Tchenw','wM3VNA3.3','TcMV')
    results <- list(Tstat=Tstat,pvalue=pvalue)
    return(results)
  }
}
