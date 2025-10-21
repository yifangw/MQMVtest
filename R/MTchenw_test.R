###########################################################
# mean-based method: MTchenw
###########################################################

#' @title mean-based method: MTchenw, including a dominant effect of heterozygous (Aa) females
#' @description
#' A function to obtain the p value and the test statistic of MTchenw for testing the mean difference of the trait value across genotypes.
#'
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate Optional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15).
#' Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' When it is NULL, 1/ sqrt(2 SampleSize) will be used (Ionita-Laza et al. 2013). Only common variants are included in the analysis.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males.
#' @param kins For the mixed data, kins is a block matrix, interpreted as the genetic relatedness among these individuals,
#' including general pedigrees and unrelated individuals.
#' For general pedigrees, it is a kinship matrix, indicating the genetic relatedness among individuals within general pedigrees,
#' calculated using the dedicated method for kinship coefficients of X chromosome provided by the “kinship2” package in R software.
#'
#' @return The p value and test statistic of MTchenw.
#' @export
#'
#' @examples
MTchenw_test <- function(Genotype,Y,Sex,Covariate=NULL,
                        missing_cutoff=0.15,
                        MAF_Cutoff=NULL,
                        MGC_Cutoff=20,
                        kins=kins){

  if (!is.null(kins) && !class(kins)[1] %in% c("matrix", "list")) {
    if (is.null(attr(class(kins), "package")))
      stop("Error: \"kins\" must be a matrix or a list.")
    else if (attr(class(kins), "package") != "Matrix")
      stop("Error: if \"kins\" is a sparse matrix, it must be created using the Matrix package.")
  }
  if (missing(Genotype)){
    stop("The Genotype input is missing.")
  }
  if (missing(Y)){
    stop("The quantitative trait input is missing.")
  }
  if (missing(Sex)){
    stop("The Sex input is missing.")
  }
  if(!all(unique(Sex) %in% c(1,2))){
    stop('Sex must be a vector of 1(males) and 2(females)')
  }
  if (length(table(Sex))==1){
    stop("Only Males or Females detected")
  }

  Sex[Sex==1] <- 0
  Sex[Sex==2] <- 1

  if(!is.null(Covariate) & (is.matrix(Covariate)|
                            is.vector(Covariate)|
                            is.data.frame(Covariate))){
    Covariate <- as.data.frame(Covariate)
    if(length(unique(c(length(Y),nrow(Covariate),nrow(Genotype),length(Sex))))!=1){
      stop("Make sure the inputs have the same length!")
    }
    Phedata <- cbind(Y,Covariate)
    Null_Model_string <- paste('Y',paste(colnames(Covariate), collapse = "+"),sep = "~")

  }else if(is.null(Covariate)){
    Phedata <- as.data.frame(Y)
    if(length(unique(c(length(Y),nrow(Genotype),length(Sex))))!=1){
      stop("Make sure the inputs have the same length!")
    }
    Null_Model_string <- 'Y~1'
  }else{
    stop("Covariate must be a vector or a matrix!")
  }
  if(is.vector(Genotype)){
    Genotype <- as.matrix(Genotype)
  }
  if(is.null(MAF_Cutoff)){
    MAF_Cutoff <- 1/sqrt(2*length(Sex))
  }
  # MAF_w <- ifelse(Sex==0,0.5,1)
  MAF_w <- ifelse(Sex==0,1,0.5)
  drop_index <- apply(Genotype,2,function(snp){
    table_snp_sex <- table(factor(interaction(snp, Sex)))
    index1 <- length(table_snp_sex) < 5
    missing_ratio <- mean(is.na(snp))
    MAF_Genotype <- mean(snp*MAF_w,na.rm=TRUE)
    index2 <- missing_ratio>=missing_cutoff
    index3 <- MAF_Genotype<MAF_Cutoff
    index4 <- any(table_snp_sex<MGC_Cutoff)
    index <- index1 | index2 | index3 | index4
  })
  Genotype <- Genotype[,!drop_index,drop=F]
  if(ncol(Genotype)==0){
    stop("No snp meets the inclusion criteria!")
  }

  myfit_string <- paste(Null_Model_string,'+snp*Sex+snp_D',sep = "") #It includes dominance effect terms


  res <- t(apply(Genotype, 2, function(snp){

    # snp <- Genotype[,1]
    # snp_D <- ifelse(Sex==0 & snp==1,1,0) #dominance effect
    snp_D <- ifelse(Sex==1 & snp==1,1,0) #dominance effect

    cluster <- matrix(0,nrow(Genotype),1)
    cluster[Sex==1 & snp==0,1] <- 1
    cluster[Sex==1 & snp!=0 & snp!=2,1] <- 2
    cluster[Sex==1 & snp==2,1] <- 3
    cluster[Sex==0 & snp==0,1] <- 4
    cluster[Sex==0 & snp!=0,1] <- 5

    df4<- data.frame(num=c(1:nrow(Genotype)),Phedata,snp,snp_D,Sex,cluster=cluster)

    #kins=GRM
    lmm_chen <- glmmkin(fixed=myfit_string,data = df4,kins=kins,id="num",groups="cluster",
                        family=gaussian(link="identity"),na.action = na.exclude)
    coef_chen <- lmm_chen$coefficients[c("snp","snp_D","snp:Sex")]
    cov_chen <- lmm_chen$cov[c("snp","snp_D","snp:Sex"),c("snp","snp_D","snp:Sex")]
    wald_chen <- t(coef_chen)%*%solve(cov_chen)%*%coef_chen
    chen_p <- pchisq(q = wald_chen,df = 3,lower.tail = FALSE) #chi-square(3)

    res <- c(wald_chen,chen_p)
    names(res) <- c('MTchenw.test','MTchenw.p')
    return(res)
  }))
  SNPnames <- rownames(res)
  if(is.null(SNPnames)){
    SNPnames <- paste('snp',seq(ncol(Genotype)),sep = "")
  }
  Tstat <- data.frame(SNP=SNPnames,res[,1,drop=F],row.names = NULL)
  colnames(Tstat) <- c('SNP','MTchenw')
  pvalue <- data.frame(SNP=SNPnames,res[,2,drop=F],row.names = NULL)
  colnames(pvalue) <- c('SNP','MTchenw')
  results <- list(Tstat=Tstat,pvalue=pvalue)
  return(results)
}


