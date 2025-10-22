
#=========================================================#
# The mean-based methods:MQXcat and MQZmax
#=========================================================#
#' @title The mean-based methods:MQXcat and MQZmax
#'
#'@description A function to obtain the p values and the test statistics of MQXcat and MQZmax for testing the mean difference
#'of the trait value across genotypes.
#'
#' @param Genotype A numeric genotype matrix with each row as a different individual and each column as a separate SNP.
#' Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele.
#' @param Y A numeric vector of a quantitative trait, such as height.
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate Optional: a vector or a matrix of covariates, such as age.
#' @param missing_cutoff Cutoff of the missing rates of SNPs (default=0.15). Any SNPs with missing rates higher than the cutoff will be excluded from the analysis.
#' @param MAF_Cutoff MAF cutoff for common vs. rare variants (default=NULL). It should be a numeric value between 0 and 0.5, or NULL.
#' When it is NULL, 1/ sqrt(2 SampleSize) will be used (Ionita-Laza et al. 2013). Only common variants are included in the analysis.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males.
#' @param kins For the mixed data, kins is a block matrix, interpreted as the genetic relatedness among these individuals,
#' including general pedigrees and unrelated individuals.
#' For general pedigrees, it is a kinship matrix, indicating the genetic relatedness among individuals within general pedigrees,
#' calculated using the dedicated method for kinship coefficients of X chromosome provided by the “kinship2” package in R software.
#'
#' @return The p values and test statistics of MQXcat and MQZmax
#' @export mean_test
#'
#' @examples mean_test(Genotype,Y,Sex,
#'                     Covariate=mixed[,"Age"],
#'                     missing_cutoff=0.15,
#'                     MAF_Cutoff=NULL,
#'                     MGC_Cutoff=20)
mean_test <- function(Genotype,Y,Sex,
                          Covariate=NULL,
                          missing_cutoff=0.15,
                          MAF_Cutoff=NULL,
                          MGC_Cutoff=20,
                          kins=NULL){
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

  Sex[Sex==1] <- 0 #Change the coding for males from 1 to 0
  Sex[Sex==2] <- 1 #Change the coding for females from 2 to 1

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
    Null_Model_string <- 'Y~1' #Y~1 means it only contains constant terms
  }else{
    stop("Covariate must be a vector or a matrix!")
  }
  if(is.vector(Genotype)){
    Genotype <- as.matrix(Genotype)
  }
  if(is.null(MAF_Cutoff)){
    MAF_Cutoff <- 1/sqrt(2*length(Sex))
  }

  #-------Exclusion criteria-------------------------#

  MAF_w <- ifelse(Sex==0,1,0.5) #Females weight-0.5，Males weight-1
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
  #-----------------------------------------------------------------------------#

  fnMatSqrtInverse = function(mA) {
    solve(expm::sqrtm(mA))
  }


  Phedataf <- Phedata[Sex==1,,drop=F]
  Phedatam <- Phedata[Sex==0,,drop=F]

  nf <- nrow(Phedataf)
  nm <- nrow(Phedatam)
  myfit_string <- paste(Null_Model_string,'+Sex+snp_AAf+snp_Af+snp_Am',sep = "")

  res <- t(apply(Genotype, 2, function(snp)
  {
    # snp <- Genotype[,1]
    # snp_AAf <- ifelse(Sex==0 & snp==2,1,0)
    # snp_Af <- ifelse(Sex==0 & snp!=0,1,0) #Genotype is AA or Aa for females
    # snp_Am <- ifelse(Sex==1 & snp!=0,1,0)
    #
    # cluster <- matrix(0,nrow(Genotype),1)
    # cluster[Sex==0 & snp==0,1] <- 1
    # cluster[Sex==0 & snp!=0 & snp!=2,1] <- 2
    # cluster[Sex==0 & snp==2,1] <- 3
    # cluster[Sex==1 & snp==0,1] <- 4
    # cluster[Sex==1 & snp!=0,1] <- 5


    snp_AAf <- ifelse(Sex==1 & snp==2,1,0)
    snp_Af <- ifelse(Sex==1 & snp!=0,1,0) #Genotype is AA or Aa for females
    snp_Am <- ifelse(Sex==0 & snp!=0,1,0)

    cluster <- matrix(0,nrow(Genotype),1)
    cluster[Sex==1 & snp==0,1] <- 1
    cluster[Sex==1 & snp!=0 & snp!=2,1] <- 2
    cluster[Sex==1 & snp==2,1] <- 3
    cluster[Sex==0 & snp==0,1] <- 4
    cluster[Sex==0 & snp!=0,1] <- 5

    # timestart <- Sys.time()
    pheno_data <- data.frame(num=c(1:nrow(Genotype)),Phedata,Genotype,snp_AAf,snp_Af,snp_Am,Sex,cluster=cluster)
    LMM <- glmmkin(fixed=myfit_string, data=pheno_data, kins=kins,id="num",groups="cluster",
                   family=gaussian(link="identity"))
    # timeend <- Sys.time()
    # duration <- timeend - timestart
    #
    # timestart <- Sys.time()
    # pheno_data <- data.frame(num=c(1:nrow(Genotype)),Y,Genotype,snp_AAf,snp_Af,snp_Am,Sex,cluster=cluster)
    # LMM2 <- glmmkin(fixed=myfit_string, data=pheno_data, kins=kins,id="num",family=gaussian(link="identity"))
    # timeend <- Sys.time()
    # duration2 <- timeend - timestart


    #var.component <- LMM$theta #“kins” is the variance components of the random effects
    coef <- LMM$coefficients #The estimated value of the fixed effect regression coefficient
    ncoef <- length(coef)
    beta <- coef[(ncoef-2):ncoef]
    var.fix <- LMM$cov #The covariance matrix of regression coefficients
    cov <- var.fix[(ncoef-2):ncoef,(ncoef-2):ncoef]
    z <- fnMatSqrtInverse(cov)%*%beta
    nresid <- length(LMM$residuals)
    pA <- pnorm(z,lower.tail = FALSE)
    #pA <- pt(z,df=nresid-ncoef,lower.tail = FALSE)
    pAf <- as.matrix(pA[-3,])
    pAm <- pA[3]
    q1_fisher <- max(-2*log(pchisq(-2*log(prod(pAf)),4,lower.tail = FALSE)*pAm), #QA
                     -2*log(pchisq(-2*log(prod(1-pAf)),4,lower.tail = FALSE)*(1-pAm))) # Qa;

    QXcat_fisher <- 2*pchisq(q1_fisher,4, lower.tail = FALSE) # p value

    # Cauchy method
    # wi <- 1/length(z)
    # t_cauchy <- max(sum(wi*tan(0.5-pA)*pi),sum(wi*tan(0.5-(1-pA))*pi))
    # p_canchy <- 0.5-atan(t_cauchy)/pi
    # QXcat_canchy<- 2*pchisq(t_cauchy,4, lower.tail = FALSE)

    zfsum <- sum(z[1:2])/sqrt(2)
    zm <- z[3]
    w1 <- (2*nf)/(nm+2*nf)  #k=1,No DC
    w2 <- nf/(nm+nf)  #k=2, complete DC
    corr <- sqrt(w1*w2)+sqrt((1-w1)*(1-w2))
    corrp <- matrix(c(1,corr,corr,1),2)
    t1 <- sqrt(w1)*zfsum+sqrt(1-w1)*zm
    t2 <- sqrt(w2)*zfsum+sqrt(1-w2)*zm
    zmaxn_fisher <- max(abs(t1),abs(t2))
    QZmax_fisher <- 1 - mvtnorm::pmvnorm(lower=-rep(zmaxn_fisher,2), upper=rep(zmaxn_fisher,2), corr=corrp)
    res <- c(q1_fisher,zmaxn_fisher,QXcat_fisher,QZmax_fisher)
    names(res) <- c('MQXcat.test','MQZmax.test','MQXcat.p','MQZmax.p')
    return(res)
  }))
  res[,3:4] <- apply(res[,3:4,drop=F],2,function(x){
    x[x>1] <- 1
    return(x)
  })
  SNPnames <- rownames(res)
  if(is.null(SNPnames)){
    SNPnames <- paste('snp',seq(ncol(Genotype)),sep = "")
  }
  Tstat <- data.frame(SNP=SNPnames,res[,1:2,drop=F],row.names = NULL)
  colnames(Tstat) <- c('SNP','MQXcat','MQZmax')
  pvalue <- data.frame(SNP=SNPnames,res[,3:4,drop=F],row.names = NULL)
  colnames(pvalue) <- c('SNP','MQXcat','MQZmax')
  results <- list(Tstat=Tstat,pvalue=pvalue)
  return(results)
}
