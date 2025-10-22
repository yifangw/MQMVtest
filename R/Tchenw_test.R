###########################################################
# mean-based method: Tchenw for unrelated individuals
###########################################################

#' @title mean-based method: Tchenw for unrelated individuals, including a dominant effect of heterozygous (Aa) females
#' Similar to the MTchen method, Tchenw is an approach designed for unrelated individuals and
#' therefore does not incorporate a kinship matrix to account for correlations between individuals.
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
#'
#' @return The p value and test statistic of Tchenw.
#' @export
#'
#' @examples Tchenw_test(Genotype,Y,Sex,
#'                       Covariate=unrelated[,"Age"],
#'                       missing_cutoff=0.15,
#'                       MAF_Cutoff=NULL,
#'                       MGC_Cutoff=20)
Tchenw_test <- function(Genotype,Y,Sex,Covariate=NULL, missing_cutoff=0.15,MAF_Cutoff=NULL,MGC_Cutoff=20){

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
  Sex[Sex==2] <- 0

  # handle Covariates
  if(!is.null(Covariate) & (is.matrix(Covariate)|
                            is.vector(Covariate)|
                            is.data.frame(Covariate))){
    Covariate <- as.data.frame(Covariate)
    if(length(unique(c(length(Y), nrow(Covariate), nrow(Genotype), length(Sex)))) != 1){
      stop("Make sure the inputs have the same length!")
    }
    Phedata <- cbind(Y, Covariate)

    covariate_names <- colnames(Covariate)
    Null_Model_string <- paste('Y ~', paste(covariate_names, collapse = " + "))
  } else if(is.null(Covariate)){
    Phedata <- as.data.frame(Y)
    if(length(unique(c(length(Y), nrow(Genotype), length(Sex)))) != 1){
      stop("Make sure the inputs have the same length!")
    }
    Null_Model_string <- 'Y ~ 1'
  } else {
    stop("Covariate must be a vector or a matrix!")
  }

  if(is.vector(Genotype)){
    Genotype <- as.matrix(Genotype)
  }

  if(is.null(MAF_Cutoff)){
    MAF_Cutoff <- 1/sqrt(2*length(Sex))
  }
  MAF_w <- ifelse(Sex==0,0.5,1)
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


  myfit_string <- paste(Null_Model_string, '+ snp*Sex + snp_D', sep="")


  myfit_string2 <- paste(Null_Model_string, '+ Sex', sep="")

  res <- t(apply(Genotype, 2, function(snp){

    if(is.null(Covariate)) {
      complete_data <- na.omit(data.frame(
        Y = Y,
        snp = snp,
        Sex = Sex
      ))
    } else {

      complete_data <- na.omit(data.frame(
        Y = Y,
        snp = snp,
        Sex = Sex,
        Covariate
      ))
    }

    # Calculate the dominant effect
    complete_data$snp_D <- ifelse(complete_data$Sex == 0 & complete_data$snp == 1, 1, 0)


    cluster <- rep(0, nrow(complete_data))
    cluster[complete_data$Sex == 0 & complete_data$snp == 0] <- 1
    cluster[complete_data$Sex == 0 & complete_data$snp != 0 & complete_data$snp != 2] <- 2
    cluster[complete_data$Sex == 0 & complete_data$snp == 2] <- 3
    cluster[complete_data$Sex == 1 & complete_data$snp == 0] <- 4
    cluster[complete_data$Sex == 1 & complete_data$snp != 0] <- 5

    complete_data$cluster <- cluster


    if(is.null(Covariate)) {
      full_formula <- as.formula("Y ~ snp*Sex + snp_D")
      null_formula <- as.formula("Y ~ Sex")
    } else {

      cov_terms <- paste(covariate_names, collapse = " + ")
      full_formula <- as.formula(paste("Y ~ snp*Sex + snp_D +", cov_terms))
      null_formula <- as.formula(paste("Y ~ Sex +", cov_terms))
    }

    # fit LMM
    lmm_Tchenw <- lm(full_formula, data = complete_data)

    # weights
    inv_var_Tchenw <- tapply(resid(lmm_Tchenw),
                           complete_data$cluster,
                           function(x){1/var(x, na.rm=TRUE)})

    # Handle the NA weight
    if(any(is.na(inv_var_Tchenw))) {
      overall_var <- var(resid(lmm_Tchenw), na.rm=TRUE)
      inv_var_Tchenw[is.na(inv_var_Tchenw)] <- 1/overall_var
    }

    # Assign weights to each observation
    w_Tchenw <- numeric(nrow(complete_data))
    for(i in 1:5) {
      w_Tchenw[complete_data$cluster == i] <- inv_var_Tchenw[i]
    }

    # LRT model
    lm_weighted <- lm(full_formula, weights = w_Tchenw, data = complete_data)
    lm_null_weighted <- lm(null_formula, weights = w_Tchenw, data = complete_data)

    # conduct LRT
    if(requireNamespace("lmtest", quietly = TRUE)) {
      Tchenw_p <- lmtest::lrtest(lm_weighted, lm_null_weighted)$Pr[2]
    } else {
      # Using anova to perform LRT
      anova_result <- anova(lm_null_weighted, lm_weighted, test="LRT")
      Tchenw_p <- anova_result$`Pr(>Chi)`[2]
    }

    res <- c(NA, Tchenw_p)
    names(res) <- c('Tchenw.test', 'Tchenw_p')
    return(res)
  }))

  SNPnames <- rownames(res)
  if(is.null(SNPnames)){
    SNPnames <- paste('snp', seq(ncol(Genotype)), sep="")
  }

  Tstat <- data.frame(SNP=SNPnames, res[,1,drop=F], row.names=NULL)
  colnames(Tstat) <- c('SNP', 'Tchenw')

  pvalue <- data.frame(SNP=SNPnames, res[,2,drop=F], row.names=NULL)
  colnames(pvalue) <- c('SNP', 'Tchenw')

  results <- list(Tstat=Tstat, pvalue=pvalue)
  return(results)
}
