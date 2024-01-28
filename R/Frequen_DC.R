#' @title A function for estimating the degree of dosage compensation at X-linked loci using the delta, Fieller's or penalized Fieller's method
#'
#' @description
#' This function is used to calculate the upper and lower limits of the delta, Fieller's or penalized Fieller's confidence interval for the ratio estimate.
#'
#' @usage
#' Frequen_DC(SNP,Y,trait_type='quantitative',Sex,Covariate=NULL,
#'            variance_heterogeneity=FALSE,MGC_Cutoff=20,method='Fieller',
#'            df_n=NULL,df_d=NULL,con_level=0.95,truncation_interval=FALSE,
#'            truncation_lower=NULL,truncation_upper=NULL,description=FALSE,
#'            H0_test=FALSE,DC_0)
#'
#'
#' @param SNP Each genotype is coded as 0, 1 or 2 for females and coded as 0 or 1 for males, indicating the number of reference allele. The length of \code{SNP} should be the same to those of \code{Y}, \code{Sex} and \code{Covariate}.
#' @param Y A numeric vector of a quantitative trait, such as human height.
#' @param trait_type A character string either being "quantitative" or "qualitative", the default value is "quantitative".
#' @param Sex A vector of the genetic sex following PLINK default coding, where males are coded as 1 and females are coded as 2.
#' @param Covariate Optional: a vector or a matrix of covariates, such as age and BMI.
#' @param variance_heterogeneity For quantitative traits, we empolyed weighted least squares to fit the model when variance_heterogeneity is set to TRUE, and ordinary least squares when variance_heterogeneity is set to FASLE.
#' @param MGC_Cutoff Cutoff for the minimum genotype count in either females or males (default=20). When the minimum genotype count for a SNP is below this specified cutoff, we will employ ordinary least squares instead of weighted least squares for genetic association analysis with quantitative traits. This is mainly to avoid inflated type I errors caused by the instability of the weights. (Yang et al., 2022; Deng et al., 2019; Soave et al., 2015).
#' @param method a character string indicating which kind of estimation methods is to be used. One of "Fieller"(default), "PenFieller" , "delta" or "all": can be abbreviated.
#' @param df_n  The degree of freedom of the numerator. The default value is NULL.
#' @param df_d  The degree of freedom of the denominator. The default value is NULL.
#' @param con_level The confidence level. Should be between 0 and 1. The default is 0.95.
#' @param truncation_interval According to the literature or known knowledge, the ratio should be within a certain interval, for example in the interval from 1 to 2, so we need to truncate the point estimate and the Fieller's confidence interval into the certain interval to get the final results. The default value is FALSE, that is, no truncation.
#' @param truncation_lower If truncation_interval=TRUE, you need to specify the truncation interval, and truncation_lower indicates the lower limit of the truncation interval. The default value is NULL.
#' @param truncation_upper If truncation_interval=TRUE, you need to specify the truncation interval, and truncation_upper indicates the upper limit of the truncation interval. The default value is NULL.
#' @param description Whether to describe the results of point_estimate and confidence interval in text. The default value is FASLE.
#' @param H0_test Whether to test the null hypothesis of no DC, incomplete DC or complete DC. The default value is FASLE.
#' @param DC_0 If H0_test=TRUE, you need to specify the null hypothesis, where DC_0 takes the values of 1 and 2 for no DC or complete DC, respectively. The default value is NULL.
#'
#' @return a list:
#' \describe{
#'    \item{point_estimate}{The point estimate of the ratio.}
#'    \item{Lower_CL}{The lower bound of the estimated interval of the Fieller's method.}
#'    \item{Upper_CL}{The upper bound of the estimated interval of the Fieller's method.}
#'    \item{CI}{The estimated interval of the Fieller's method.}
#'    \item{length_CI}{The length of the estimated confidence interval.}
#'    \item{CI_type}{The CI of the Fieller's method can be divided into 'continuous interval','discontinuous interval', 'empty set', 'noninformative interval' or 'point'.}
#'    \item{H0_test}{If H0_test=TRUE, returns the Z-score test statistic and corresponding p-value under the null hypothesis.}
#' }
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' data("Graves_Meta_analysis")
#' Frequen_DC(Graves_Meta_analysis$rs3827440,Graves_Meta_analysis$Graves_disease,
#'            trait_type='qualitative',Graves_Meta_analysis$Sex,method='all',
#'            con_level=0.95,truncation_interval=TRUE,truncation_lower=1,truncation_upper=2,
#'            description=TRUE,H0_test=TRUE,DC_0=2)
#'
#' @references Yang ZY, Liu W, Yuan YX, et al. Robust association tests for quantitative traits on the X chromosome. \emph{Heredity}, 2022, \strong{129}: 244–256.
#' @references Deng WQ, Mao S, Kalnapenkis A, et al. Analytical strategies to include the X-chromosome in variance heterogeneity analyses: evidence for trait-specific polygenic variance structure. \emph{Genetic Epidemiology}, 2019, \strong{43}: 815-830.
#' @references Soave D, Corvol H, Panjwani N, et al. A joint location-scale test improves power to detect associated SNPs, gene sets, and pathways. \emph{The American Journal of Human Genetics}, 2015, \strong{97}: 125–138.
#'
#'
Frequen_DC <- function(SNP,Y,trait_type='quantitative',Sex,Covariate=NULL,
                       variance_heterogeneity=FALSE,MGC_Cutoff=20,method='Fieller',
                       df_n=NULL,df_d=NULL,con_level=0.95,truncation_interval=FALSE,
                       truncation_lower=NULL,truncation_upper=NULL,description=FALSE,
                       H0_test=FALSE,DC_0){
  if (missing(SNP)){
    stop("The SNP input is missing.")
  }
  if (missing(Y)){
    stop("The trait input is missing.")
  }
  if (missing(Sex)){
    stop("The Sex input is missing.")
  }
  if(!all(unique(Sex) %in% c(1,2))){
    stop('Sex must be a vector of 1(males) and 2(females)')
  }
  if(!all(c(is.vector(SNP),is.vector(Y),is.vector(Sex)))){
    stop("'SNP', 'Y' and 'Sex' must be vectors!")
  }
  if(length(table(Sex))==1){
    stop("Only Males or Females detected")
  }
  Sex[Sex==2] <- 0
  if(!is.null(Covariate) & (is.matrix(Covariate)|
                            is.vector(Covariate)|
                            is.data.frame(Covariate))){
    Covariate <- as.data.frame(Covariate)
    if(length(unique(c(length(Y),nrow(Covariate),length(SNP),length(Sex))))!=1){
      stop("Make sure the inputs have the same length!")
    }
    Phedata <- cbind(Y,Covariate)
    Null_Model_string <- paste('Y',paste(colnames(Covariate), collapse = "+"),sep = "~")
  }else if(is.null(Covariate)){
    if(length(unique(c(length(Y),length(SNP),length(Sex))))!=1){
      stop("Make sure the inputs have the same length!")
    }
    Phedata <- as.data.frame(Y)
    Null_Model_string <- 'Y~1'
  }else{
    stop("Covariate must be a vector or a matrix!")
  }

  myfit_string <- paste(Null_Model_string,'+Sex+SNP_AAf+SNP_Af+SNP_m',sep = "")
  SNP_AAf <- ifelse(SNP==2 & Sex==0,1,0)
  SNP_Af <- ifelse(SNP!=0 & Sex==0,1,0)
  SNP_m <- ifelse(SNP!=0 & Sex==1,1,0)
  if(trait_type=='quantitative'){
    myfit0 <- lm(as.formula(myfit_string),na.action = na.exclude,data = Phedata)
    ncoef <- length(myfit0$coefficients)
    beta_hat <- coefficients(myfit0)[(ncoef-2):ncoef]
    cov_hat <- vcov(myfit0)[(ncoef-2):ncoef,(ncoef-2):ncoef]

    if(variance_heterogeneity==TRUE){
      #Males have two genotypes, females have three genotypes
      SEX_SNP <- factor(interaction(Sex,SNP))

      #Grouping (five groups) to calculate the inverse of the residual variance
      inv_var <- tapply(resid(myfit0),SEX_SNP,function(x){1/var(x,na.rm=T)})
      if(length(inv_var)!=5 | any(table(inv_var)<MGC_Cutoff)){
        #When there are no five groups or the minimum genotype count is less than MGC_Cutoff, we use the ordinary least squares regression
        myfit <- myfit0
      } else {
        names_inv_var <- names(inv_var)
        w_fit <- ifelse(SEX_SNP==names_inv_var[1],inv_var[1],
                        ifelse(SEX_SNP==names_inv_var[2],inv_var[2],
                               ifelse(SEX_SNP==names_inv_var[3],inv_var[3],
                                      ifelse(SEX_SNP==names_inv_var[4],inv_var[4],inv_var[5]))))
        #the weighted least square method
        myfit <- lm(as.formula(myfit_string),weights=w_fit,na.action = na.exclude,data = Phedata)
      }
      beta_hat <- coefficients(myfit)[(ncoef-2):ncoef]
      cov_hat <- vcov(myfit)[(ncoef-2):ncoef,(ncoef-2):ncoef]
    }
  } else if(trait_type=='qualitative'){
    myfit <- glm(as.formula(myfit_string),family = binomial(link = "logit"),
                 na.action = na.exclude,data = Phedata)
    ncoef <- length(myfit$coefficients)
    beta_hat <- coefficients(myfit)[(ncoef-2):ncoef]
    cov_hat <- vcov(myfit)[(ncoef-2):ncoef,(ncoef-2):ncoef]
  } else {
    stop("trait_type should be 'quantitative' or 'qualitative'")
  }

  mu_n <- beta_hat[3]
  mu_d <- (beta_hat[1]+beta_hat[2])/2
  var_n <- cov_hat[3,3]
  var_d <- (cov_hat[1,1]+cov_hat[2,2]+2*cov_hat[1,2])/4
  cov_nd <- (cov_hat[1,3]+cov_hat[2,3])/2
  rho <- cov_nd/sqrt(var_n*var_d)
  parameters <- c(mu_n,mu_d,var_n,var_d,cov_nd,rho)
  if(any(is.na(parameters))| any(is.infinite(parameters))){
    stop("Can't obtain the estimate due to missing or infinite values!")
  }

  if(method=='Fieller'){
    DC_Fieller <- Fieller(mu_n,mu_d,var_n,var_d,rho,
                          df_n,df_d,con_level,truncation_interval,
                          truncation_lower,truncation_upper,description,
                          H0_test,DC_0)
    return(DC_Fieller)
  } else if(method=='PenFieller'){
    DC_PenFieller <- PenFieller(mu_n,mu_d,var_n,var_d,rho,
                                df_n,df_d,con_level,truncation_interval,
                                truncation_lower,truncation_upper,description,
                                H0_test,DC_0)
    return(DC_PenFieller)
  } else if(method=='delta'){
    DC_delta <- delta(mu_n,mu_d,var_n,var_d,rho,
                      con_level,truncation_interval,
                      truncation_lower,truncation_upper,description,
                      H0_test,DC_0)
    return(DC_delta)
  } else if(method=='all'){
    DC_Fieller <- Fieller(mu_n,mu_d,var_n,var_d,rho,
                          df_n,df_d,con_level,truncation_interval,
                          truncation_lower,truncation_upper,description=FALSE,
                          H0_test,DC_0)
    DC_PenFieller <- PenFieller(mu_n,mu_d,var_n,var_d,rho,
                                df_n,df_d,con_level,truncation_interval,
                                truncation_lower,truncation_upper,description=FALSE,
                                H0_test,DC_0)
    DC_delta <- delta(mu_n,mu_d,var_n,var_d,rho,
                      con_level,truncation_interval,
                      truncation_lower,truncation_upper,description=FALSE,
                      H0_test,DC_0)
    return(list(DC_Fieller=DC_Fieller,DC_PenFieller=DC_PenFieller,DC_delta=DC_delta))
  } else {
    stop('method must be one of "Fieller" (default), "PenFieller", "delta" or "all" (showing the results of the above three methods)!')
  }
}
