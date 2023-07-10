#' @title A function for obtaining the confidence interval of a ratio by the penalized Fieller's method
#' @description This function is used to calculate the upper and lower limits of the penalized Fieller's confidence interval for the ratio estimate.
#' @usage PenFieller(mu_n,mu_d,var_n,var_d,rho=0,df_n=NULL,df_d=NULL,con_level=0.95,
#'                  truncation_interval=FALSE,truncation_lower=NULL,
#'                  truncation_upper=NULL,description=FALSE,H0_test=FALSE,DC_0)
#'
#' @param mu_n  The estimated mean of the numerator.
#' @param mu_d  The estimated mean of the denominator.
#' @param var_n A positive value gives the estimated variance of the numerator.
#' @param var_d A positive value gives the estimated variance of the denominator.
#' @param rho A value between -1 and 1 that represents the estimated correlation coefficient of the numerator and the denominator.
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
#' @return a data.frame (H0_test=FALSE) or list (H0_test=TRUE):
#' \describe{
#'    \item{point_estimate}{The point estimate of the ratio.}
#'    \item{Lower_CL}{The lower bound of the estimated interval of the penalized Fieller's method.}
#'    \item{Upper_CL}{The upper bound of the estimated interval of the penalized Fieller's method.}
#'    \item{CI}{The estimated interval of the penalized Fieller's method.}
#'    \item{length_CI}{The length of the estimated confidence interval.}
#'    \item{CI_type}{The CI of the penalized Fieller's method can be divided into 'continuous interval', 'empty set', 'noninformative interval' or 'point'.}
#'    \item{H0_test}{If H0_test=TRUE, returns the Z-score test statistic and corresponding p-value under the null hypothesis.}
#' }
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' data(Graves_Meta_analysis)
#' attach(Graves_Meta_analysis)
#' SNP_AAf <- ifelse(rs3827440==2 & Sex==2,1,0)
#' SNP_Af <- ifelse(rs3827440!=0 & Sex==2,1,0)
#' SNP_m <- ifelse(rs3827440!=0 & Sex==1,1,0)
#' myfit <- glm(Graves_disease~Sex+SNP_AAf+SNP_Af+SNP_m,
#'             family = binomial(link = "logit"),
#'             na.action = na.exclude)
#' detach(Graves_Meta_analysis)
#' ncoef <- length(myfit$coefficients)
#' beta_hat <- coefficients(myfit)[(ncoef-2):ncoef]
#' cov_hat <- vcov(myfit)[(ncoef-2):ncoef,(ncoef-2):ncoef]
#'
#' mu_n <- beta_hat[3]
#' mu_d <- (beta_hat[1]+beta_hat[2])/2
#' var_n <- cov_hat[3,3]
#' var_d <- (cov_hat[1,1]+cov_hat[2,2]+2*cov_hat[1,2])/4
#' cov_nd <- (cov_hat[1,3]+cov_hat[2,3])/2
#' rho <- cov_nd/sqrt(var_n*var_d)
#'
#' PenFieller(mu_n,mu_d,var_n,var_d,rho,df_n=NULL,df_d=NULL,con_level=0.95,
#'            truncation_interval=FALSE,truncation_lower=NULL,truncation_upper=NULL,
#'            description=FALSE,H0_test=TRUE,DC_0=c(1,2))
#'
#' @references Wang, P.; Xu, S.; Wang, Y. X.; et al. Penalized Fieller's confidence interval for the ratio of bivariate normal means. Biometrics 2021, 77, 1355-1368.
#'
#'
PenFieller <- function(mu_n,mu_d,var_n,var_d,rho=0,df_n=NULL,df_d=NULL,con_level=0.95,
                       truncation_interval=FALSE,truncation_lower=NULL,truncation_upper=NULL,description=FALSE,
                       H0_test=FALSE,DC_0){
  if(!all(sapply(list(mu_n,mu_d,var_n,var_d,rho,con_level,truncation_interval,description),length)==1) &
     !all(sapply(list(df_n,df_d,truncation_lower,truncation_upper),length) %in% c(0,1))){
    stop('All input variables should be scalars!')
  }
  if(any(c(class(mu_n),class(mu_d),class(var_n),class(var_d),class(rho))!="numeric")){
    mu_n <- as.numeric(mu_n);mu_d <- as.numeric(mu_d)
    var_n <- as.numeric(var_n);var_d <- as.numeric(var_d);rho <- as.numeric(rho)
  }
  if(any(is.na(c(mu_n,mu_d,var_n,var_d,rho)))){
    stop('The value of mu_n,mu_d,var_n,var_d or rho is missing!')
  }
  if(var_n<=0 | var_d<=0) stop("Variance should be a positive number!")
  if(abs(rho)>=1) stop("Correlation coefficient should be between -1 and 1!")
  if(any(is.na(c(df_n,df_d))) |
     any(is.infinite(c(df_n,df_d))) |
     ((!inherits(df_n,'NULL')) && (df_n%%1!=0 || df_n<=0)) |
     ((!inherits(df_d,'NULL')) && (df_d%%1!=0 || df_d<=0))){
    stop("Degree of freedom should be a positive integer!")
  }
  if(con_level<=0 | con_level>=1) stop("Confidence level should be between 0 and 1!")
  if(truncation_interval==TRUE){
    if(is.null(truncation_lower) | is.null(truncation_upper)){
      stop("You need to specify the lower and upper bounds of the truncation interval!")
    } else if(truncation_lower>=truncation_upper){
      stop("'truncation_upper' should be greater than 'truncation_lower'!")
    }
  }

  if(mu_d==0){
    warning("The 'mu_d' is 0, please change the order of the 'mu_n' and 'mu_d', or use the Fieller method to estimate the confidence interval. The estimated results of the Fieller method will be shown here.")
    results <- Fieller(mu_n,mu_d,var_n,var_d,rho,
                       df_n,df_d,con_level,truncation_interval,
                       truncation_lower,truncation_upper,description)
    return(results)
  } else {
    alpha <- 1-con_level
    dhat <- mu_n/mu_d #the point estimate of DC
    cov_nd <- rho*sqrt(var_n)*sqrt(var_d) #the estimated covariance of the numerator and the denominator
    #calculate the (1 − alpha/2) upper quantile of a standard normal distribution or t distribution
    if (inherits(df_n,'NULL')||inherits(df_d,'NULL')) {
      t_alpha <- qnorm(1-alpha/2) #standard normal distribution
    } else {
      if (df_n==df_d) {
        df <- df_n
      } else {
        #correct for the degrees of freedom
        df <- (var_n+dhat^2*var_d)^2/(var_n^2/df_n+dhat^2*var_d^2/df_d)
      }
      #t distribution
      t_alpha <- qt(1-alpha/2,df)
    }
    c_p <- t_alpha^2/4
    p_mu_d <- mu_d/2+sign(mu_d)*sqrt(mu_d^2/4+c_p*var_d)
    w <- p_mu_d/(2*p_mu_d-mu_d)
    dt <- mu_n/p_mu_d
    p_mu_n <- w^(-1)*mu_n
    p_var_d <- w^2*var_d
    p_var_n <- w^(-2)*var_n-4*(w^(-1)-1)*dt*cov_nd+4*(1-w)^2*dt^2*var_d
    p_cov_nd <- cov_nd-2*w*(1-w)*dt*var_d

    results <- Fieller(p_mu_n,p_mu_d,p_var_n,p_var_d,rho=p_cov_nd/sqrt(p_var_n*p_var_d),
                       df_n=df_n,df_d=df_d,con_level=con_level,truncation_interval=truncation_interval,
                       truncation_lower=truncation_lower,truncation_upper=truncation_upper,description=description,
                       H0_test=H0_test,DC_0=DC_0)
    return(results)
  }

}
