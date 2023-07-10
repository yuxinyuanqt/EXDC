#' @title A function for obtaining the confidence interval of a ratio by the delta method
#' @description This function is used to calculate the upper and lower limits of the  delta confidence interval for the ratio estimate.
#' @usage delta(mu_n,mu_d,var_n,var_d,rho=0,con_level=0.95,
#'             truncation_interval=FALSE,truncation_lower=NULL,
#'             truncation_upper=NULL,description=FALSE,H0_test=FALSE,DC_0)
#'
#' @param mu_n  The estimated mean of the numerator.
#' @param mu_d  The estimated mean of the denominator.
#' @param var_n A positive value gives the estimated variance of the numerator.
#' @param var_d A positive value gives the estimated variance of the denominator.
#' @param rho A value between -1 and 1 that represents the estimated correlation coefficient of the numerator and the denominator.
#' @param con_level The confidence level. Should be between 0 and 1. The default is 0.95.
#' @param truncation_interval According to the literature or known knowledge, the ratio should be within a certain interval, for example in the interval from 1 to 2, so we need to truncate the point estimate and the delta's confidence interval into the certain interval to get the final results. The default value is FALSE, that is, no truncation.
#' @param truncation_lower If truncation_interval=TRUE, you need to specify the truncation interval, and truncation_lower indicates the lower limit of the truncation interval. The default value is NULL.
#' @param truncation_upper If truncation_interval=TRUE, you need to specify the truncation interval, and truncation_upper indicates the upper limit of the truncation interval. The default value is NULL.
#' @param description Whether to describe the results of point_estimate and confidence interval in text. The default value is FASLE.
#' @param H0_test Whether to test the null hypothesis of no DC, incomplete DC or complete DC. The default value is FASLE.
#' @param DC_0 If H0_test=TRUE, you need to specify the null hypothesis, where DC_0 takes the values of 1 and 2 for no DC or complete DC, respectively. The default value is NULL.
#'
#' @return a data.frame (H0_test=FALSE) or list (H0_test=TRUE):
#' \describe{
#'    \item{point_estimate}{The point estimate of the ratio.}
#'    \item{Lower_CL}{The lower bound of the estimated interval of the delta method.}
#'    \item{Upper_CL}{The upper bound of the estimated interval of the delta method.}
#'    \item{CI}{The estimated interval of the delta method.}
#'    \item{length_CI}{The length of the estimated confidence interval.}
#'    \item{CI_type}{The CI of the delta method can be divided into 'continuous interval', 'empty set', 'noninformative interval' or 'point'.}
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
#' delta(mu_n,mu_d,var_n,var_d,rho,con_level=0.95,
#'      truncation_interval=FALSE,truncation_lower=NULL,truncation_upper=NULL,
#'      description=FALSE,H0_test=TRUE,DC_0=c(1,2))
#'
#'
#' @references Wang, P.; Zhang, Y.; Wang, BQ; et al. A statistical measure for the skewness of X chromosome inactivation based on case-control design. BMC Bioinformatics 2019, 20, 11.
#'

delta <- function(mu_n,mu_d,var_n,var_d,rho=0,con_level=0.95,truncation_interval=FALSE,
                  truncation_lower=NULL,truncation_upper=NULL,description=FALSE,
                  H0_test=FALSE,DC_0){
  if(!all(sapply(list(mu_n,mu_d,var_n,var_d,rho,con_level,truncation_interval,description),length)==1) &
     !all(sapply(list(truncation_lower,truncation_upper),length) %in% c(0,1))){
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
  if(con_level<=0 | con_level>=1) stop("Confidence level should be between 0 and 1!")
  if(truncation_interval==TRUE){
    if(is.null(truncation_lower) | is.null(truncation_upper)){
      stop("You need to specify the lower and upper bounds of the truncation interval!")
    } else if(truncation_lower>=truncation_upper){
      stop("'truncation_upper' should be greater than 'truncation_lower'!")
    }
  }
  if(mu_d==0){
    warning("The 'mu_d' is 0, the delta method cannot produce meaningful interval estimate, please change the order of the 'mu_n' and 'mu_d', or use the Fieller method to estimate the confidence interval.")
    results <- data.frame(point_estimate=mu_n/mu_d,
                          Lower_CL=NaN,
                          Upper_CL=NaN,
                          CI=NaN,
                          length_CI=NaN,
                          CI_type=NaN)
    return(results)
  } else {
    alpha <- 1-con_level
    dhat <- mu_n/mu_d #the point estimate of DC
    cov_nd <- rho*sqrt(var_n)*sqrt(var_d) #the estimated covariance of the numerator and the denominator
    #calculate the (1 − alpha/2) upper quantile of a standard normal distribution
    z_alpha <- qnorm(1-alpha/2) #standard normal distribution

    #the estimated SE of dhat
    se_dhat <- sqrt(var_n/mu_d^2+var_d*mu_n^2/mu_d^4-(2*mu_n*cov_nd)/mu_d^3)

    Upper_origin <- dhat+z_alpha*se_dhat
    Lower_origin <- dhat-z_alpha*se_dhat
    length_origin <- Upper_origin-Lower_origin
    if(is.finite(Lower_origin) | is.finite(Upper_origin)){
      CI_type_origin <- 'continuous interval'
      CI_origin <- ifelse(is.finite(Lower_origin) & is.finite(Upper_origin),
                          paste('[',Lower_origin,', ',Upper_origin,']',sep=''),
                          ifelse(is.infinite(Lower_origin),
                                 paste('(',Lower_origin,', ',Upper_origin,']',sep=''),
                                 paste('[',Lower_origin,', ',Upper_origin,')',sep='')))

    } else{
      CI_type_origin <- 'noninformative interval'
      CI_origin <- '(-Inf, Inf)'
    }


    if(truncation_interval==TRUE){
      #truncate the origin point estimate
      dhat_delta <- ifelse(dhat>=truncation_lower & dhat<=truncation_upper,dhat,
                             ifelse(dhat<truncation_lower,truncation_lower,truncation_upper))

      if(CI_type_origin=='continuous interval'){
        if(Lower_origin>truncation_upper | Upper_origin<truncation_lower){
          Upper_delta <- Lower_delta <- NA
          length_delta <- 0
          CI_type <- 'empty set'
          CI_delta <- 'empty set'
        } else {
          Upper_delta <- min(Upper_origin,truncation_upper)
          Lower_delta <- max(Lower_origin,truncation_lower)
          length_delta <- Upper_delta-Lower_delta
          CI_type <- 'continuous interval'
          CI_delta <- paste('[',Lower_delta,', ',Upper_delta,']',sep='')
          if(length_delta==0){ #Lower_delta=Upper_delta
            CI_type <- 'point'
            CI_delta <- paste('{',Lower_delta,'}',sep='')
          } else if(Upper_delta==truncation_upper & Lower_delta==truncation_lower){
            CI_type <- 'noninformative interval'
          }
        }
      } else {
        Lower_delta <- truncation_lower
        Upper_delta <- truncation_upper
        length_delta <- Upper_delta-Lower_delta
        CI_type <- 'noninformative interval'
        CI_delta <- paste('[',Lower_delta,', ',Upper_delta,']',sep='')
      }
    }

    if(truncation_interval==FALSE){
      if(description==TRUE){
        cat("The point estimate of DC is ",dhat, ".\n",sep='')
        cat("The confidence interval of DC is ",CI_origin, ".\n\n",sep='')
      }
      results <- data.frame(point_estimate=dhat,Lower_CL=Lower_origin,
                            Upper_CL=Upper_origin,CI=CI_origin,
                            length_CI=length_origin,CI_type=CI_type_origin)
      rownames(results) <- 'origin'
    } else{
      if(description==TRUE){
        cat("The origin point estimate of DC is ",dhat, ".\n",sep='')
        cat("The origin confidence interval of DC is ",CI_origin, ".\n",sep='')
        cat("The truncated point estimate of DC is ",dhat_delta, ".\n",sep='')
        cat("The truncated confidence interval of DC is ",CI_delta, ".\n\n",sep='')
      }

      results <- data.frame(point_estimate=c(dhat,dhat_delta),
                            Lower_CL=c(Lower_origin,Lower_delta),
                            Upper_CL=c(Upper_origin,Upper_delta),
                            CI=c(CI_origin,CI_delta),
                            length_CI=c(length_origin,length_delta),
                            CI_type=c(CI_type_origin,CI_type))
      rownames(results) <- c('origin','truncation')
    }

    if(H0_test==TRUE){
      Z_delta <- (dhat-DC_0)/se_dhat
      P_delta <- 2*pnorm(abs(Z_delta),0,1,lower.tail = FALSE)  # P-value for the Fieller's method
      H0_test <- data.frame(Z_score=Z_delta,P_value=P_delta)
      row.names(H0_test) <- paste('DC_0=',DC_0,sep='')
      results <- list(Estimation_result=results,
                      H0_test=H0_test)
    }

    return(results)
  }

}
