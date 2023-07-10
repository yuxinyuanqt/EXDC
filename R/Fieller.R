#' @title A function for obtaining the confidence interval of a ratio by the Fieller's method
#' @description This function is used to calculate the upper and lower limits of the Fieller's confidence interval for the ratio estimate.
#' @usage Fieller(mu_n,mu_d,var_n,var_d,rho=0,df_n=NULL,df_d=NULL,con_level=0.95,
#'               truncation_interval=FALSE,truncation_lower=NULL,
#'               truncation_upper=NULL,description=FALSE,H0_test=FALSE,DC_0)
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
#' Fieller(mu_n,mu_d,var_n,var_d,rho,df_n=NULL,df_d=NULL,con_level=0.95,
#'        truncation_interval=FALSE,truncation_lower=NULL,truncation_upper=NULL,
#'        description=FALSE,H0_test=TRUE,DC_0=c(1,2))
#'
#' @references Wang, P.; Zhang, Y.; Wang, BQ; et al. A statistical measure for the skewness of X chromosome inactivation based on case-control design. BMC Bioinformatics 2019, 20, 11.
#'
#'
Fieller <- function(mu_n,mu_d,var_n,var_d,rho=0,df_n=NULL,df_d=NULL,con_level=0.95,
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

  #The numbers A_Fieller, B_Fieller and C_Fieller are the coefficients of the equation and may be distinguished by respectively calling them, the quadratic coefficient, the linear coefficient and the constant coefficient or free term.
  #A_Fieller*DC^2+B_Fieller*DC+C_Fieller<=0
  A_Fieller <- mu_d^2-t_alpha^2*var_d
  B_Fieller <- 2*(cov_nd*t_alpha^2-mu_n*mu_d)
  C_Fieller <- mu_n^2-t_alpha^2*var_n
  #delta_Fieller: the discriminant of the quadratic equation
  delta_Fieller <- B_Fieller^2-4*A_Fieller*C_Fieller
  CI_type_origin <- NULL

  #According to the signs of A_Fieller and delta_Fieller, the CI of the Fieller’s method can be divided into 'continuous interval', 'discontinuous interval', 'empty set', 'noninformative interval' or 'point'
  if(A_Fieller!=0){
    if(delta_Fieller>0){
      if(A_Fieller>0){
        Upper_origin <- (-B_Fieller+sqrt(delta_Fieller))/(2*A_Fieller)
        Lower_origin <- (-B_Fieller-sqrt(delta_Fieller))/(2*A_Fieller)
        length_origin <- Upper_origin-Lower_origin
        CI_type_origin <- 'continuous interval'
        CI_origin <- paste('[',Lower_origin,', ',Upper_origin,']',sep='')
      } else {
        Upper_origin <- (-B_Fieller-sqrt(delta_Fieller))/(2*A_Fieller)
        Lower_origin <- (-B_Fieller+sqrt(delta_Fieller))/(2*A_Fieller)
        length_origin <- Inf
        CI_type_origin <- 'discontinuous interval'
        CI_origin <- paste('(-Inf, ',Lower_origin,'] U [',Upper_origin,', Inf)',sep='')
      }
    } else if(delta_Fieller<0){
      if(A_Fieller>0){
        Upper_origin <- Lower_origin <- NA
        length_origin <- 0
        CI_type_origin <- 'empty set'
        CI_origin <- 'empty set'
      } else {
        Upper_origin <- Inf
        Lower_origin <- -Inf
        length_origin <- Inf
        CI_type_origin <- 'noninformative interval'
        CI_origin <- '(-Inf, Inf)'
      }
    } else { #delta_Fieller=0
      if(A_Fieller>0){
        Upper_origin <- Lower_origin <- -B_Fieller/(2*A_Fieller)
        length_origin <- 0
        CI_type_origin <- 'point'
        CI_origin <- paste('{',Lower_origin,'}',sep='')
      } else if(A_Fieller<0){
        Upper_origin <- Inf
        Lower_origin <- -Inf
        length_origin <- Inf
        CI_type_origin <- 'noninformative interval'
        CI_origin <- '(-Inf, Inf)'
      }
    }
  } else { #A_Fieller=0
    #the quadratic equation degenerates to be B_Fieller*DC+C_Fieller<=0
    if(B_Fieller<0){
      Lower_origin <- -C_Fieller/B_Fieller
      Upper_origin <- Inf
      CI_origin <- paste('[',Lower_origin,', ',Upper_origin,')',sep='')
    } else if(B_Fieller>0){
      Lower_origin <- -Inf
      Upper_origin <- -C_Fieller/B_Fieller
      CI_origin <- paste('(',Lower_origin,', ',Upper_origin,']',sep='')
    } else{
      stop('Both the quadratic coefficient and the linear coefficient of the quadratic equation are 0!')
    }
    length_origin <- Inf
    CI_type_origin <- 'continuous interval' #(-Inf, a] or [b, Inf)
  }

  if(truncation_interval==TRUE){
    #truncate the origin point estimate
    dhat_Fieller <- ifelse(dhat>=truncation_lower & dhat<=truncation_upper,dhat,
                           ifelse(dhat<truncation_lower,truncation_lower,truncation_upper))

    if(CI_type_origin=='continuous interval'){
      if(Lower_origin>truncation_upper | Upper_origin<truncation_lower){
        Upper_Fieller <- Lower_Fieller <- NA
        length_Fieller <- 0
        CI_type <- 'empty set'
        CI_Fieller <- 'empty set'
      } else {
        Upper_Fieller <- min(Upper_origin,truncation_upper)
        Lower_Fieller <- max(Lower_origin,truncation_lower)
        length_Fieller <- Upper_Fieller-Lower_Fieller
        CI_type <- 'continuous interval'
        if(Upper_Fieller==truncation_upper & Lower_Fieller==truncation_lower){
          CI_type <- 'noninformative interval'
        }
        CI_Fieller <- paste('[',Lower_Fieller,', ',Upper_Fieller,']',sep='')
        if(length_Fieller==0){ #Lower_Fieller=Upper_Fieller
          CI_type <- 'point'
          CI_Fieller <- paste('{',Lower_Fieller,'}',sep='')
        }
      }
    } else if(CI_type_origin=='discontinuous interval'){
      if(Lower_origin<truncation_lower & Upper_origin>truncation_upper){ #no intersection
        Upper_Fieller <- Lower_Fieller <- NA
        length_Fieller <- 0
        CI_type <- 'empty set'
        CI_Fieller <- 'empty set'
      } else if(Lower_origin>=truncation_lower & Upper_origin<=truncation_upper){
        Lower_Fieller <- Lower_origin
        Upper_Fieller <- Upper_origin
        length_Fieller <- (Lower_Fieller-truncation_lower)+(truncation_upper-Upper_Fieller)
        if(Lower_Fieller==truncation_lower & truncation_upper==Upper_Fieller){
          CI_type <- 'point' #two point
          CI_Fieller <- paste('{',Lower_Fieller,', ',Upper_Fieller,'}',sep='')
        } else if(Lower_Fieller==truncation_lower){
          CI_type <- 'discontinuous interval' #{a} U [b, c]
          CI_Fieller <- paste('{',Lower_Fieller,'} U [',
                              Upper_Fieller,', ',truncation_upper,']',sep='')
        } else if(truncation_upper==Upper_Fieller){
          CI_type <- 'discontinuous interval' #[a,b] U {c}
          CI_Fieller <- paste('[',truncation_lower,', ',Lower_Fieller,'] U {',
                              Upper_Fieller,'}',sep='')
        } else {
          #[a,b] U [c,d]
          CI_type <- 'discontinuous interval' #truncated discontinuous interval
          CI_Fieller <- paste('[',truncation_lower,', ',Lower_Fieller,'] U [',
                              Upper_Fieller,', ',truncation_upper,']',sep='')
        }
      } else{ #After truncating, it becomes a continuous interval or a point
        if(Lower_origin >= truncation_lower){ #Upper_origin>truncation_upper
          Lower_Fieller <- truncation_lower
          Upper_Fieller <- min(Lower_origin,truncation_upper)
        } else{
          #Lower_origin < truncation_lower & Upper_origin<=truncation_upper
          Lower_Fieller <- max(Upper_origin,truncation_lower)
          Upper_Fieller <- truncation_upper
        }
        length_Fieller <- Upper_Fieller-Lower_Fieller
        if(length_Fieller==0){
          CI_type <- 'point'
          CI_Fieller <- paste('{',Lower_Fieller,'}',sep='')
        } else if(Upper_Fieller==truncation_upper & Lower_Fieller==truncation_lower){
          CI_type <- 'noninformative interval'
          CI_Fieller <- paste('[',Lower_Fieller,', ',Upper_Fieller,']',sep='')
        } else{
          CI_type <- 'continuous interval'
          CI_Fieller <- paste('[',Lower_Fieller,', ',Upper_Fieller,']',sep='')
        }
      }
    } else if(CI_type_origin=='noninformative interval'){
      Lower_Fieller <- truncation_lower
      Upper_Fieller <- truncation_upper
      length_Fieller <- Upper_Fieller-Lower_Fieller
      CI_type <- 'noninformative interval'
      CI_Fieller <- paste('[',Lower_Fieller,', ',Upper_Fieller,']',sep='')
    } else if(CI_type_origin=='point'){
      #The point is not in [truncation_lower, truncation_upper]
      if((Lower_origin-truncation_lower)*(Lower_origin-truncation_upper)>0){
        Lower_Fieller <- Upper_Fieller <- NA
        length_Fieller <- 0
        CI_type <- 'empty set'
        CI_Fieller <- 'empty set'
      }
    } else{
      #The point is in [truncation_lower, truncation_upper]
      Lower_Fieller <- Upper_Fieller <- Lower_origin
      length_Fieller <- 0
      CI_type <- 'point'
      CI_Fieller <- 'point'
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
      cat("The truncated point estimate of DC is ",dhat_Fieller, ".\n",sep='')
      cat("The truncated confidence interval of DC is ",CI_Fieller, ".\n\n",sep='')
    }

    # description <- data.frame(description=c(paste("The origin point estimate of DC is ",dhat,sep=''),
    #                                         paste("The origin confidence interval of DC is ",CI_origin,sep=''),
    #                                         paste("The truncated point estimate of DC is ",dhat_Fieller,sep=''),
    #                                         paste("The truncated confidence interval of DC is ",CI_Fieller,sep='')))
    results <- data.frame(point_estimate=c(dhat,dhat_Fieller),
                          Lower_CL=c(Lower_origin,Lower_Fieller),
                          Upper_CL=c(Upper_origin,Upper_Fieller),
                          CI=c(CI_origin,CI_Fieller),
                          length_CI=c(length_origin,length_Fieller),
                          CI_type=c(CI_type_origin,CI_type))
    rownames(results) <- c('origin','truncation')
  }

  if(H0_test==TRUE){
    Z_Fieller <- (mu_n-DC_0*mu_d)/sqrt(var_n+DC_0^2*var_d-2*DC_0*cov_nd)
    P_Fieller <- 2*pnorm(abs(Z_Fieller),0,1,lower.tail = FALSE)  # P-value for the Fieller's method
    H0_test <- data.frame(Z_score=Z_Fieller,P_value=P_Fieller)
    row.names(H0_test) <- paste('DC_0=',DC_0,sep='')
    results <- list(Estimation_result=results,
                    H0_test=H0_test)
  }

  return(results)
}
