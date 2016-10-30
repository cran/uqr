#' @title Unconditional Quantile Regression
#'
#'@description Returns an object of class \code{urq}. that represents an Unconditional Quantile Regression Fit
#'@usage urq(formula,data,tau=NULL,kernel=NULL)
#'@param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.
#'@param data a \code{dataframe} in which to interpret the variables named in the formula 
#'@param tau the quantile(s) to be estimated, this must be a number (or a vector of numbers) strictly between 0 and 1.
#'@param kernel a character string giving the smoothing kernel to be used. This must match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian".
#'@details This function returns a Recentered Influence Function regression of given quantiles as proposed by Firpo, S., Fortin, N. M., & Lemieux, T. (2009).
#'@export 
#'@import
#'stats
#'@references Firpo, S., Fortin, N. M., & Lemieux, T. (2009). Unconditional quantile regressions. Econometrica, 77(3), 953-973.
#'@keywords NULL
#'@seealso \code{\link{density},\link[uqr]{urqCI}}
#'@return NULL
#'@examples 
#'data(engel)
#'formula = foodexp ~ income
#' rifreg=urq(formula,data = engel)
urq <- function(formula,data,tau=NULL,kernel=NULL) {
  #call <- match.call()
  #mf <- match.call(expand.dots = FALSE)
  
  if(is.null(kernel)) kernel<-"gaussian"
  if(is.null(tau)) tau<-1:9/10
  indicator<-function(condition) ifelse(condition,1,0)
  c1=NULL
  for (i in 1:length(tau)){
    #which(is.na(data[,idx.dep]))->miss
    #data2=data[-miss,]
    idx.dep=which(colnames(data)==all.vars(formula)[1])
    
    
    quantile(x=data[,idx.dep],probs=tau)->q
    density(data[,idx.dep],kernel=kernel)->f
    approx(f$x, f$y, q)$y->fq
    RIF=q[i]+((tau[i]-indicator(data[,idx.dep]<q[i]))/fq[i])
    data2=data
    data2[,idx.dep]=RIF
    lm=lm(formula,data2,weights=NULL)
    c=coef(lm)
    c1=rbind(c1,c)
    #print(i)
    data2=NULL
  }
  RIF=t(c1)
  colnames(RIF)<-paste("tau=",tau)
  
  fit=list (coefficients=RIF,tau=tau,formula=formula)
  fit$call <- sys.call()
  class(fit) <- "urq"
  #return(call)
  return(fit)
}