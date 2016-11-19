#' @title Unconditional Quantile Regression
#'
#'@description Function Not intended for user. Returns an object of class "urq" that represents an Unconditional Quantile Regression Fit. 
#'@usage urqb(data,tau,formula,kernel=NULL,cluster=cluster)
#'@param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.
#'@param data a data.frame in which to interpret the variables named in the formula 
#'@param tau the quantile(s) to be estimated, this must be a number (or a vector of numbers) strictly between 0 and 1.
#'@param kernel a character string giving the smoothing kernel to be used. This must match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian".
#'@param cluster column name of variable to be used in order to obtain cluster robust standard errors.
#'@import gtools Hmisc
#'@keywords NULL
#'@export
#'@seealso \code{\link{density},\link{urq}}
#'@return NULL
#'@examples NULL
urqb <- function(data=data,tau=tau,formula=formula,kernel=NULL,cluster=cluster) {
  #library(Hmisc)
  if(is.null(kernel)) kernel<-"gaussian"
  #if(is.null(tau)) tau<-1:9/10
  indicator<-function(condition) ifelse(condition,1,0)
  c1=NULL
  for (i in 1:length(tau)){
    #which(is.na(data[,idx.dep]))->miss
    #data2=data[-miss,]
    formula=as.formula(formula)
    as.data.frame(data)
    idx.dep=which(colnames(data)==all.vars(formula)[1])
    data2=data
    
    
    wtd.quantile(data2[,1],tau,weights=data$wts,normwt=TRUE)->q
    density(data[,1],kernel=kernel,weights=data$wts)->f
    approx(f$x, f$y, q)$y->fq
    RIF=q[i]+((tau[i]-indicator(data2[,idx.dep]<q[i]))/fq[i])
    
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
  class(fit) <- "urq"
  return(fit)
}