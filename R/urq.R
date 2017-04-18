#' @title Unconditional Quantile Regression
#'
#'@description Returns an object of class \code{urq}. that represents an Unconditional Quantile Regression Fit
#'@usage urq(formula,data,tau=NULL,kernel=NULL,cre=NULL,id=NULL)
#'@param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.
#'@param data a \code{dataframe} in which to interpret the variables named in the formula 
#'@param tau the quantile(s) to be estimated, this must be a number (or a vector of numbers) strictly between 0 and 1.
#'@param kernel a character string giving the smoothing kernel to be used. This must match one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine", with default "gaussian".
#'@param cre The CRE formula (right hand side only) is a specification of the variables in the CRE component. These are possibly endogenous variables (in the sense that they are affected by the fixed effects) and must be time-varying. If left empty, a cross-sectional analysis is performed.
#'@param id defines the structure of the panel.
#'@details This function returns a Recentered Influence Function regression of given quantiles as proposed by Firpo, S., Fortin, N. M., & Lemieux, T. (2009). Panel data analysis is performed extending the correlated random effects (CRE)  model by Mundlak (1978) and Chamberlain (1984) to an unconditional quantile regression framework. See Abrevaya and Dahl (2008) and Bache et al (2011) for more details.
#'@export 
#'@import
#'stats
#'@references Firpo, S., Fortin, N. M., & Lemieux, T. (2009). Unconditional quantile regressions. Econometrica, 77(3), 953-973.
#'@references Mundlak, Y. 1978. On the pooling of time series and cross section data. Econometrica 46: 69-85.
#'@references Chamberlain G (1984) Panel Data. In: Griliches Z, Intriligator MD (eds) Handbook of Econometrics, vol 2, Elsevier Science B. V., pp 1247-1318
#'@references Abrevaya, Jason and Christian M. Dahl. 2008. The effects of birth inputs on birthweight. Jounal of Business and Economic Statistics. 26-4. Pages 379-397.
#'@references Bache, Stefan Holst; Christian M. Dahl; Johannes Tang Kristensen. 2011. Headlights on tobacco road to low birthweight - Evidence from a battery of quantile regression estimators and a heterogeneous panel.
#'@keywords NULL
#'@seealso \code{\link{density},\link[uqr]{urqCI}}
#'@return NULL
#'@examples 
#' ### example for cross-sectional data ###
#' 
#'data(engel)
#'formula = foodexp ~ income
#' rifreg=urq(formula,data = engel)
#' 
#' ### example for panel data ###
#' 
#' data(trust)
#'formula=Trust_in_the_ECB~Trust_in_the_EU+Trust_in_National_Government
#'cre=~Trust_in_the_EU+Trust_in_National_Government
#'rif=urq(formula,data=trust,cre=cre,id="countryname")
urq <- function(formula,data,tau=NULL,kernel=NULL,cre=NULL,id=NULL) {
  #call <- match.call()
  #mf <- match.call(expand.dots = FALSE)
  
  if(is.null(kernel)) kernel<-"gaussian"
  if(is.null(tau)) tau<-1:9/10
  indicator<-function(condition) ifelse(condition,1,0)
  c1=NULL
  idx.dep=which(colnames(data)==all.vars(formula)[1])
  quantile(x=data[,idx.dep],probs=tau)->q
  density(data[,idx.dep],kernel=kernel)->f
  approx(f$x, f$y, q)$y->fq
  for (i in 1:length(tau)){
    #which(is.na(data[,idx.dep]))->miss
    #data2=data[-miss,]


    
    

    RIF=q[i]+((tau[i]-indicator(data[,idx.dep]<q[i]))/fq[i])
    data2=data
    data2[,idx.dep]=RIF
    
    if (!is.null(cre)){
      if(is.null(id)){print("id variable is missing")}
      vars=all.vars(cre)
      new=rep(NA,length(vars))
      old=all.vars(formula)[2:length(all.vars(formula))]
      dep=all.vars(formula)[1]
      for(i in 1:length(vars)){

        data2[,dim(data2)[2]+1] <- ave(data[,which(colnames(data2)==vars[i])], group=data2[,which(colnames(data2)==id)])
        new[i]=colnames(data2)[dim(data2)[2]]<-paste(vars[i],sep="","_M")
      }
      formula_cre=paste(dep, "~", paste(c(old,new), collapse=" + "))
      lm=lm(formula_cre,data2,weights=NULL)
      c=coef(lm)
      c1=rbind(c1,c)
      #print(i)
      data2[,idx.dep]=data[,idx.dep]
    }
    
  else{
    lm=lm(formula,data2,weights=NULL)
  c=coef(lm)
  c1=rbind(c1,c)
#print(i)
  data2=NULL}
    
    
  }
  RIF=t(c1)
  colnames(RIF)<-paste("tau=",tau)
  if (!is.null(cre)){formula=formula_cre}
  if (!is.null(cre)){
  data=data2
    #data=data.frame(data,data2[,dim(data)[2]+1])
  #colnames(data[,idx.dep])=all.vars(formula)[1]
  }
  
  fit=list (coefficients=RIF,tau=tau,formula=formula,id=id,data=data)
  fit$call <- sys.call()
  class(fit) <- "urq"
  #return(call)
  return(fit)
}