#' @title Inference for Unconditional Quantile Regression 
#'
#'@description Returns a summary list for an Unconditional Quantile Regression Fit.
#'@usage urqCI(urq,R=20,seed=NULL,colour=NULL,confidence=NULL,graph=TRUE,cluster=NULL,BC=FALSE)
#'@param urq an object of class \code{urq}.
#'@param R the number of bootstrap replications to be used.
#'@param seed random number generator.
#'@param colour colour of plot: default is lightblue.
#'@param confidence significance level.
#'@param BC plot option: If set to \code{TRUE}, Bias-Corrected Bootstrap confidence bands are plotted (black dashed lines), along with the bootstrap median (orange dashed line).
#'@param graph boolean, if \code{TRUE} a graph is produced. At least two quantiles are needed for plot to work.
#'@param cluster column name of variable to be used in order to obtain cluster robust standard errors and confidence intervals.
#'@details This function provides standard errors and confidence intervals for the Recentered Influence Function regression fit \code{urq}. If the cluster option is used, standard errors are cluster robust according to the variable supplied by the user, otherwise observations are assumed to be iid. 
#'Inference is obtained though a bayesian bootstrap drawing observation (or cluster) weights from a Dirichlet distribution. 
#'If the option graph is TRUE, then a quantile plot is provided showing estimates and confidence intervals (t approximation) or Bias-Corrected (BC) intervals. Confidence intervals using the BC percentile method typically require 1000 or more replications.
#' @references Rubin, D. B. (1981). The bayesian bootstrap. The annals of statistics, 9(1), 130-134.
#' @references Efron, B. and R. J. Tibshirani. Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. Statistical science (1986): 54-75.
#'@export 
#'@import
#'graphics grDevices
#'@keywords NULL
#'@seealso \code{\link{urq}}
#'@return NULL
#'@examples 
#'### example for cross-sectional data ###
#'
#'data(engel)
#'formula=foodexp ~ income
#'rifreg=urq(formula=formula,data=engel)
#'summary=urqCI(urq = rifreg,R = 10,graph = TRUE,seed = 1234)
#'
#'### example for panel data ###
#'
#'data(trust)
#'formula=Trust_in_the_ECB~Trust_in_the_EU+Trust_in_National_Government
#'cre=~Trust_in_the_EU+Trust_in_National_Government
#'rif=urq(formula,data=trust,cre=cre,id="countryname")
#'summary=urqCI(urq = rif,R = 10,graph = TRUE,seed = 1234,cluster="countryname")

urqCI=function(urq,R=20,seed=NULL,colour=NULL,confidence=NULL,graph=TRUE,cluster=NULL,BC=FALSE){
  if (is.null(seed)) {
    seed <- runif(1, 0, .Machine$integer.max)
  }
  
  
  set.seed(seed)
  boot=matrix(ncol=dim(urq$coefficients)[1]*length(urq$tau),nrow=R)
  
  #data=urq$data
  cat("Bootstrapping\n")
  #print(R)
  for(i in 1:R){
    #bayesian bootstrap
    if (is.null(cluster)){urq$data$wts=as.numeric(rdirichlet(1, rep(1, dim(urq$data)[1])))}
    #frequentist bootstrap
    # if (is.null(cluster)){
    #   freq=as.numeric(rdirichlet(1, rep(1, dim(urq$data)[1])))
    #   urq$data$wts=freq=freq/sum(freq)}
    else{
      idx.cluster=which(colnames(urq$data)==cluster)
      clusters=names(table(urq$data[,idx.cluster]))
      urq$data=urq$data[order(urq$data[,idx.cluster]),] 
      freq=as.numeric(rdirichlet(1, rep(1, length(clusters))))
      freq=rep(freq,times=table(urq$data[,idx.cluster]))
      urq$data$wts=freq/sum(freq)
    }
    #data=data,tau=tau,formula=formula,kernel=NULL,cluster=cluster
    boot[i,]=as.numeric(urqb(data=urq$data,tau=urq$tau,formula=urq$formula,kernel=NULL,cluster=cluster)$coefficients)
    #print(i)
    #weights=NULL
    if (!i %% 50)cat(".")
  }
  
  if (R<=50) cat("Number of replications is small, please consider increasing it.")
  names.urq=rownames(urq$coefficients)
  if(is.null(confidence)) confidence<-0.95
  alpha=1-confidence
  se=apply(boot,2,sd,na.rm=TRUE)
  
  
  #colnames(boot)=paste(rep(rownames(urq),times=dim(urq)[2]),rep(tau,each=dim(urq)[1]))
  colnames(boot)=paste(rep(rownames(urq$coefficients),times=dim(urq$coefficients)[2]))
  if(is.null(colour)) colour<-"lightblue"
  mar <- c(2, 2, 2, 1.6) #margins for plot
  
  par(mfrow=n2mfrow(dim(urq$coefficients)[1]-1),mar=mar)
  #par(mfrow=c(rep(floor(sqrt(dim(urq$coefficients)[1])),2)),mar=mar)
  par(oma=c(0,0,3,0))
  
  se=apply(boot,2,sd,na.rm=TRUE)
  urq.vector=as.numeric(urq$coefficients)
  tP <- 2 * pt(-abs(urq.vector/se), dim(boot)[1] - 1)
  lower <- urq.vector + qt(alpha/2, dim(boot)[1] - 1) * se
  upper <- urq.vector + qt(1 - alpha/2, dim(boot)[1] - 1) * se
  recap=cbind(urq.vector,lower,upper,tP,se,rep(urq$tau,each=dim(urq$coefficients)[1]))
  colnames(recap)=c("Coef","t Lower","t Upper","P>|t|","Std.Err.","tau")
  
  ######BCA
  sapply(1:dim(boot)[2], function(x) qnorm(mean(boot[,x] <= urq.vector[x],na.rm=TRUE)))->z0
  #sapply(1:dim(jack)[2], function(x)sum((mean(jack[,x])-jack[,x])^3)/(6*sum((mean(jack[,x])-jack[,x])^2)^(3/2)))->a
  a=rep(0,length(z0))
  
  q.lb=sapply(1:dim(boot)[2], function(x) pnorm(z0[x]+(z0[x]+qnorm(alpha/2))/(1-a[x]*(z0[x]+qnorm(alpha/2)))))
  q.ub=sapply(1:dim(boot)[2], function(x) pnorm(z0[x]+(z0[x]+qnorm(1-alpha/2))/(1-a[x]*(z0[x]+qnorm(1-alpha/2)))))
  type=7
  sapply(1:dim(boot)[2], function(x) quantile(boot[,x],q.lb[x],na.rm=TRUE,type=type))->BClower
  sapply(1:dim(boot)[2], function(x) quantile(boot[,x],q.ub[x],na.rm=TRUE,type=type))->BCupper
  
  recap=cbind(recap,BClower,BCupper,apply(boot,2,median,na.rm=TRUE))
  
  cairo_pdf(height = 7, width = 10,onefile=T)
  if (isTRUE(graph) && isTRUE(length(urq$tau)>1)){
    for(i in c(2:dim(urq$coefficients)[1])){
      cond=names.urq[i]
      recap[which(rownames(recap) == cond),]->cfi
      #recap[grep(cond, rownames(recap),fixed=TRUE), ]->cfi  
      #cfi=cbind(cfi[,1],cfi[,2],cfi[,3])
      
      #plot(rep(tau, 2), c(cfi[, 2], cfi[, 3]),xaxt = "n" ,type = "n",main=cond,xlab="",ylab="",ylim=c(min=min(coef(rq)[cond,],coef(canay)[cond,],coef(canay.fd)[cond,],cfi[,2],cfi[,8],urq[cond,]),max=max(coef(rq)[cond,],cfi[,3],cfi[,9],urq[cond,],coef(canay)[cond,],coef(canay.fd)[cond,])))
      plot(rep(urq$tau, 2), c(cfi[, 2], cfi[, 3]),xaxt = "n" ,type = "n",main=cond,xlab="",ylab="",ylim=c(min=min(cfi,na.rm=TRUE),max=max(cfi,na.rm=TRUE)))
      axis(1, at=urq$tau, labels=urq$tau,cex.axis=.8,las=2)
      
      polygon(c(urq$tau, rev(urq$tau)),c(cfi[,2],rev(cfi[,3])),col=colour,border=NA)
      points(urq$tau,cfi[,1],type="b",lty="longdash",cex = 0.5,pch = 20,col="blue")
      

      abline(h=0)
      if (isTRUE(BC)){
        points(urq$tau,cfi[,7],type="b",lty="longdash",cex = 0.5,pch = 20,col="black")
        points(urq$tau,cfi[,8],type="b",lty="longdash",cex = 0.5,pch = 20,col="black")
        #plots median of bootstrap distribution
        points(urq$tau,cfi[,9],type="b",lty="longdash",cex = 0.5,pch = 20,col="orange")
      }
    }
    dev.off()
  }
  
  
  fit=list (results=recap,bootstrap=boot)
  class(fit) <- "urqCI"
  return(fit)
  
  
}