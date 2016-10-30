#' @title Unconditional Quantile Regression
#'
#'@description Returns a summary list for an Unconditional Quantile Regression Fit.
#'@usage urqCI(urq,data,R=20,seed=NULL,colour=NULL,confidence=NULL,graph=TRUE,cluster=NULL)
#'@param urq an object of class \code{urq}.
#'@param data a \code{dataframe} in which to interpret the variables named in the \code{formula} .
#'@param R the number of bootstrap replications to be used.
#'@param seed random number generator.
#'@param colour colour of plot: default is lightblue.
#'@param confidence significance level.
#'@param graph boolean, if \code{TRUE} a graph is displayed. At least two quantiles are needed for plot to work.
#'@param cluster column name of variable to be used in order to obtain cluster robust standard errors and confidence intervals.
#'@details This function provides standard errors and confidence intervals for the Recentered Influence Function regression fit \code{urq}. If the cluster option is used, standard errors are cluster robust according to the variable supplied by the user, otherwise observations are assumed to be iid. 
#'Inference is obtained though a bayesian bootstrap drawing observation (or cluster) weights from a Dirichlet distribution. 
#'If the option graph is TRUE, then a quantile plot is provided showing estimates and confidence intervals (t approximation [polygon] and percentile bootstrap [dashed lines]).
#' @references Rubin, D. B. (1981). The bayesian bootstrap. The annals of statistics, 9(1), 130-134.
#'@export 
#'@import
#'graphics grDevices
#'@keywords NULL
#'@seealso \code{\link{urq}}
#'@return NULL
#'@examples data(engel)
#'formula=foodexp ~ income
#'rifreg=urq(formula=formula,data=engel)
#'summary=urqCI(urq = rifreg,data = engel,R = 100,graph = TRUE,seed = 1234)

urqCI=function(urq,data,R=20,seed=NULL,colour=NULL,confidence=NULL,graph=TRUE,cluster=NULL){
  if (is.null(seed)) {
    seed <- runif(1, 0, .Machine$integer.max)
  }
  
  
  set.seed(seed)
  boot=matrix(ncol=dim(urq$coefficients)[1]*length(urq$tau),nrow=R)
  
  cat("Bootstrapping\n")
  #print(R)
  for(i in 1:R){
    if (is.null(cluster)){data$wts=as.numeric(rdirichlet(1, rep(1, dim(data)[1])))}
    else{
      idx.cluster=which(colnames(data)==cluster)
      clusters=names(table(data[,idx.cluster]))
      data=data[order(data[,idx.cluster]),] 
      freq=as.numeric(rdirichlet(1, rep(1, length(clusters))))
      freq=rep(freq,times=table(data[,idx.cluster]))
      data$wts=freq/sum(freq)
      
    }
    
    boot[i,]=as.numeric(urqb(urq$formula,data=data,tau=urq$tau,cluster=cluster)$coefficients)
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
  recap=cbind(urq.vector,t(apply(boot,2,quantile,c(alpha/2,1-alpha/2),na.rm=TRUE)),lower,upper,tP,se,rep(urq$tau,each=dim(urq$coefficients)[1]))
  colnames(recap)=c("Coef","Perc. Lower","Perc. Upper","t Lower","t Upper","P>|t|","Std.Err.","tau")
  
  if (isTRUE(graph) && isTRUE(length(urq$tau)>1)){
    for(i in c(2:dim(urq$coefficients)[1])){
      cond=names.urq[i]
      recap[which(rownames(recap) == cond),]->cfi
      #recap[grep(cond, rownames(recap),fixed=TRUE), ]->cfi  
      #cfi=cbind(cfi[,1],cfi[,2],cfi[,3])
      
      #plot(rep(tau, 2), c(cfi[, 2], cfi[, 3]),xaxt = "n" ,type = "n",main=cond,xlab="",ylab="",ylim=c(min=min(coef(rq)[cond,],coef(canay)[cond,],coef(canay.fd)[cond,],cfi[,2],cfi[,8],urq[cond,]),max=max(coef(rq)[cond,],cfi[,3],cfi[,9],urq[cond,],coef(canay)[cond,],coef(canay.fd)[cond,])))
      plot(rep(urq$tau, 2), c(cfi[, 2], cfi[, 3]),xaxt = "n" ,type = "n",main=cond,xlab="",ylab="",ylim=c(min=min(cfi[,1:5],na.rm=TRUE),max=max(cfi[,1:5],na.rm=TRUE)))
      axis(1, at=urq$tau, labels=urq$tau,cex.axis=.8,las=2)
      
      polygon(c(urq$tau, rev(urq$tau)),c(cfi[,4],rev(cfi[,5])),col=colour,border=NA)
      points(urq$tau,cfi[,1],type="b",lty="longdash",cex = 0.5,pch = 20,col="blue")
      points(urq$tau,cfi[,2],type="l",lty="longdash",cex = 0.5,pch = 20,col="black")
      points(urq$tau,cfi[,3],type="l",lty="longdash",cex = 0.5,pch = 20,col="black")
      abline(h=0)
    }
  }
  class(recap)="urqCI"
  return(recap)
  
}