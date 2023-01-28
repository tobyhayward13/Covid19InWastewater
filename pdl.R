
lagged<-function(x,lag=1){
  if (lag==0) return(x)
  n<-length(x)
  c(rep(NA,lag),x[-( (n-lag+1):n)])
}

pdlweights<-function(lag,degree,tiedown=c(F,F)){
  if (any(tiedown)) stop("Tiedown not working")
  contr.poly(lag,contrasts = FALSE)[,tiedown[1]+(1:degree)]
}

# pdl<-function(x,lag,degree,tiedown=c(F,F)){
#   xlags<-as.matrix(sapply((1:lag)-1,function(i) lagged(x,i)))
#   pnom<-pdlweights(lag,degree,tiedown)
#   z<-xlags%*%pnom
#   attr(z,"call")<-match.call()
#   attr(z,"weights")<-pnom
#   z
# }

pdl<-function(xlags,lag,degree,tiedown=c(F,F)){
  m<-NCOL(xlags)
  if (lag>m)
    stop("more lags than variables")
  else
    xlags<-xlags[,1:lag,drop=FALSE]
  pnom<-pdlweights(lag,degree,tiedown)
  z<-xlags%*%pnom
  attr(z,"call")<-match.call()
  attr(z,"weights")<-pnom
  z
}


##
## redoes attr(modelmatrix,"assign") in the nice S-PLUS 3.4 format
##
attrassign<-function (object, ...) UseMethod("attrassign")

attrassign.lm<-function(lmobj){
  attrassign(model.matrix(lmobj),terms(lmobj))}

attrassign.default<-function(mmat,tt){
  if (!inherits(tt,"terms"))
    stop("need terms object")
  aa<-attr(mmat,"assign")
  if (is.null(aa))
    stop("argument is not really a model matrix")
  ll<-attr(tt,"term.labels")
  if (attr(tt,"intercept")>0)
    ll<-c("(Intercept)",ll)
  aaa<-factor(aa,labels=ll)
  split(order(aa),aaa)
}

##
## Fits polynomial distributed lag models
##

pdlglm<-function(formula,data,family,action="na.omit",...){
  fformula<-formula
  rm(formula)
  tt<-terms(fformula,specials="pdl")
  sp<-untangle.specials(tt,"pdl",2:10)	
  if (length(sp$terms)!=0) 
    stop("Can't handle interactions with pdl() terms")
  glmfit<-match.call()
  glmfit[[1]]<-as.name("glm")
  glmfit<-eval(glmfit,sys.frame(sys.parent()))
  sp<-untangle.specials(tt,"pdl")
  ass<-attrassign(glmfit)
  beta<-coef(glmfit)
  p<-length(beta)
  slop<-0
  for (term in sp$vars){
    wterm<-parse(text=term)[[1]]
    wterm[[1]]<-as.name("pdlweights")
    if (!is.null(wterm$x)) {
      xn<-wterm$x
      wterm$x<-NULL
    }else{
      xn<-wterm[[2]]
      wterm[[2]]<-NULL
    }
    wts<-eval(wterm,sys.frame(sys.parent()))
    which<-ass[[term]]+slop
    before<-(1:p)[(1:p)<min(which)]
    after<-(1:p)[(1:p)>max(which)]
    nb<-names(beta)
    wnames<-paste("lag(",xn,",",(1:nrow(wts))-1,")",sep="")
    beta<-c(beta[before],beta[which]%*%t(wts),beta[after])
    names(beta)<-c(nb[before],wnames,nb[after])
    slop<-slop+nrow(wts)-length(which)
    p<-p+nrow(wts)-length(which)
  }
  glmfit$lagcoef<-beta
  glmfit$ass<-ass
  glmfit$sp<-sp$vars 
  class(glmfit)<-c("pdlglm",class(glmfit))
  glmfit
}



plot.pdlglm<-function(pdlobj,xlab="lag",ylab="Coefficient",...){
  for (term in pdlobj$sp){
    pp<-pdlobj$ass[[term]]
    betas<-pdlobj$coef[pp]
    wterm<-parse(text=term)[[1]]
    wterm[[1]]<-as.name("pdlweights")
    if (!is.null(wterm$x)) {
      xn<-wterm$x
      wterm$x<-NULL
    }else{
      xn<-wterm[[2]]
      wterm[[2]]<-NULL
    }
    wts<-eval(wterm,sys.frame(sys.parent()))
    betas<-betas%*%t(wts)
    if (is.null(ylim)) 
      yylim<-range(c(0,betas))
    else 
      yylim<-ylim
    plot(0:((length(betas)-1)),betas,type="h",xlab=xlab,ylab=ylab,main=term,ylim=yylim,...)
    abline(h=0,lty=2)
  }       
  invisible(NULL)
}

## Martin Maechler's print.glm function from R
print.pdlglm<- function (x, digits = max(3, .Options$digits - 3), na.print = "", 
                         ...) 
{
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Coefficients:\n")
  print.default(round(x$lagcoef, digits), print.gap = 2)
  cat("\nDegrees of Freedom:", x$df.null, "Total; ", x$df.residual, 
      "Residual\n")
  cat("Null Deviance:", format(signif(x$null.deviance, 
                                      digits)), "\n")
  cat("Residual Deviance:", format(signif(x$deviance, digits)), 
      "\n")
  invisible(x)
}

summary.pdlglm<-function (object, dispersion = NULL, correlation = FALSE, na.action = na.omit) 
{	
  ans<-NextMethod()
  covmat<-ans$cov.unscaled
  slop<-0
  for (term in object$sp){
    wterm<-parse(text=term)[[1]]
    wterm[[1]]<-as.name("pdlweights")
    if (!is.null(wterm$x)) {
      xn<-wterm$x
      wterm$x<-NULL
    }else{
      xn<-wterm[[2]]
      wterm[[2]]<-NULL
    }
    wts<-eval(wterm,sys.frame(sys.parent()))
    which<-object$ass[[term]]+slop
    p<-ncol(covmat)
    before<-(1:p)[(1:p)<min(which)]
    after<-(1:p)[(1:p)>max(which)]
    c11<-covmat[before,before,drop=F]
    c12<-covmat[before,which,drop=F]
    c13<-covmat[before,after,drop=F]
    c22<-covmat[which,which,drop=F]
    c23<-covmat[which,after,drop=F]
    c33<-covmat[after,after,drop=F]
    covmat<-wts%*%c22%*%t(wts)
    if (length(before)!=0) {
      d12<-c12%*%t(wts)
      covmat<-rbind(cbind(c11,d12),
                    cbind(t(d12),covmat))
    }	
    if (length(after)!=0){
      d23<-wts%*%c23
      if (length(before)!=0){
        covmat<-rbind(cbind(covmat, rbind(c13,d23)),
                      cbind(t(c13),t(d23), c33))
      }else{	
        covmat<-rbind(cbind(covmat, d23),
                      cbind(t(d23),c33))
      }
    }
    slop<-slop+nrow(wts)-length(which)
    
  }
  ans$cov.unscaled<-covmat
  if(correlation)
    ans$correlation<-covmat/sqrt(outer(diag(covmat),diag(covmat)))
  else
    ans$correlation<-NULL
  ans$coefficients<-cbind(object$lagcoef,sqrt(diag(covmat)),object$lagcoef/sqrt(diag(covmat)))
  dimnames(ans$coefficients)[[2]]<-c("Value","Std. Error","t value")
  ans
}



## Terry Therneau's function from survival4
untangle.specials<-function (tt, special, order = 1) 
{
  #
  # There was a change in the specials, so action depends on your relea
  #   of S
  #
  spc <- attr(tt, "specials")[[special]]
  ###<TSL>    facs <- attr(tt, 'factor')
  if (length(spc) == 0) 
    return(vars = character(0), terms = numeric(0))
  facs <- attr(tt, "factors")
  ###</TSL>
  fname <- dimnames(facs)
  if ((attr(terms(y ~ zed(x), specials = "zed"), "specials"))$zed == 
      1) {
    # old style
    if (any(order > 1)) 
      warning("Can't select specials based on the order of  terms")
    list(vars = (fname[[2]])[spc], terms = spc)
  }
  else {
    ff <- apply(facs[spc, , drop = F], 2, sum)
    list(vars = (fname[[1]])[spc], terms = seq(ff)[ff & 
                                                     match(attr(tt, "order"), order, nomatch = 0)])
  }
}
