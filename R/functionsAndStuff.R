####################################################################################
###
###GET THE ROBUST VERSION of the INFORMATION MATRIX
###
####################################################################################

getI.robust2  <- function(beta,data3,Yname,Xnames,IDname){
  betac2         <- matrix(beta,ncol=1)
  XS             <- c(Xnames,"strata")
  strata         <- unique(data3[,"strata"])
  ID             <- unlist(strsplit(strata,"\\."))[c(1:length(strata))%%2==1]

  XA             <- data3[data3[,Yname]==1,XS]
  XA             <- as.matrix(XA[match(strata,XA[,"strata"]),Xnames])

  XU             <- data3[data3[,Yname]==0,XS]
  XU             <- as.matrix(XU[match(strata,XU[,"strata"]),Xnames])

  XDIF           <- as.matrix(XA-XU)
  XDIF           <- as.matrix(XA-XU)
  mu             <- exp(XA%*%betac2)/(exp(XA%*%betac2)+exp(XU%*%betac2))
  mu.mat         <- matrix(mu,nrow=nrow(XA),ncol=ncol(XA),byrow=FALSE)
  Xbar           <- mu.mat*XA+(1-mu.mat)*XU
  EE2            <- (XA-Xbar)
  EE             <- EE2 - matrix(colMeans(EE2),nrow=nrow(EE2),ncol=ncol(EE2),byrow=TRUE)
  p              <- ncol(XA)
  meat           <- matrix(0,nrow=p,ncol=p)
  uID            <- unique(ID)
  n              <- length(uID)
  colSums2       <- function(mat) {out = mat; if(nrow(mat)>1) out <- colSums(mat); return(out)}
  I              <- matrix(nrow=p,ncol=p)
  for (i in 1:p) for (j in 1:p) I[i,j] <- sum(XDIF[,i]*XDIF[,j]*mu*(1-mu))
    for (i in 1:n){
    IDw <- c(1:nrow(XDIF))[ID==uID[i]]
    if (p == 1)                   EEi   <- matrix(sum(EE[IDw,]),nrow=1,ncol=1)
    if (length(IDw) >  1 & p > 1) EEi   <- matrix(colSums2(mat=EE[IDw,]),ncol=1)
    if (length(IDw) == 1 & p > 1) EEi   <- matrix(EE[IDw,],ncol=1)
    meat <- meat + EEi%*%t(EEi)
  }
  bread  <- I
  breadi <- solve(bread)
  solve(t(breadi)%*%meat%*%breadi)
}


####################################################################################
###
###GET THE Degees of Freedom (ONLY USED WHEN CHOOSING TUNING PARAMETER)
###
####################################################################################
fake.function <- function(){
n    <-  10
Y    <- rep(c(1,0,1,0),n)
ID   <- rep(c(1:n),each=4)
X1   <- rnorm(4*n)
X2   <- rnorm(4*n)
Yname <- "Y"
Xnames <- c("X1","X2")
IDname   <- "ID"
formula2 <- "Y ~ X1 + X2"
data  <- cbind(Y,X1,X2,ID)
}


#' Conditional Logistic Regression with Robust Variance
#'
#' Performs conditional logistic regression with a robust variance estimator
#'
#' @param formula Similar to a formula for glm.
#' @param id.set  The ID variable for matched sets.
#' @param data    The dataset
#' @return coef   The coefficeints (i.e. betas) and their standard error
#' @return cov    The covariance matrix of the coefficeints (useful for wald tests)
#' @export
#' @examples
#' ###A simple Example
#' n    <-  10
#' Y    <- rep(c(1,0,1,0),n)
#' ID   <- rep(c(1:n),each=4)
#' X1   <- rnorm(4*n)
#' X2   <- rnorm(4*n)
#' data  <- cbind(Y,X1,X2,ID)
#' clogitRV(Y~X1+X2,~ID,data)
#' ##
#' ##
#' ##
#' ##An example with a time-varying coefficient
#' ##simData is an example dataset included in the package
#' ?simData
#'
#' ###Run the Model (note: z is just a covariate)
#' simData$x1 <- simData$t   * simData$x
#' simData$x2 <- simData$t^2 * simData$x
#' model      <- clogitRV(Y~x+x1+x2+z,~ID, data=simData)
#' beta       <- model$coef[-4,1]
#' cov        <- model$cov[-4,-4]
#'
#' ##Draw the picture
#' library(ggplot2)
#' t2       <- seq(min(simData$t),max(simData$t),length.out=100)
#' t2M      <- cbind(1,t2,t2^2)
#' beta2     <- t2M %*% beta
#' se2       <- sqrt(diag(t2M%*%cov%*%t(t2M)))
#' plotData  <- as.data.frame(t2,beta2,se2)
#' ggplot(data=plotData,aes(t2,beta2)) +
#'   geom_line() +
#'   geom_point()+
#'   geom_errorbar(aes(ymin=beta2-1.96*se2, ymax=beta2+1.96*se2))+
#'   xlab("Time Until Diaganosis")+
#'   ylab("log(RR)")
#'
#' ##Get the p-values
#' ##P-value for overall effect of the biomarker
#' 1-pchisq(t(beta)%*%solve(cov)%*%beta,length(beta))
#' ##P-value for non-linear effect of the biomarker
#' 1-pchisq(t(beta[-1])%*%solve(cov[-1,-1])%*%beta[-1],length(beta[-1]))
#'
clogitRV=function(formula,id.set,data){

  require(tidyverse)
  require(survival)

  ### Get the variable names
  f1     <- strsplit(as.character(formula),"~")
  f2     <- strsplit(as.character(id.set),"~")
  Yname  <- f1[[2]][1]
  Xnames <- str_trim(unlist(strsplit(f1[[3]][[1]],"\\+")))
  IDname <- f2[[2]][1]


  ### Expand the data
  data    <- as.data.frame(data)
  datab   <- data[complete.cases(data[,c(Yname,Xnames)]),]
  n       <- nrow(datab)
  Y       <- datab[,Yname]
  setID   <- datab[,IDname]
  usetID  <- unique(setID)
  nsets   <- length(usetID)

  newData  <- 0
  for (set in 1:nsets){
    caseRows <- c(1:n)[setID==usetID[set] & Y == 1]
    ncase    <- length(caseRows)
    contRows <- c(1:n)[setID==usetID[set] & Y == 0]
    ncont    <- length(contRows)
    if (ncase > 0 & ncont >0){
        newData  <- newData + 1
        rowsWant <- c(rep(caseRows,each=ncont),rep(contRows,ncase))
        strata   <- paste0(usetID[set],".",rep(c(1:(ncase*ncont)),2))
        data2    <- cbind(datab[rowsWant,],strata)
        if (newData == 1) data3 <- data2
        if (newData >  1) data3 <- rbind(data3,data2)
    }
  }

  ###Run clogit on the expanded dataset
  formula2 <- paste0(formula[2],formula[1],formula[3])
  command  <- paste0("clogit(",formula2,"+strata(strata),data=data3)")
  tm       <- eval(parse(text=command))
  beta     <- tm$coef


  ####Get the Robust standard error
  V        <- solve(getI.robust2(beta,data3,Yname,Xnames,IDname))
  se       <- sqrt(diag(V))
  coef     <- cbind(beta,se)
  cov      <- V
  return(list(coef=coef,cov=cov))
}

############################################################################################
#
#' A Simulated Dataset
#'
#' A simulated dataset for showing a time-varying effect
#'
#' @format A data frame with 100 matched case/control pairs with two measurements per individual
#' \describe{
#' \item{Y}{Case/Control Status}
#' \item{x}{biomarker}
#' \item{t}{time between biomarker measurement and case diagnosis}
#' \item{z}{a covariate}
#' }
"simData"

