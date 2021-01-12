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


#' Add together two numbers.
#'
#' This is a description
#'
#' @param formula Similar to a formula for logistic regression.
#' @param id.set  The ID variable for matched sets.
#' @return The coefficeints, standard errors, and the covariance matrix for the estimates.
#' @export
#' @examples
#' n    <-  10
#' Y    <- rep(c(1,0,1,0),n)
#' ID   <- rep(c(1:n),each=4)
#' X1   <- rnorm(4*n)
#' X2   <- rnorm(4*n)
#' data  <- cbind(Y,X1,X2,ID)
#' clogitRV(Y~X1+X2,~ID,data)

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




