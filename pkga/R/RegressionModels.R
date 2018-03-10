#' linmod
#'
#' Linear Regression between x and y. Returns data.frame to be displayed in UI
#' containing relationship (i.e linear), model (ie WTonCL), estimate and pval.
#' @param x Numeric vector 
#' @param y Numeric vector 
#' @param xlab Character string. X variable label.
#' @param ylab Character string. Y variable label.
#' @export

linmod <- function (x,y,xlab="X",ylab="Y"){
  cor <- lm(y~ x)
  return(c(R="Linear",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}

#' expmod
#'
#' Linear Regression between x and log(y). Returns data.frame to be displayed in UI
#' containing relationship (i.e exp), model (ie WTonCL), estimate and pval.
#' @param x Numeric vector 
#' @param y Numeric vector 
#' @param xlab Character string. X variable label.
#' @param ylab Character string. Y variable label.
#' @export

expmod <- function (x,y,xlab="X",ylab="Y"){
  if (any(is.nan(log(y)))| any(is.infinite(log(y)))) return( NULL)
  cor <- lm(log(y)~ x)
  return(c(R="Exp",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}


#' powmod
#'
#' Linear Regression between log(x) and log(y). Returns data.frame to be displayed in UI
#' containing relationship (i.e power), model (ie WTonCL), estimate and pval.
#' @param x Numeric vector 
#' @param y Numeric vector 
#' @param xlab Character string. X variable label.
#' @param ylab Character string. Y variable label.
#' @export

powmod <- function (x,y,xlab="X",ylab="Y"){
  if (any(is.nan(log(y))) | any(is.nan(log(x)))| any(is.infinite(log(x)))| any(is.infinite(log(y)))) return( NULL)
  cor <- lm(log(y)~ log(x))
  return(c(R="Power",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}

#' center
#'
#' Center by median by subtraction
#' @param x Numeric vector 
#' @export

center <- function (x){x-median(x)}

#' center2
#'
#' Center by median by division
#' @param x Numeric vector 
#' @export

center2 <- function (x){x/median(x)}