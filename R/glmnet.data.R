#' @name glmnet.data
#' @title To Acquire The Glmnet Data for Visualization
#' @description Return a list, which contains the basic data needed for drawing, 
#'     to customize the drawing.
#'
#' @param fit The result of \link{glmnet}.
#' @param xvar Observation type, which will be x-axis mapping.
#'     You can choose "norm","lambda" or "dev".
#' @import glmnet
#' @return A list about the data of glmnet
#' @export
#'
#' @examples
glmnet.data=function(fit,xvar=c("norm","lambda","dev")){
  beta <- fit$beta
  lambda <- fit$lambda
  df <- fit$df
  dev <- fit$dev.ratio
  which=nonzeroCoef(beta)
  nwhich=length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  beta=as.matrix(beta[which,,drop=FALSE])
  xvar=match.arg(xvar)
  switch(xvar,
         "norm"={
           index=if(missing(norm))apply(abs(beta),2,sum)else norm
         },
         "lambda"={
           index=log(lambda)
         },
         "dev"= {
           index=dev
         }
  )
  dat <- cbind(t(beta),index)
  result <- list(dat =dat,index=index,df=df)
  return(result)
}

nonzeroCoef = function (beta, bystep = FALSE) { 
  nr=nrow(beta)
  if (nr == 1) {#degenerate case
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
        1
        else NULL)
    else {
      if (any(abs(beta) > 0)) 
        1
      else NULL
    }
  }
  else {
    beta=abs(beta)>0 # this is sparse
    which=seq(nr)
    ones=rep(1,ncol(beta))
    nz=as.vector((beta%*%ones)>0)
    which=which[nz]
    if (bystep) {
      if(length(which)>0){
        beta=as.matrix(beta[which,,drop=FALSE])
        nzel = function(x, which) if (any(x)) 
          which[x]
        else NULL
        which=apply(beta, 2, nzel, which)
        if(!is.list(which))which=data.frame(which)# apply can return a matrix!!
        which
      }
      else{
        dn=dimnames(beta)[[2]]
        which=vector("list",length(dn))
        names(which)=dn
        which
      }
      
    }
    else which
  }
}

