#' @name glmplot
#' @title Visualization of Lasso
#' @description Using ggplot2 to plot more elegant figure.
#' 
#' @param fit The result of \link{glmnet}.
#' @param xvar Observation type, which will be x-axis mapping.
#'     You can choose "norm","lambda" or "dev".
#' @param mode A number of Figure mode. 1 or 2.
#' @param colours Vector of colours to use for n-colour gradient.
#'     see \link{scale_color_gradientn}
#'
#' @return ggplot objest
#' @export
#' @import ggplot2
#' @import ggprism
#' @import stats
#' @import tidyr
#' @import glmnet
#' @examples
#' \dontrun{
#' library(survival)
#' library(glmnet)
#' library(glmplot)
#' data(exprset)
#' data(cliDat)
#' x <- as.matrix(exprset)
#' y <- as.matrix(Surv(cliDat[,1],cliDat[,2]))
#' fit <- glmnet(x,             
#'               y,             
#'               family = "cox",  
#'               alpha = 1,       
#'               nlambda = 1000)   
#' 
#' glmplot(fit,xvar = "dev",mode = 1,
#'    colours=c("#3B4992FF","#EE0000FF","#008B45FF","#631879FF","#008280FF","#BB0021FF"))
#' }
glmplot <- function(fit,xvar,mode,colours=NULL){
  gene=NULL;Coefficient=NULL
  #Step1.设置xlab和approx.f
  switch(xvar, #"norm","lambda","dev"
         "norm"={
           xlab <- "L1 Norm"
           approx.f <- 1},
         "lambda"={
           xlab <- "Log Lambda"
           approx.f <- 0},
         "dev"={
           xlab <- "Fraction Deviance Explained"
           approx.f <- 1
         }
  )
  #注意：
  #xvar为"norm"时，xlab为"L1 Norm"，approx的f参数为1（向下取数）
  # "lambda"时xlab为"Log Lambda"，approx的f参数为0（向上取数）
  # "dev"时，xlab为"Fraction Deviance Explained"，approx的f参数为1（向下取数）
  
  #Step2.获取绘图数据
  result <- glmnet.data(fit,xvar=xvar)
  index <- result$index;df <- result$df
  ggplot_dat <- result$dat%>%
    as.data.frame()%>%
    gather(key = "gene",value = "Coefficient",-index)
  
  #Step3.绘图
  xvalue <- max(abs(range(ggplot_dat$index)))
  if(xvalue>1){
    seqnum <- 0.5
  }else{
    if(xvalue>0.1){
      seqnum <- 0.05
    }else{
      seqnum <- 0.005
    }
  }
  xbreaks <- seq(min(ggplot_dat$index),max(ggplot_dat$index),seqnum)
  yvalue <- as.character(round(range(ggplot_dat$Coefficient),2))
  yvalue1 <- yvalue[1];yvalue2 <- yvalue[2]
  if(as.numeric(substr(yvalue1,nchar(yvalue1),nchar(yvalue1)))<=5){
    ylimit1 <- as.numeric(gsub(".$","5",yvalue1))
  }else{
    ylimit1 <- as.numeric(round(round(as.numeric(yvalue1),1),2))
  }
  if(as.numeric(substr(yvalue2,nchar(yvalue2),nchar(yvalue2)))<=5){
    ylimit2 <- as.numeric(gsub(".$","5",yvalue2))
  }else{
    ylimit2 <- as.numeric(round(round(as.numeric(yvalue2),1),2))
  }
  sec_xlabel <- approx(x=index,y=df,xout=xbreaks,rule=2,method="constant",f=approx.f)$y
  ggplot_dat$gene <- factor(ggplot_dat$gene,levels = unique(ggplot_dat$gene))
  
  if(mode==1){
    p <- ggplot(data = ggplot_dat,mapping = aes(x=index,y=Coefficient,color=gene))+
      geom_smooth(size=1.5,method = "loess",se = FALSE)+
      scale_x_continuous(breaks = xbreaks,
                         sec.axis = dup_axis(labels =  sec_xlabel),
                         guide = "prism_offset"
      )+
      scale_y_continuous(limits = c(ylimit1,ylimit2),
                         guide = "prism_offset")+
      theme_classic()+
      theme(legend.title = element_blank(),
            axis.title.x.top = element_blank())+
      scale_color_manual(values = colours)+
      xlab(xlab)
  }
  
  if(mode==2){
    a <- result$dat
    class(a[nrow(a),1:ncol(a)-1])  
    annotate_x <- rep(a[nrow(a),ncol(a)],ncol(a)-1)
    annotate_y <- a[nrow(a),1:ncol(a)-1]
    annotate_labels <- names(a[nrow(a),1:ncol(a)-1])
    xbreaks2 <- seq(min(ggplot_dat$index),max(ggplot_dat$index)+seqnum,seqnum)
    sec_xlabel2 <- approx(x=index,y=df,xout=xbreaks2,rule=2,method="constant",f=approx.f)$y
    p <- ggplot(data = ggplot_dat,mapping = aes(x=index,y=Coefficient,color=gene))+
      geom_line(size=1.5)+
      annotate("text",
               x=annotate_x,
               y=annotate_y,
               label=annotate_labels,
               hjust=0,size=3.5,fontface="bold",
               color=colours
      )+
      scale_x_continuous(limits = c(min(xbreaks2),max(xbreaks2)),
                         breaks = xbreaks2,
                         sec.axis = dup_axis(labels =  sec_xlabel2),
                         guide = "prism_offset"
      )+
      scale_y_continuous(limits = c(ylimit1,ylimit2),guide = "prism_offset")+
      theme_classic()+
      theme(legend.position = "none",axis.title.x.top = element_blank())+
      scale_color_manual(values = colours)+
      xlab(xlab)
  }
  return(p)
}

