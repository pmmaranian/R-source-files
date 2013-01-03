library(foreign)
library(gmodels)
library(xtable)
library(ggplot2)
library(XLConnect)
library(reshape)

NmeanSD <- function(x) {
  x<-x[!is.na(x)]
  N<-as.character(length(x))
  mean<-as.character(round(mean(x),2))
  SD<-as.character(round(sd(x),2))
  return(c(N,paste(mean," (",SD,")",sep="")))
}

tblcontvars <- function(y,ds,d) {
  x <- ds[[pmatch(y,names(ds))]]
  D <- ds[[pmatch(d,names(ds))]]
  A <- unlist(with(ds,unlist(tapply(x,D,NmeanSD))))
  if (length(levels(D))>2) p <- round(if (shapiro.test(x)[[2]] < .05) with(ds,kruskal.test(x~D))$p.value else summary(with(ds,aov(x~D)))[[1]][1,5],4)
  else p <- round(if (shapiro.test(x)[[2]] < .05) with(ds,wilcox.test(x~D))$p.value else with(ds,wilcox.test(x~D))$p.value,4)
  if (p < .0001) p <- "<.0001"
  return(c(y,A,p))
}
#tblcontvars("ankle.Sw",wc,"BMIcat2")

tblcatvars <- function(y,ds,d) {
  x <- ds[[pmatch(y,names(ds))]]
  D <- factor(ds[[pmatch(d,names(ds))]])
  Ns <- with(ds,table(x,D))
  props <- round(100*prop.table(Ns,2),1)
  v <- paste(as.vector(t(Ns))," (",as.vector(t(props)),"\\%)",sep="")
  m <- c()
  for (i in 1:length(v)) 
    m <- c(m,"",v[i])
  top <- c()
  tots <- apply(Ns,2,sum)
  for (i in 1:with(ds,length(levels(D))))
    top <- c(top,tots[i],"")
  ft <- try(fisher.test(Ns),silent=T)
  if (class(ft)=="try-error" | tots[1]==0 | tots[2]==0) top <- c(top,"NA")
  else top <- c(top,as.character(round(ft$p.value,4)))
  m1 <- cbind(matrix(m,nrow=length(levels(x)),byrow=T),rep("",length(levels(x))))
  m1 <- rbind(top,m1)
  m1 <- cbind(c(paste("\\textbf{",y,"}",sep=""),paste("\\hspace*{5mm}",levels(x),sep="")),m1)
  rownames(m1) <- NULL
  return(m1)
}
#tblcatvars("totulcer",rbind(CRISS,three.RCTs.spl[[1]]),"ds")

maketable <- function(Y,DS,d) {
  X <- DS[[pmatch(Y,names(DS))]]
  D <- DS[[pmatch(d,names(DS))]]
  if (is.factor(X)) {
    try1 <- try(tblcatvars(Y,DS,d),silent=T)
    return(try1)
    if (class(try1)=="try-error")
      return(matrix("",nrow=length(levels(X)),ncol=2*length(levels(D))+2))
    else return(try1)
  }
  else {
    try2 <- try(tblcontvars(Y,DS,d),silent=T)
    if (class(try2)=="try-error")
      return(matrix("",nrow=1,ncol=2*length(levels(D))+2))
    else return(try2)
  }
}
#maketable("race",test.ds,"ds")


tblcatvars.1d <- function(y,ds) {
  x <- ds[[pmatch(y,names(ds))]]
#  D <- ds[[pmatch(d,names(ds))]]
  Ns <- with(ds,table(x))
  props <- round(100*prop.table(Ns),1)
  v <- paste(as.vector(t(Ns))," (",as.vector(t(props)),"\\%)",sep="")
  c2 <- rbind("",cbind(matrix(v,nrow=length(levels(x)))))
  c1 <- c(sum(Ns),rep("",length(levels(x))))
  m1 <- cbind(c(paste("\\textbf{",y,"}",sep=""),paste("\\hspace*{5mm}",levels(x),sep="")),c1,c2)
  rownames(m1) <- NULL
  colnames(m1) <- NULL
  return(m1)
}
#tblcatvars.1d("sex",combined)

tblcontvars.1d <- function(y,ds) {
  x <- ds[[pmatch(y,names(ds))]]
  return(c(y,NmeanSD(x)))
}
#tblcontvars.1d("age",combined)

maketable.1d <- function(x,ds) {
  X <- ds[[pmatch(x,names(ds))]]
  if(is.factor(X)) return(tblcatvars.1d(x,ds))
  else return(tblcontvars.1d(x,ds))
}
