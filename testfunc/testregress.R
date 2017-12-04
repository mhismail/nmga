#read in continuous covariates
cotab<-readLines("cotab1")
cotab <- if(substr(cotab[1],1,5)=="TABLE") cotab[-1]
cotab <- strsplit(cotab,fixed=T)


cotab<-read.table("cotab1",skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)
catab<-read.table("catab1",skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)
patab <- read.table("patab1",skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)

ETAS <- names(patab)[grepl("ETA",names(patab))]
CO<- names(cotab)[!(grepl("ID",names(cotab)))]
CA <- names(catab)[!(grepl("ID",names(catab)))]

etatab <- patab[c(1,which(grepl("ETA",names(patab))))]

alltab <- right_join(cotab,catab)%>%right_join(etatab)

cor <- list()
cor <- lm (alltab[,3:4]~alltab[,2])

summary(cor)
coef(summary(cor))[, 2]

k = 0
corlist <- list()
for (i in ETAS){
  ETAi <- which(names(alltab)==i)
  for (j in CO){
    k=k+1
    PARAMi <- which(names(alltab)==j)
    cor <- lm (alltab[,ETAi]~alltab[,PARAMi])
    corlist[[k]]<-(c(Relationship="Linear",
            Model=paste0(j,"on",i),
            Estimate =round(as.numeric(cor[[1]][2]),3), 
            STDError = round(as.numeric(coef(summary(cor))[, 2][2] ),3)))
    }
} 

allcor<-data.frame()
allcor<- rbind(allcor, do.call(rbind, corlist))



center <- function (x){x-median(x)}
center2 <- function (x){x/median(x)}

linmod <- function (x,y,xlab="X",ylab="Y"){
  cor <- lm(y~ x)
  return(c(R="Linear",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
  
}
expmod <- function (x,y,xlab="X",ylab="Y"){
  if (any(is.nan(log(y)))) return( NULL)
  cor <- lm(log(y)~ x)
  return(c(R="Exp",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}
powmod <- function (x,y,xlab="X",ylab="Y"){
  if (any(is.nan(log(y))) | any(is.nan(log(x)))) return( NULL)
  cor <- lm(log(y)~ log(x))
  return(c(R="Power",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}
linmod(center(xpdb@Data$WT),xpdb@Data$CL)
expmod(center(xpdb@Data$WT),xpdb@Data$CL)
powmod(center2(xpdb@Data$WT),xpdb@Data$CL)
xpdb<-xpose.data(xpose.runno)

k=0
all <- list()
for (i in xpdb@Prefs@Xvardef$parms){
  y <- select_(xpdb@Data,i)
  for (j in xpdb@Prefs@Xvardef$covariates){
    x <- select_(xpdb@Data,j)
    if(!is.factor(x[,1])){
    k=k+1
    a<-rbind(
    linmod(center(x[,1]),y[,1],xlab=names(x),ylab=names(y)),
    expmod(center(x[,1]),y[,1],xlab=names(x),ylab=names(y)),
    powmod(center2(x[,1]),y[,1],xlab=names(x),ylab=names(y)))
    all[[k]] <- a
    }
  }
}

allcor<-data.frame()
allcor<- rbind(allcor, do.call(rbind, all))

