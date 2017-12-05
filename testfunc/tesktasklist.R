library(stringr)
library(future)

a <- future({
b<-1
while(length(b)>0){
  a <- shell("tasklist /v",intern=T)
  b<- as.vector(na.omit(str_extract(a,"mod[0-9]* - execute")))
}}) %plan% multiprocess

resolved(a)

