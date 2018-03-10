#' Copy Control Stream 
#'
#' This function copies a control stream to a new directory
#' @param copy A character string containing the path to control stream.
#' @param copyTo A character string containing the directory to copy to.
#' @param control Optional. Control stream template if model control stream does not exist yet. 
#' @param alltokens Optional. data frame containing tokens.
#' @param allmods  Optional. Data frame containing all model combinations.
#' @export

CopyModel <- function (copy,copyTo,control=NULL,alltokens=NULL,allmods=NULL){
  
  if(!is.null(copyTo)){
    CheckThenCreate(copy,control,alltokens,allmods)
    if(copyTo %in% list.dirs(paste0(getwd(),"/models/"),recursive = F,full.names = F)){
      file.copy(copy, paste0("models/",copyTo), recursive=TRUE)
      
      mod <-read.csv(paste0(copy,"/mod.csv"),as.is=T)
      i <- mod$Number[1]
      a <-mutate(mod,Path=paste0("models/",copyTo,"/mod",i,"/mod.ctl"))
      write.csv(a,paste0("models/",copyTo,"/mod",i,"/mod.csv"),row.names = F)
    }
  }
}

#' CheckThenCreate
#'
#' This function checks if a model control stream exists, and creates it based on control stream template and tokens if it doesn't
#' @param path A character string containing the path to control stream.
#' @param control Optional. Control stream template if model control stream does not exist yet. 
#' @param alltokens Optional. data frame containing tokens.
#' @param allmods  Optional. Data frame containing all model combinations.
#' @export

CheckThenCreate<- function (path,control=NULL,alltokens=NULL,allmods = NULL )({
  if (!dir.exists(path)){
    mod<-as.numeric(sub(".*mod([0-9]+).*","\\1",path))
    x <- control
    x <- gsub("\\$DATA\\s*", "\\$DATA ../../",x) 
    x<-CreateModel(x,alltokens,allmods[mod,],names(allmods))
    dir.create(paste0("models/All/mod",mod))
    writeLines(x,paste0("models/All/mod",mod,"/mod.ctl"))
    
    maxtheta <- gsub("[THETA\\(\\)]", "", regmatches(x, gregexpr("THETA\\(.*?\\)", x))[[1]])
    maxtheta<-max(as.numeric(c(maxtheta,0)),na.rm=T)
    
    maxeta <- gsub("[ETA\\(\\)]", "", regmatches(x, gregexpr("\\bETA\\(.*?\\)", x))[[1]])
    maxeta<-max(as.numeric(c(maxeta,0)),na.rm=T)
    
    maxeps <- gsub("[EPS\\(\\)]", "", regmatches(x, gregexpr("\\bEPS\\(.*?\\)", x))[[1]])
    maxeps <-max(as.numeric(c(maxeps,0)),na.rm=T)
    
    a <-cbind(data.frame(Number=mod,Path=paste0("models/All/mod",mod,"/mod.ctl")),OFV="",Fitness="",S="",C="",NTHETA=maxtheta,NETA=maxeta,NEPS = maxeps, allmods[mod,])
    
    
    write.csv(a,paste0("models/All/mod",mod,"/mod.csv"),row.names = F)
  }else(return(TRUE))
  
})


#' CreateModel
#'
#' Creates syntactically correct control stream based on template and tokens
#' CreateModel takes a control stream, a token data frame (used as a look up)
#' to translate phenotypes to model syntax, a model phenotype, and tokengroups
#' it returns a new control stream
#' @param controlstream A character string of control stream template.
#' @param alltokens Data frame containing tokens.
#' @param phenotype Data frame containing specific model features (tokensets).
#' @param tokengroups Data frame containing general model features (tokengroups).
#' @export

CreateModel<- function(controlstream,alltokens,phenotype,tokengroups){ 
  x<-controlstream 
  
  for (j in 1:length(phenotype)){
    selectedTokenGroup <- tokengroups[j] 
    selectedTokenSet <- phenotype[1,j] 
    
    tokengrouptext <- paste0("\\{",selectedTokenGroup,"\\}") 
    tokenlist <- alltokens[(alltokens$tokengroup==selectedTokenGroup & alltokens$tokenset==selectedTokenSet),3]# filter(alltokens,tokengroup==selectedTokenGroup,tokenset==selectedTokenSet)$token 
    tokenlist[tokenlist=="N/A"]<-"" 
    numberedtokengroups <- as.numeric(gsub(paste0("\\{",selectedTokenGroup,":|\\}"),  
                                           "", 
                                           regmatches(x, gregexpr(paste0("\\{",selectedTokenGroup,":([0-9])\\}"), x))[[1]])) 
    
    for (k in numberedtokengroups){ 
      x<- gsub(paste0("\\{",selectedTokenGroup,":",k,"\\}"),tokenlist[k],x) 
    } 
    
    
    maxtheta <- gsub("[THETA\\(\\)]", "", regmatches(x, gregexpr("THETA\\(.*?\\)", x))[[1]]) 
    maxtheta<-max(as.numeric(c(maxtheta,0)),na.rm=T) 
    maxeta <- gsub("[ETA\\(\\)]", "", regmatches(x, gregexpr("\\bETA\\(.*?\\)", x))[[1]]) 
    maxeta<-max(as.numeric(c(maxeta,0)),na.rm=T) 
    maxeps <- gsub("[EPS\\(\\)]", "", regmatches(x, gregexpr("\\bEPS\\(.*?\\)", x))[[1]]) 
    maxeps<-max(as.numeric(c(maxeps,0)),na.rm=T) 
    
    
    for (k in 1:length(tokenlist)){ 
      
      x<-sub(tokengrouptext,tokenlist[k],x) 
      
    }
    if (grepl("THETA\\(\\[([0-9])\\]\\)",x)){
      x<-gsubfn("THETA\\(\\[([0-9])\\]\\)", x~ paste0("THETA(",as.numeric(x)+maxtheta,")"),x,backref = -1,engine="R") 
    }
    if(grepl("ETA\\(\\[([0-9])\\]\\)",x)){
      x<-gsubfn("ETA\\(\\[([0-9])\\]\\)", x~ paste0("ETA(",as.numeric(x)+maxeta,")"),x,backref = -1,engine="R")
    }
    if(grepl("EPS\\(\\[([0-9])\\]\\)",x)){
      x<-gsubfn("EPS\\(\\[([0-9])\\]\\)", x~ paste0("EPS(",as.numeric(x)+maxeps,")"),x,backref = -1,engine="R")
    }
    
    maxeta <- gsub("[ETA\\(\\)]", "", regmatches(x, gregexpr("\\bETA\\(.*?\\)", x))[[1]]) 
    maxeta<-max(as.numeric(c(maxeta,0)),na.rm=T) 
    print(maxeta)
  } 
  x<- gsub("ALLETAS",paste0("ETA",1:maxeta,collapse = " "),x)
  return(x) 
} 
