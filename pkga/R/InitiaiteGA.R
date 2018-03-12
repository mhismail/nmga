
#' InitiateGA
#'
#' Initiate the genetic algorithm on defined model features
#' @param nmods Numeric. The number of all possible models
#' @param nindiv Numeric. The number of individuals per generation
#' @param control Optional. Control stream template if model control stream does not exist yet. 
#' @param alltokens Optional. data frame containing tokens.
#' @param allmods  Optional. Data frame containing all model combinations.
#' @param seed Numeric. Seed number for reproducibility 
#' @export

InitiateGA <- function (nmods,nindiv,control=NULL,alltokens=NULL,allmods=NULL,seed=1){
  if (nindiv > nmods) nindiv = nmods
  set.seed(seed)
  pop <- sample(1:nmods,nindiv,replace=F)
  popmods <- paste0("models/All/mod",pop)
  for (i in popmods){
    CheckThenCreate(i,control,alltokens,allmods)
  }
  
  homepath <-getwd()
  for (i in popmods){
    relpath <- i
    if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
      RunModel(paste0(homepath,'/',relpath),basename(relpath))
    }
  }
  
  write.csv(data.frame(Generation=1,Inidvidual=pop,Dir=popmods), "GAgens.csv",row.names=F)
}


#' NextGA
#'
#' Determine the next generation of individuals
#' @param mods Character vector of paths to candidate models
#' @param control Control stream template if model control stream does not exist yet. 
#' @param alltokens Data frame containing tokens.
#' @param allmods  Optional. Data frame containing all model combinations.
#' @export

NextGA <- function (alltokens,control=NULL,allmods=NULL){
  GAgens <- read.csv("GAgens.csv")
  maxgen<-max(GAgens$Generation)
  currentgen <- maxgen+1
  if((currentgen%%5!=0 & currentgen%%5!=1)|currentgen<5){ # generation is not downhill gen
  mods <- filter(GAgens, Generation==maxgen)
  mods <- paste0("./",mods$Dir)
  a<-selection(mods)
  b<-crossover(a,2)
  c <- mutation(b,alltokens)
  
  zz=gzfile('Allmodsresults.csv.gz')   
  allmodsresults <- read.csv(zz)
  modlookup <- allmodsresults[-2:-9]
  names(modlookup)[1]<- "X"
  print(head(modlookup))
  newmodnumbers <- merge(modlookup,c)$X
  newmods <-  paste0("models/All/mod",newmodnumbers)
  print(currentgen)
    if (currentgen!=1){ # insert elites
      mods_last_gen <- filter(GAgens, Generation==(currentgen-1))
      mods_last_gen <- paste0("./",mods_last_gen$Dir)
      
      allmods1<-data.frame()
      a<-lapply(mods_last_gen,RetrieveResultsEach)
      allmods1<- rbind(allmods1, do.call(rbind, a))
      
      best <- unique(filter(allmods1, Fitness<=nthmin(allmods1$Fitness,2))$Number)
      best_path <-  paste0("models/All/mod",best)
      print(newmods)
      replace <- sample(1:length(newmods),2,replace = F)

      newmods[replace]<-best_path
      
    }

  for (i in newmods){
    CheckThenCreate(i,control,alltokens,allmods)
  }
  
  
  homepath <-getwd()
  for (i in unique(newmods)){
    relpath <- i
    if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
      RunModel(paste0(homepath,'/',relpath),basename(relpath))
    }
  }
  
  
  
  write.table(data.frame(Generation=maxgen+1,Inidvidual=newmodnumbers,Dir=newmods),
              "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
  return(currentgen) #return new max \gen to update buttons
  }
  if(currentgen%%5==0){ # downhill gen
    mods <- filter(GAgens, Generation==max(Generation))
    mods <- paste0("./",mods$Dir)
    
    allmods1<-data.frame()
    a<-lapply(mods,RetrieveResultsEach)
    allmods1<- rbind(allmods1, do.call(rbind, a))
    
    best <- filter(allmods1, Fitness==min(allmods1$Fitness))
    newmods <- c(downhill(best[1,-(1:9)],allmods))
    newmodnumbers<- sub("models/All/mod","",newmods)

    
    for (i in newmods){
      CheckThenCreate(i,control,alltokens,allmods)
    }
    
    
    homepath <-getwd()
    for (i in unique(newmods)){
      relpath <- i
      if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
        RunModel(paste0(homepath,'/',relpath),basename(relpath))
      }
    }
    
    maxgen<-max(read.csv("GAgens.csv")$Generation)
    
    write.table(data.frame(Generation=maxgen+1,Inidvidual=newmodnumbers,Dir=newmods),
                "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
    return(currentgen) #return new max \gen to update buttons

  }
  
  if(currentgen%%5==1 & currentgen>5){ # generation after downhill

    downhill <- filter(GAgens, Generation==currentgen-1)
    downhill <- paste0("./",downhill$Dir)
    allmods1<-data.frame()
    a<-lapply(downhill,RetrieveResultsEach)
    allmods1<- rbind(allmods1, do.call(rbind, a))
    
    best <- filter(allmods1, Fitness==min(allmods1$Fitness))$Number[1]
    best_path <-  paste0("models/All/mod",best)
    

    mods <- filter(GAgens, Generation==currentgen-2)
    mods <- paste0("./",mods$Dir)
    a<-selection(mods)
    b<-crossover(a,2)
    c <- mutation(b,alltokens)
    
    
    zz=gzfile('Allmodsresults.csv.gz')   
    allmodsresults <- read.csv(zz)
    modlookup <-allmodsresults[-2:-9]
    names(modlookup)[1]<- "X"
    
    newmodnumbers <- merge(modlookup,c)$X
    newmods <-  paste0("models/All/mod",newmodnumbers)
    print(newmods)
    replace <- sample(1:length(newmods),1)
    newmods[replace]<-best_path
    
    for (i in newmods){
      CheckThenCreate(i,control,alltokens,allmods)
    }
    
    
    homepath <-getwd()
    for (i in unique(newmods)){
      relpath <- i
      if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
        RunModel(paste0(homepath,'/',relpath),basename(relpath))
      }
    }
    
    
    
    write.table(data.frame(Generation=maxgen+1,Inidvidual=newmodnumbers,Dir=newmods),
                "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
    return(currentgen) #return new max \gen to update buttons
    
  }
  
  
  
}

downhill <- function (mod,allmods=NULL){

  
  similarmods <- which(apply(allmods, 1, function(x) sum( x ==mod)>=(dim(allmods)[2]-1)))
  newmods<- paste0("models/All/mod",similarmods)

 
}

nthmin <- function (x,n){
  length <- length(x)
  sort(x,partial=n)[n]
}

# 
# RandomNextGA <- function (nmods,nindiv){
#   if (nindiv > nmods) nindiv = nmods
#   pop <- sample(1:nmods,nindiv,replace=F)
#   popmods <- paste0("models/All/mod",pop)
#   
#   homepath <-getwd()
#   for (i in popmods){
#     relpath <- i
#     unlink(paste0(homepath,'/',relpath,"/results"),recursive = TRUE)
#     RunModel(paste0(homepath,'/',relpath),basename(relpath))
#   }
#   maxgen<-max(read.csv("GAgens.csv")$Generation)
#   
#   write.table(data.frame(Generation=maxgen+1,Inidvidual=pop,Dir=popmods),
#               "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
#   return(maxgen+1) #return new max \gen to update buttons
# }