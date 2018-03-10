
#' Selection
#'
#' Performs tournament style selection
#' @param candidates Character vector of paths to candidate models
#' @return Returns character vector of paths to selected models
#' @export

selection <- function(candidates){
  pool <- c()
  
  for (i in 1:length(candidates)){
    if (file.exists(paste0(candidates[i],"/results/raw_results_mod.csv"))){
      mod <-  read.csv(paste0(candidates[i],"/results/raw_results_mod.csv"))
      modres <-read.csv(paste0(candidates[i],"/mod.csv"))
      modobj <- modres$Fitness
    }else modobj = NA
    if (is.na(modobj)) modobj = 9999999
    
    competitormod <- sample(candidates[-i],1)
    if (file.exists(paste0(competitormod,"/results/raw_results_mod.csv"))){
      competitor <-  read.csv(paste0(competitormod,"/results/raw_results_mod.csv"))
      cmodres <- read.csv(paste0(competitormod,"/mod.csv"))
      competitorobj <- cmodres$Fitness
    }else competitorobj = NA
    if (is.na(competitorobj)) competitorobj = 9999999
    
    if (modobj <=competitorobj) {pool <- c(pool,candidates[i])} else pool <- c(pool,competitormod)
    
  }
  return (pool)
}



#' crossover
#'
#' Performs n-point crossover on selected individuals
#' @param selectedmods Character vector of paths to selected models
#' @return Returns data.frame of model phenotypes
#' @export
#' 
crossover<- function(selectedmods,ncrossoverpoints){
  
  #Split population into 2 equal parts
  n <- length(selectedmods)
  sample <-sample(n,n/2)
  mate1 <- selectedmods [sample]
  mate2 <- selectedmods [-sample]
  print(selectedmods)
  
  offspring <- list()
  
  for (i in 1:(n/2)){
    mate1i <-  read.csv(paste0(mate1[i],"/mod.csv"),as.is = T)[-(1:9)] #drop col 1:9, phenotype is 9+
    mate2i <-  read.csv(paste0(mate2[i],"/mod.csv"),as.is = T)[-(1:9)]
    
    #Create two off springs for each pair and store in list
    crossoverpoints <- sort(sample(0:length(mate1i),ncrossoverpoints))
    if(max(crossoverpoints)<length(mate1i)){
      offspring[[i*2-1]] <- data.frame(mate1i[0:crossoverpoints[1]],mate2i[(crossoverpoints[1]+1):crossoverpoints[2]],mate1i[(crossoverpoints[2]+1):length(mate1i)])
      offspring[[i*2]] <- data.frame(mate2i[0:crossoverpoints[1]],mate1i[(crossoverpoints[1]+1):crossoverpoints[2]],mate2i[(crossoverpoints[2]+1):length(mate1i)])
      
    }else{
      offspring[[i*2-1]] <- data.frame(mate1i[0:crossoverpoints[1]],mate2i[(crossoverpoints[1]+1):crossoverpoints[2]])
      offspring[[i*2]] <- data.frame(mate2i[0:crossoverpoints[1]],mate1i[(crossoverpoints[1]+1):crossoverpoints[2]])
      
    }
    
  }
  
  #list of offspring to df
  alloffspring <- data.frame()
  alloffspring<- rbind(alloffspring, do.call(rbind, offspring))
  print(alloffspring)
  return(alloffspring)
  
  
}

#' Mutation
#'
#' Performs mutation
#' @param crossoveredmods Data frame of models from crossover step
#' @param alltokens Data frame containing tokens.
#' @param options List containing mutation options (pmutation)
#' @return Returns data.frame of model phenotypes 
#' @export
mutation <- function (crossoveredmods, alltokens, options= list(pmutation=.07)){
  pmutation <- options$pmutation
  for (i in 1:dim(crossoveredmods)[2]){
    gene <- names(crossoveredmods)[i]
    for (j in  1:dim(crossoveredmods)[1]){
      if(runif(1)<=pmutation){
        crossoveredmods[j,i]<- sample(unique(filter(alltokens,tokengroup==gene)$tokenset),1)
      }
    }
  }
  mutatedmods <-crossoveredmods
  return(mutatedmods)
}
