
#' RetrieveResults
#'
#' Retrieve results if results directory exists. returns dataframe with model number,
#' OFV, run success status, and covariance success status
#' This version of retrieveresults is used for the "Allmods" section
#' @param directory Character string containing project directory
#' @export

RetrieveResults <- function(directory){
  if (file.exists(paste0(directory,"/results/raw_results_mod.csv"))){
    results<- read.csv(paste0(directory,"/results/raw_results_mod.csv"),as.is = T,stringsAsFactors = F)
    results <- data.frame(Number=sub(".*mod([0-9]+).*","\\1",directory),
                          OFV=results$ofv[1],
                          S=results$minimization_successful[1],
                          C=results$covariance_step_successful[1])
    return(results)
    
    
  }
}

# 
#' RetrieveResultsEach
#'
#' Retrieve results if results directory exists. returns dataframe with model number,
#' OFV, run success status, and covariance success status
#' This version of retrieve results is used in the GA portion of the interface,
#' as well as user created subdirectories outside of "All".
#' @param directory Character string containing project directory
#' @export

RetrieveResultsEach <- function(directory){
  pheni <- read.csv(paste0(directory,"/mod.csv"),as.is = T,stringsAsFactors = F)
  if (file.exists(paste0(directory,"/results/raw_results_mod.csv"))){
    results<- read.csv(paste0(directory,"/results/raw_results_mod.csv"),as.is = T,stringsAsFactors = F)
    
    pheninew <- mutate(pheni, OFV=results$ofv[1],
                       S=results$minimization_successful[1],
                       C=results$covariance_step_successful[1],
                       Fitness=fitness(cbind(results,pheni)))
    if (!(identical(pheninew,pheni))) write.csv(pheninew,paste0(directory,"/mod.csv"),row.names = F)
    pheni <- pheninew
    
  }
  return(pheni)
}
