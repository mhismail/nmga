

#' AllCombos
#'
#' Create all possible combinations of tokens sets. alltokens argument must be 
#'  a dataframe with $tokengroup, $tokenset, $token
#' @param alltokens A data.frame containing model features tokens
#' @export

AllCombos <- function (alltokens){
  tokenset1 <- alltokens
  
  #create column of group concatenated with set 
  tokenset1$groupset <- paste0(tokenset1$tokengroup,"-",tokenset1$tokenset) 
  
  #group by token group and find lenght of each group
  tokenset1_grouped <- group_by(tokenset1,tokengroup)
  nmodspergroup <- summarize(tokenset1_grouped,n = length(unique(groupset)))
  
  #total number of mods is product of number of tokens sets in each group
  total_mods_n <- prod(nmodspergroup$n)
  
  #initialize empty list
  list <- list()
  
  #create list of vectors, one for each token group
  for (i in 1:length(unique(tokenset1$tokengroup))){
    list[[i]]<- unique(filter(tokenset1,tokengroup==unique(alltokens$tokengroup)[i])$tokenset)
  }
  names(list)<- unique(alltokens$tokengroup)
  
  # create data frame of every possible combination of token sets
  allmods <-do.call(expand.grid,list)
  
  
  return (allmods)
}