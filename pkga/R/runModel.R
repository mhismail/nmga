#' RunModel
#'
#' This function executes mod.ctl in the path specified via PsN.
#' @param directory A character string containing the path to folder in which mod.ctl is located
#' @param model Optional. A character string providing title of running process
#' @export

RunModel <- function(directory,model=NULL){
  shell(paste0(' start cmd /k  "cd ',directory,' &title ', model,'&execute -directory=results -clean=3 mod.ctl'),wait = T)
  write.table(data.frame(directory,Sys.time()),
              "modelruntemp.csv", sep = ",",row.names = F,col.names = F,  append = T)
}