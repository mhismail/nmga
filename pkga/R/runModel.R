#' RunModel
#'
#' This function executes mod.ctl in the path specified via PsN.
#' @param directory A character string containing the path to folder in which mod.ctl is located
#' @param model Optional. A character string providing title of running process
#' @export

RunModel <- function(directory, model = NULL, wait = T, exit = T){
  if(!dir.exists(directory)) {
    alert(paste("Directory", directory, "not found")) 
    return()
  }
  
  if(!file.exists(paste0(directory, "/mod.ctl"))) {
    alert(paste("mod.ctl not found in", directory)) 
    return()
  }
  
  
  if (exit == T) close <- " & exit" else close <- ""
  
  shell(paste0(' start cmd /k  "cd ',
               directory,
               ' &title ', 
               model,
               '&execute -directory=results -clean=3 mod.ctl', 
               close),
        wait = wait)
  
  write.table(data.frame(directory, Sys.time()),
              "modelruntemp.csv", 
              sep = ",",
              row.names = F,
              col.names = F,  
              append = T)
}
