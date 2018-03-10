#' p1cmt
#'
#' Create 1 compartment diagram
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#' 

p1cmt <- function(compartment,cmtn,ncmts){
  grViz(paste0(' digraph  {
               graph [layout = circo]
               
               node [margin=0 width=0.5 shape=circle]
               ',compartment,' ->  ',compartment,'2  [xlabel=" K',cmtn,ncmts+1,'"];
               ',compartment,'2 -> ',compartment,'   [xlabel="K',ncmts+1,cmtn,'"];
}'))}

#' p2cmt
#'
#' Create 2 compartment diagram
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#' 

p2cmt <- function(compartment,cmtn,ncmts){
  
  grViz(paste0(' digraph  {
               forcelabels=true;
               
               graph [layout = circo]
               
               node [shape=circle]
               ',compartment,' ->  ',compartment,'2  [xlabel=" K',cmtn,ncmts+1,'"];
               ',compartment,'2 -> ',compartment,'   [xlabel="K',ncmts+1,cmtn,'"];
               ',compartment,' ->  ',compartment,'3  [headlabel="K',cmtn,ncmts+2,'", labeldistance=3, labelangle=-20];
               ',compartment,'3 -> ',compartment,'   [headlabel="K',ncmts+2,cmtn,'", labeldistance=2.4, labelangle=-30];
               ',compartment,'3 -> ',compartment,'2 [style="invis"]
               
}'))}


#' p3cmt
#'
#' Create 3 compartment diagram
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#' 
p3cmt <- function(compartment,cmtn,ncmts){
  
  grViz(paste0(' digraph  {
               forcelabels=true;
               
               graph [layout = circo]
               
               node [shape=circle]
               ',compartment,' ->  ',compartment,'2  [xlabel="K',cmtn,ncmts+1,'"];
               ',compartment,'2 -> ',compartment,'   [xlabel="K',ncmts+1,cmtn,'"];
               ',compartment,' ->  ',compartment,'3  [headlabel="K',cmtn,ncmts+2,'", labeldistance=3, labelangle=-20];
               ',compartment,'3 -> ',compartment,'   [headlabel="K',ncmts+2,cmtn,'", labeldistance=2, labelangle=-50];
               ',compartment,' ->  ',compartment,'4  [headlabel="K',cmtn,ncmts+3,'", labeldistance=4, labelangle=-50];
               ',compartment,'4 -> ',compartment,'   [headlabel="K',ncmts+3,cmtn,'", labeldistance=3, labelangle=-20];
               
}'))}

#' p1cmttext
#'
#' Create 1 compartment text strings to be inserted in template
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#' 

p1cmttext <- function(compartment,cmtn,ncmts){
  a<- list()
  a[[1]]<-paste0("COMP(",substr(compartment,1,7),"2)")
  a[[2]]<-paste0(
    "TVQ",cmtn,ncmts+1," = THETA([1])
    Q",cmtn,ncmts+1,"=TVQ",cmtn,ncmts+1,"
    TVV",cmtn,ncmts+1,"=THETA([2])
    V",cmtn,ncmts+1,"=TVV",cmtn,ncmts+1,"
    K",cmtn,ncmts+1,"=Q",cmtn,ncmts+1,"/V
    K",ncmts+1,cmtn,"=Q",cmtn,ncmts+1,"/V",cmtn,ncmts+1)
  a
}



#' p2cmttext
#'
#' Create 2 compartment text strings to be inserted in template
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#'

p2cmttext <- function(compartment,cmtn,ncmts){
  a<- list()
  
  a[[1]]<-paste0("COMP(",substr(compartment,1,7),"2)\nCOMP(",substr(compartment,1,7),"3)")
  a[[2]]<-paste0(
    "TVQ",cmtn,ncmts+1," = THETA([1])
    Q",cmtn,ncmts+1,"=TVQ",cmtn,ncmts+1,"
    TVV",cmtn,ncmts+1,"=THETA([2])
    V",cmtn,ncmts+1,"=TVV",cmtn,ncmts+1,"
    K",cmtn,ncmts+1,"=Q",cmtn,ncmts+1,"/V
    K",cmtn,ncmts+1,"=Q",cmtn,ncmts+1,"/V",cmtn,ncmts+1,"
    TVQ",cmtn,ncmts+2," = THETA([3])
    Q",cmtn,ncmts+2,"=TVQ",cmtn,ncmts+2,"
    TVV",cmtn,ncmts+2,"=THETA([4])
    V",cmtn,ncmts+2,"=TVV",cmtn,ncmts+2,"
    K",cmtn,ncmts+2,"=Q",cmtn,ncmts+2,"/V
    K",ncmts+2,cmtn,"=Q",cmtn,ncmts+2,"/V",cmtn,ncmts+2)
  a
}


#' p3cmttext
#'
#' Create 3 compartment text strings to be inserted in template
#' @param compartment Character string containing name of center compartment
#' @param cmtn Numeric. Compartment number of center compartment
#' @param ncmts Numeric. Number of compartments already present
#' @export 
#'


p3cmttext <- function(compartment,cmtn,ncmts){
  a<- list()
  
  a[[1]]<-paste0("COMP(",substr(compartment,1,7),"2)\nCOMP(",substr(compartment,1,7),"3)")
  a[[2]]<-paste0(
    "TVQ",cmtn,ncmts+1," = THETA([1])
    Q",cmtn,ncmts+1,"=TVQ",cmtn,ncmts+1,"
    TVV",cmtn,ncmts+1,"=THETA([2])
    V",cmtn,ncmts+1,"=TVV",cmtn,ncmts+1,"
    K",cmtn,ncmts+1,"=Q",cmtn,ncmts+1,"/V
    K",cmtn,ncmts+1,"=Q",cmtn,ncmts+1,"/V",cmtn,ncmts+1,"
    TVQ",cmtn,ncmts+2," = THETA([3])
    Q",cmtn,ncmts+2,"=TVQ",cmtn,ncmts+2,"
    TVV",cmtn,ncmts+2,"=THETA([4])
    V",cmtn,ncmts+2,"=TVV",cmtn,ncmts+2,"
    K",cmtn,ncmts+2,"=Q",cmtn,ncmts+2,"/V
    K",ncmts+2,cmtn,"=Q",cmtn,ncmts+2,"/V",cmtn,ncmts+2,"
    TVQ",cmtn,ncmts+3," = THETA([5])
    Q",cmtn,ncmts+3,"=TVQ",cmtn,ncmts+3,"
    TVV",cmtn,ncmts+3,"=THETA([6])
    V",cmtn,ncmts+3,"=TVV",cmtn,ncmts+3,"
    K",cmtn,ncmts+3,"=Q",cmtn,ncmts+3,"/V
    K",ncmts+3,cmtn,"=Q",cmtn,ncmts+3,"/V",cmtn,ncmts+3)
  a
}