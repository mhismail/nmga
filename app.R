#STABLE-test

#check if packages are installed, and if not install them
if (!require("shiny")) install.packages("shiny")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyBS")) install.packages("shinyBS")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("DT")) install.packages("DT")
if (!require("xpose4")) install.packages("xpose4")
if (!require("gsubfn")) install.packages("gsubfn")
if (!require("DiagrammeR")) install.packages("DiagrammeR")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("foreach")) install.packages("foreach")
if (!require("doSNOW")) install.packages("doSNOW")
if (!require("future")) install.packages("future")
# if (!require("git2r")) install.packages("git2r")





#Set plot theme
theme_set(theme_bw())

# Define functions --------------------------------------------------------

create.model<- function(controlstream,alltokens,phenotype,tokengroups){ 
  x<-controlstream 
  for (j in 1:length(phenotype)){ 
    selectedTokenGroup <- tokengroups[j] 
    selectedTokenSet <- phenotype[j] 
    
    tokengrouptext <- paste0("\\{",selectedTokenGroup,"\\}") 
    tokenlist <- filter(alltokens,tokengroup==selectedTokenGroup,tokenset==as.character(selectedTokenSet))$token 
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
    
    for (k in 1:length(tokenlist)){ 
      
      
      x<-sub(tokengrouptext,tokenlist[k],x) 
      x<-gsubfn("THETA\\(\\[([0-9])\\]\\)", x~ paste0("THETA(",as.numeric(x)+maxtheta,")"),x) 
      
      x<-gsubfn("ETA\\(\\[([0-9])\\]\\)", x~ paste0("ETA(",as.numeric(x)+maxeta,")"),x) 
    } 
  } 
  return(x) 
} 

# Function Name : run.model-------------------------------------------------------------------------------------
# Description: Start execution of NONMEM run via PsN
# Arguments: directory(required)-directory in which "mod.ctl" is located 
#            model(optional)- model name (to be used for identifying in running task) 
# Example: run.model("models/All/mod1","mod1")

run.model <- function(directory,model=NULL){
  shell(paste0(' start cmd /k  "cd ',directory,' &title ', model,'&execute -directory=results -clean=3 mod.ctl'),wait = T)
  write.table(data.frame(directory,Sys.time()),
              "modelruntemp.csv", sep = ",",row.names = F,col.names = F,  append = T)
}



# Function Name : allcombos-------------------------------------------------------------------------------------
# Description: create all possible combinations of tokens sets. alltokens argument must be 
#              a dataframe with $tokengroup, $tokenset, $token
# Arguments: alltokens (required) - dataframe with tokengroup, tokenset, and token columns
# Example: 
# alltokens <- read.csv("alltokens.csv",as.is=T)
# allcombos (alltokens)
# Output example:

# tokengroup1 | tokengroup2 | tokengroup3
# ------------+-------------+-------------
# tokenset11  | tokenset 21 | tokenset31
# ------------+-------------+-------------
# tokenset12  | tokenset 21 | tokenset31
# etc.

allcombos <- function (alltokens){
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


# Function Name : selection-------------------------------------------------------------------------------------
# Description: Performs tournament selection on group of models, returns models in tabular format
# Arguments: candidates (required) - character vector of paths to models
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

crossover<- function(selectedmods,ncrossoverpoints){
  
  #Split population into 2 equal parts
  n <- length(selectedmods)
  sample <-sample(n,n/2)
  mate1 <- selectedmods [sample]
  mate2 <- selectedmods [-sample]

  offspring <- list()
  
  for (i in 1:(n/2)){
    mate1i <-  read.csv(paste0(mate1[i],"/mod.csv"),as.is = T)[-(1:8)] #drop col 1:8, phenotype is 8+
    mate2i <-  read.csv(paste0(mate2[i],"/mod.csv"),as.is = T)[-(1:8)]
    
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
  return(alloffspring)
  
  
}

mutation <- function (crossoveredmods, alltokens, options= list(pmutation=.05)){
  pmutation <- options$pmutation
  for (i in 1:dim(crossoveredmods)[2]){
    gene <- names(crossoveredmods)[i]
    for (j in  1:dim(crossoveredmods)[1]){
      if(runif(1)<=pmutation){
        print(alltokens)
        crossoveredmods[j,i]<- sample(unique(filter(alltokens,tokengroup==gene)$tokenset),1)
      }
    }
  }
  mutatedmods <-crossoveredmods
  return(mutatedmods)
}

initiateGA <- function (nmods,nindiv){
  
  if (nindiv > nmods) nindiv = nmods
  pop <- sample(1:nmods,nindiv,replace=F)
  popmods <- paste0("models/All/mod",pop)
  
  homepath <-getwd()
  for (i in popmods){
    relpath <- i
    if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
      run.model(paste0(homepath,'/',relpath),basename(relpath))
    }
  }
  
  write.csv(data.frame(Generation=1,Inidvidual=pop,Dir=popmods), "GAgens.csv",row.names=F)
}

nextGA <- function (mods,alltokens){
  a<-selection(mods)
  b<-crossover(a,2)
  c <- mutation(b,alltokens)
  modlookup <- read.csv("allmods.csv")
  newmodnumbers <- merge(modlookup,c)$X
  newmods <-  paste0("models/All/mod",newmodnumbers)

  homepath <-getwd()
  for (i in unique(newmods)){
    relpath <- i
    if(!dir.exists(paste0(homepath,'/',relpath,"/results"))){
      run.model(paste0(homepath,'/',relpath),basename(relpath))
    }
    }
  
  maxgen<-max(read.csv("GAgens.csv")$Generation)
  
  write.table(data.frame(Generation=maxgen+1,Inidvidual=newmodnumbers,Dir=newmods),
              "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
  return(maxgen+1) #return new max \gen to update buttons
  

}

randomnextGA <- function (nmods,nindiv){
  if (nindiv > nmods) nindiv = nmods
  pop <- sample(1:nmods,nindiv,replace=F)
  popmods <- paste0("models/All/mod",pop)
  
  homepath <-getwd()
  for (i in popmods){
    relpath <- i
    unlink(paste0(homepath,'/',relpath,"/results"),recursive = TRUE)
    run.model(paste0(homepath,'/',relpath),basename(relpath))
  }
  maxgen<-max(read.csv("GAgens.csv")$Generation)
  
  write.table(data.frame(Generation=maxgen+1,Inidvidual=pop,Dir=popmods),
              "GAgens.csv", sep = ",",row.names = F,col.names = F,  append = T)
  return(maxgen+1) #return new max \gen to update buttons
}

retrieveresults <- function(directory){
  print(directory)
  print(file.exists(paste0(directory,"/results/raw_results_mod.csv")))
  if (file.exists(paste0(directory,"/results/raw_results_mod.csv"))){
    results<- read.csv(paste0(directory,"/results/raw_results_mod.csv"),as.is = T,stringsAsFactors = F)
    results <- data.frame(Number=sub(".*mod([0-9]+).*","\\1",directory),
                          OFV=results$ofv[1],
                    S=results$minimization_successful[1],
                    C=results$covariance_step_successful[1])
    return(results)
    

    }
}


p1cmt <- function(compartment,cmtn,ncmts){
  grViz(paste0(' digraph  {
               graph [layout = circo]
               
               node [margin=0 width=0.5 shape=circle]
               ',compartment,' ->  ',compartment,'2  [xlabel=" K',cmtn,ncmts+1,'"];
               ',compartment,'2 -> ',compartment,'   [xlabel="K',ncmts+1,cmtn,'"];
}'))}

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
  if (any(is.nan(log(y))) | any(is.nan(log(x)))| any(is.infinite(log(x)))) return( NULL)
  cor <- lm(log(y)~ log(x))
   return(c(R="Power",Model=paste0(xlab,"on",ylab),Estimate=round(as.numeric(cor[[1]][2]),3),Pval=signif(summary(cor)$coefficients[,4][2] ,3)))
}

center <- function (x){x-median(x)}
center2 <- function (x){x/median(x)}

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  
  useShinyjs(), 

  actionButton("viewmod","View Models"),
  column(6,
    column(12,offset=2,
    actionButton("loadproj","Choose directory"),
    actionButton("saveproj","Save Project")),
    tabsetPanel(id="tabs",tabPanel("Control Stream",
    textAreaInput("ace",label = NULL,width="100%",cols=NULL,value="Select a directory with a single .ctl file")),
    tabPanel("Preview",
             textAreaInput("previewtext",label = NULL,width="100%",cols=NULL,value="")),
    tabPanel("Data",DT::dataTableOutput("data")))),
  column(6,
          column(12,
                 align="center",
                 actionButton("openGA","Initiate Genetic Algorithm"),
                 actionButton("GAsettings",icon("gear","fa-2x")),
                 radioGroupButtons("tokentype",NULL,choices = c("Covariate","ETA","EPS","Structure", "Custom"),selected = NULL)),
         column(12,
          column(4,
                 awesomeRadio(inputId = "tokengroupinput","Token Group", choices = c("")),
                 actionButton("deletetokengroup","Delete",style="background-color: #b71616; color: white;")
                 ),
          column(4,
                 awesomeRadio(inputId = "tokensetinput","Token Set", choices = c("")),
                 actionButton("preview","Preview"),
                 actionButton("deletetokenset","Delete",style="background-color: #b71616; color: white;")
                 
          ),
          column(4,
                 awesomeRadio(inputId = "tokeninput","Token", choices = c(""))
          ),
          column(12,
         div(id = "tokenedit"
         ),
         div(id = "Covariate",
             class = "tokentype",
             selectizeInput(
               "cov",
               "Select Covariate",
               NULL,
               selected = NULL,
               multiple = F,
               options = list(maxOptions =
                                100,create=T)
             ),
                checkboxGroupButtons("covtypes","Select Covariate relationships",choices= c("None", "Linear","Power","Exponential","Proportional") ),
                column(6,
                materialSwitch(inputId = "center", label = "Center Covariate (Median)", status = "success"),
             
                actionButton("addcovs","Add Selected Token Sets")),
             column(6,"Required {Token Group} in:", HTML("<h4 id='covrequiredpk'>$PK</h4><h4 id='covrequiredtheta'>$THETA</h4>"),class="required-holder")
          ),
          div(id = "ETA",
              class = "tokentype",
              checkboxGroupButtons("etatypes","Select IIV relationships",choices= c("None","Normal","Logarithmic","Normal (Proportional)") ),
              actionButton("addetas","Add Selected Token Sets")

          ),
          div(id = "EPS",
              class = "tokentype",
              checkboxGroupButtons("epstypes","Select Residual Error Models",choices= c("Additive","Proportional","Additive + Proportional") ),
              actionButton("addeps","Add Selected Token Sets")
          ),
         div(id = "Structure",
             class = "tokentype",
             radioGroupButtons("strtypes","Select Structural Component",choices= c("Compartments","Tlag","Michaelis-Menten") ),
             fluidRow(column(5,
             selectizeInput(
               "selectcmt",
               "Select Compartment",
               NULL,
               selected = NULL,
               multiple = F,
               options = list(maxOptions =
                                100,create=T)
             ),
             radioGroupButtons("ncmt","Additional Compartments",choices= c("1","2","3") ),
             actionButton("addstr","Add Selected Token Sets")),
             column(7,
             grVizOutput("distributiondiagram")))
         ))), class="token-div-holder"
  ),
  HTML('
        <script src="https://code.jquery.com/ui/1.12.0/jquery-ui.min.js"
       integrity="sha256-eGE6blurk5sHj+rmkfsGYeKyZx3M4bG+ZlFyA7Kns7E="
       crossorigin="anonymous"></script>
       <link rel="stylesheet" href="https://code.jquery.com/ui/1.12.0/themes/smoothness/jquery-ui.css">'),
  includeScript("www/jquery.highlighttextarea.js"),
  includeScript("www/drag.js"),
  includeCSS("www/style.css"),
  includeCSS("www/jquery.highlighttextarea.min.css"),
  includeCSS("www/font-awesome-4.7.0/css/font-awesome.min.css")
    )


# Server ------------------------------------------------------------------

server <- function(input, output,session) {

  
  invalidate <- reactiveValues(nextGA=1,future=1)
  observe({
    copyTo <- input$draggedfile[1]
    copy <- input$draggedfile[2]

    if(!is.null(copyTo)){
      if(copyTo %in% list.dirs(paste0(getwd(),"/models/"),recursive = F,full.names = F)){
    file.copy(copy, paste0("models/",copyTo), recursive=TRUE)
    
    mod <-read.csv(paste0(copy,"/mod.csv"),as.is=T)
    print(mod)
    i <- mod$Number[1]
    a <-mutate(mod,Path=paste0("models/",copyTo,"/mod",i,"/mod.ctl"))
    print(a)
    write.csv(a,paste0("models/",copyTo,"/mod",i,"/mod.csv"),row.names = F)
      }
  }
  })
  
  # Define globals
  alltokens <-data.frame(tokengroup= character(0), tokenset= character(0), token = character(0), stringsAsFactors=FALSE)
  GAProgress <- data.frame (Generation = numeric(0), Fitness = numeric(0))
  f <- NULL
  
  
  onclick("loadproj",{

    setwd(choose.dir("M:\\Users\\mhismail-shared\\Rprogramming\\genetic algorithm\\"))
    updateTextAreaInput(session,inputId = "ace",value=paste(readLines(list.files(pattern = ".*.ctl")),collapse = "\n")) 
    updateTextAreaInput(session,inputId = "preview",value=paste(readLines(list.files(pattern = ".*.ctl")),collapse = "\n")) 
    if (file.exists("alltokens.csv")){
      alltokens <- read.csv("alltokens.csv",as.is=T)
    }else({alltokens <-data.frame(tokengroup= character(0), tokenset= character(0), token = character(0), stringsAsFactors=FALSE)})
    })

# Saving and loading projects ---------------------------------------------

  # Save project on click ---------------------------------------------------
  onclick("saveproj",{
    write.csv(alltokens,"alltokens.csv",row.names=F)
    writeLines(input$ace,list.files(pattern = ".*.ctl")[1])
  })
  
  

  # Hide compartment numbers if Compartments isnt selected ------------------
  observe({
    if(input$strtypes!="Compartments"){
      hide("selectcmt")
      hide("ncmt")
    }else(show("ncmt"))
  })
# END ---------------------------------------------------------
  
  


# Parsing control stream --------------------------------------------------

  observe({
    x <- input$ace
    
    # Parse for tokens inside {} and add to tokengroupinput UI
    # choices is tokens found in control stream
    # choices2 is tokens found nested in other tokens
    choices<- gsub("\\{|:(.*)|\\}", "", regmatches(x, gregexpr("\\{.*?\\}", x))[[1]])
    choices2 <- gsub("\\{|:(.*)|\\}", "", regmatches(toString(alltokens$token), gregexpr("\\{.*?\\}", toString(alltokens$token)))[[1]])
    choices <- unique(c(choices,choices2))
    if(length(choices)==0){choices=c("")}
    updateAwesomeRadio(session,inputId = "tokengroupinput",choices =choices,selected = choices[1])
    
    # Parse control stream, search for "$INPUT", set covariate choices --------
    # Choices will be words following $INPUT, filtering out NONMEM reserved ---
    a<- strsplit(x,"\n")[[1]]
    if(length(strsplit(a[startsWith(a,"$INPUT")]," "))>0){
    covlist <- strsplit(a[startsWith(a,"$INPUT")]," ")[[1]][-1]
    covlist <- covlist[!(covlist%in%c("DATE=DROP","DATE","RATE","CMT","ID","TIME","AMT",
                                      "SS","ADDL","II","DV","MDV","EVID","DUR"))]
    updateSelectizeInput(session,"cov",choices = covlist)
    }
    
    # Parse control stream for compartments present
    b <- gsub(".*\\$MODEL(.*?)\\$.*","\\1",x)
    c<- regmatches(b, gregexpr("COMP.*\\(.*?\\)", b))[[1]]
    d <- gsub(".*\\((.*)\\).*","\\1",c)
    e<- gsub("\\,","",d)
    updateSelectizeInput(session,"selectcmt",choices = e)
    
    

  },priority = 100)
# END ---------------------------------------------------------------------
  
 
# Parse control stream to see if required tokens are present and --------
#       change UI accordingly by adding/removie css class ---------------
  
  observe({
    x <- input$ace
    
    #check if required tokens are present
    selectedtokengroup <- paste0("\\{",input$tokengroupinput,"\\}")
    selectedtokengroupnumbered <-  paste0("\\{",input$tokengroupinput,":[0-9]\\}")
    
    numberedtokengroups <- as.numeric(gsub(paste0("\\{",selectedtokengroup,":|\\}"), 
                                           "",
                                           regmatches(x, gregexpr(paste0("\\{",selectedtokengroup,":([0-9])\\}"), x))[[1]]))
      
    for (k in numberedtokengroups){
      x<- gsub(paste0("\\{",selectedtokengroup,":",k,"\\}"),tokenlist[k],x)
    }
    
    PKblock <- gsub(".*\\$PK(.*?)\\$.*","\\1",x)
    THETAblock <- gsub(".*\\$THETA(.*?)\\$.*","\\1",x)
    
    if(grepl(selectedtokengroup,PKblock)|grepl(selectedtokengroupnumbered,PKblock)){
      removeClass("covrequiredpk","red")
      addClass("covrequiredpk","green")
    }else{
      removeClass("covrequiredpk","green")
      addClass("covrequiredpk","red")}
    
    if(grepl(selectedtokengroup,THETAblock)|grepl(selectedtokengroupnumbered,THETAblock)){
      removeClass("covrequiredtheta","red")
      addClass("covrequiredtheta","green")
    }else{
      removeClass("covrequiredtheta","green")
      addClass("covrequiredtheta","red")}
    
#     fun<-paste0("$('#ace').highlightTextarea({words: ['{', '}'],
#                                                           id: 'highlighted2'});")
#     
#     runjs(fun)
    
  },priority = 100)
  
# END ---------------------------------------------------------------------
  
  
  
# filter tokenset list for selected token group ---------------------------
  observe({
    if(!is.null(input$tokengroupinput)){
    selectedTokenGroup <-input$tokengroupinput
    tokensetlist <- filter(alltokens,tokengroup==selectedTokenGroup)$tokenset%>%as.character()

    hide("tokeneditdiv",anim=T) #close editing tab on selection of new tokengroup
    updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
        }
     },priority = -10)
# END ---------------------------------------------------------------------
  
  
# check if cov is in tokengroup and set to selected value if match is found ----------------------
  observe({
    selectedTokenGroup <- input$tokengroupinput
    x <-input$ace
    a<- strsplit(x,"\n")[[1]]
    if(length(strsplit(a[startsWith(a,"$INPUT")]," "))>0){
      covlist <- strsplit(a[startsWith(a,"$INPUT")]," ")[[1]][-1]
      covlist <- covlist[!(covlist%in%c("DATE=DROP","DATE","RATE","CMT","ID","TIME","AMT",
                                        "SS","ADDL","II","DV","MDV","EVID","DUR"))]
      bestguess <- covlist[str_detect(selectedTokenGroup,covlist)][which(nchar(covlist[str_detect(selectedTokenGroup,covlist)])==max(nchar(covlist[str_detect(selectedTokenGroup,covlist)])))]
      updateSelectizeInput(session,"cov",selected = bestguess)
    }
  })
# END ---------------------------------------------------------------------
  
  
# filter token list for selected token set ---------------------------
  observe({
    if(!is.null(input$tokengroupinput)&!is.null(input$tokensetinput)){
    selectedTokenGroup <- input$tokengroupinput
    selectedTokenSet <- input$tokensetinput
    tokenlist <- filter(alltokens,tokengroup==selectedTokenGroup,tokenset==selectedTokenSet)$token%>%as.character()
    updateAwesomeRadio(session,inputId = "tokeninput",choices =tokenlist,selected =unique(tokenlist)[1])
      }
    },priority = -11)
# END ---------------------------------------------------------------------
  
# highlight selected tokengroup -------------------------------------------
# highlight package used: https://github.com/garysieling/jquery-highlighttextarea
  onclick("tokengroupinput",
          {
          updateTabsetPanel(session,"tabs",selected = "Control Stream")

          runjs(" $('#ace').highlightTextarea('destroy')")
          fun<-paste0("$('#ace').highlightTextarea({
                                                      words: [{
                                                        color: '#FFFF00',
                                                        words: ['{",input$tokengroupinput,"}',
                                                                '{",input$tokengroupinput,":[0-9]}']
                                                      }
                                                      , {
                                                        color: '#ff1b0f',
                                                        words: []
                                                      }]
                                                      });")
          
#           words:[{
#             words: ['{",input$tokengroupinput,"}',
#                     '{",input$tokengroupinput,":[0-9]}'],
#             color: '#FFFF00'
#           },{words: ['{','}'],
#             color: '#FFFF00'}])
          
#           words: [{
#             color: '#ADF0FF',
#             words: ['Lorem ipsum', 'vulputate']
#           }, {
#             color: '#FFFF00',
#             words: ['Donec']
#           }]


          delay(1,runjs(fun)) #needs delay for tokengroup change to register first
          })
  
# END ---------------------------------------------------------------------
  
  
# show how token set appears in control stream ----------------------------
  onclick("preview",
          {
            x <- input$ace
            selectedTokenGroup <- input$tokengroupinput
            selectedTokenSet <- input$tokensetinput
            tokengrouptext <- paste0("\\{",input$tokengroupinput,"\\}")
            tokenlist <- filter(alltokens,tokengroup==selectedTokenGroup,tokenset==selectedTokenSet)$token%>%as.character()
            tokenlist[tokenlist=="N/A"]<-""
            for (i in 1:length(tokenlist)){
              x<-sub(tokengrouptext,tokenlist[i],x)
            }
            x<- gsub("\\{[^//}]*\\}","",x)
            updateTextAreaInput(session,"previewtext",value=x)
            words <-paste0("'",tokenlist,"',",collapse="")
            words<-substr(words, 1, nchar(words)-1)

            updateTabsetPanel(session,"tabs",selected = "Preview")
            })
# END ---------------------------------------------------------------------
  
  
  
# display selected token type div (ie covariate, ETA,EPS,etc.) ------------
  onclick("tokentype",
          {
           hide(selector =  ".tokentype")
           showElement(id=input$tokentype)
          })
# END ---------------------------------------------------------------------
  
  
  
# Delete functions -------------------------------------------------------
  onclick("deletetokengroup",
          {
            x <- input$ace
            x<- gsub(paste0("\\{",input$tokengroupinput,"\\}"),"",x)
            alltokens<- filter(alltokens,tokengroup!=input$tokengroupinput)
            updateTextAreaInput(session,"ace",value=x)
          })
  
  onclick("deletetokenset",
          {
            alltokens<- filter(alltokens,!(tokengroup==input$tokengroupinput & tokenset==input$tokensetinput))
            tokensetlist <- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
            updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
            
          })
# END ---------------------------------------------------------------------
  
  

# Add to existing tokengroup ----------------------------------------------

onclick("addtokensetedit",{

  
  tokenToAdd<-data.frame(tokengroup=input$tokengroupinput,tokenset="New Token",token=c("N/A","N/A"))
  
  # If tokenset exists already, replace it
  alltokens <- filter(alltokens,!(tokengroup%in%tokenToAdd$tokengroup & tokenset%in%tokenToAdd$tokenset))
  alltokens<-rbind(alltokens,tokenToAdd)
  
  tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
  updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
  
  
  tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
  updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)

  selectedtokengroup <- input$tokengroupinput
  removeUI(selector=".tokeneditdiv",immediate=T) #.tokeneditdiv class is tokenset panel div
  
  #only display panel if tokensets exist for selected tokengroup
  if (length(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])>0){
    insertUI(selector="#tokenedit",where="beforeEnd",
             ui = {
               column(12, align="center",
                      textInput("edittokengroupname","Token Group Name",value = selectedtokengroup,width = "33%"),
                      lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])), function(i) {
                        column(12,
                               textInput(paste0("tokenset",i),"Token Set Name",value =unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i] ,width = "33%"),
                               lapply(1:length(alltokens$token[alltokens$tokengroup==input$tokengroupinput & alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]]),function(j){
                                 column(6,
                                        textAreaInput(paste0("edittoken",i,j),paste("Token",j),value = alltokens$token[alltokens$tokengroup==input$tokengroupinput & 
                                                                                                                         alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]][j],
                                                      resize="vertical"))}),
                               class="tokeneditdiv")}),
                      column(12,
                             actionButton("addtokensetedit","Add Token Set"),
                             actionButton("addtokenedit","Add Token"),
                             actionButton("savetokenedit",label="Save")),class="tokeneditdiv", id="tokeneditdiv")
             }, session = session, immediate=T)
  }
  })
  
# END---------------------------------------------------------------------

# Add token to existing tokenset ----------------------------------------------

onclick("addtokenedit",{
  
  tokensets <- unique(filter(alltokens,tokengroup==input$tokengroupinput)$tokenset)
  tokenToAdd<-data.frame(tokengroup=input$tokengroupinput,tokenset=tokensets,token="N/A")
  
  alltokens<-rbind(alltokens,tokenToAdd)
  
  tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
  updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
  
  
  tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
  updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)
  
  selectedtokengroup <- input$tokengroupinput
  removeUI(selector=".tokeneditdiv",immediate=T) #.tokeneditdiv class is tokenset panel div
  
  #only display panel if tokensets exist for selected tokengroup
  if (length(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])>0){
    insertUI(selector="#tokenedit",where="beforeEnd",
             ui = {
               column(12, align="center",
                      textInput("edittokengroupname","Token Group Name",value = selectedtokengroup,width = "33%"),
                      lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])), function(i) {
                        column(12,
                               textInput(paste0("tokenset",i),"Token Set Name",value =unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i] ,width = "33%"),
                               lapply(1:length(alltokens$token[alltokens$tokengroup==input$tokengroupinput & alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]]),function(j){
                                 column(6,
                                        textAreaInput(paste0("edittoken",i,j),paste("Token",j),value = alltokens$token[alltokens$tokengroup==input$tokengroupinput & 
                                                                                                                         alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]][j],
                                                      resize="vertical"))}),
                               class="tokeneditdiv")}),
                      column(12,
                             actionButton("addtokensetedit","Add Token Set"),
                             actionButton("addtokenedit","Add Token"),
                             actionButton("savetokenedit",label="Save")),class="tokeneditdiv", id="tokeneditdiv")
             }, session = session, immediate=T)
  }
})

# END---------------------------------------------------------------------
  
  
# Create and add covariate tokens  ----------------------------------------
  
  onclick("addcovs",
          {tokenTypes<- data.frame(tokengroup = input$tokengroupinput,tokenset = input$covtypes, stringsAsFactors=FALSE)
          x <- input$ace
          a<- strsplit(x,"\n")[[1]]
          covlist <- strsplit(a[startsWith(a,"$INPUT")]," ")[[1]][-1]
          
          if(any(covlist==input$cov)){
          values<-as.numeric(datafile()[which(covlist==input$cov)][[1]])
          mean <- median(values,na.rm=T)%>%round(2)
  
          if (input$center==T){
            
            lookup <- data.frame(tokenset=rep(c("None","Linear","Power","Exponential","Proportional"),each=2),
                                 token= c("N/A",
                                          "N/A",
                                          paste0("+(",input$cov,"-",mean,")*THETA([1])"),
                                          "(-10,.1,10)",
                                          paste0("*(",input$cov,"/",mean,")**THETA([1])"),
                                          "(-5,.5,10)",
                                          paste0("*exp((",input$cov,"-",mean,")*THETA([1]))"),
                                          "(-50,.1,50)",
                                          paste0("*(1+(",input$cov,"-",mean,")*THETA([1]))"),
                                          "(-100,1,100)"), stringsAsFactors=FALSE)
            
          }
          }
          
          if (input$center==F){
            
            lookup <- data.frame(tokenset=rep(c("None","Linear","Power","Exponential","Proportional"),each=2),
                                 token= c("N/A",
                                          "N/A",
                                          paste0("+(",input$cov,")*THETA([1])"),
                                          "(-100,.01,100)",
                                          paste0("*",input$cov,"**THETA([1])"),
                                          "(-100,1,100)",
                                          paste0("*exp(",input$cov,")*THETA([1])"),
                                          "(-100,.01,100)",
                                          paste0("*(1+",input$cov,"*THETA([1]))"),
                                          "(-100,1,100)"), stringsAsFactors=FALSE)
            
          }





          tokenToAdd<-merge(tokenTypes,lookup,by="tokenset")%>%select(tokengroup,tokenset,token)
          
          # If tokenset exists already, replace it
          alltokens <- filter(alltokens,!(tokengroup%in%tokenToAdd$tokengroup & tokenset%in%tokenToAdd$tokenset))
          alltokens<-rbind(alltokens,tokenToAdd)
          
          tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
          updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
          

          tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
          updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)
          hide("Covariate",anim=T)
          }
          )
  
# END ---------------------------------------------------------------------
  
  
# Create and add eta tokens  ----------------------------------------
  
  onclick("addetas",
          {
          tokenTypes<- data.frame(tokengroup = input$tokengroupinput,tokenset = input$etatypes)
          lookup <- data.frame(tokenset=rep(c("None","Normal","Logarithmic","Normal (Proportional)"),each=2),
                               token= c("N/A",
                                        "N/A",
                                        paste0("+ETA([1])"),
                                        "(0.01)",
                                        paste0("*exp(ETA([1]))"),
                                        "(0.01)",
                                        paste0("*(1+ETA([1]))"),
                                        "(0.01)"), stringsAsFactors=FALSE)
          tokenToAdd<-merge(tokenTypes,lookup,by="tokenset")%>%select(tokengroup,tokenset,token)
          alltokens <- filter(alltokens,!(tokengroup%in%tokenToAdd$tokengroup & tokenset%in%tokenToAdd$tokenset))
          
          alltokens<-rbind(alltokens,tokenToAdd)
          
          tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
          updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
          
          
          tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
          updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)
          hide("ETA",anim=T)
          }
  )
# END ---------------------------------------------------------------------
  
  
# Create and add eps tokens  ----------------------------------------
  
  onclick("addeps",
          {
            tokenTypes<- data.frame(tokengroup = input$tokengroupinput,tokenset = input$epstypes)
            lookup <- data.frame(tokenset=rep(c("Additive","Proportional","Additive + Proportional"),each=2),
                                 token= c(paste0("Y=F+EPS([1])"),
                                          "(0.01)",
                                          paste0("Y=F*(1+EPS([1]))"),
                                          "(0.01)",
                                          paste0("Y=F*(1+EPS([1]))+EPS([2])"),
                                          "(0.01) \n (0.01)"), stringsAsFactors=FALSE)
            tokenToAdd<-merge(tokenTypes,lookup,by="tokenset")%>%select(tokengroup,tokenset,token)
            alltokens <- filter(alltokens,!(tokengroup%in%tokenToAdd$tokengroup & tokenset%in%tokenToAdd$tokenset))
            
            alltokens<-rbind(alltokens,tokenToAdd)
            
            tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
            updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
            
            
            tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
            updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)
            
            hide("EPS",anim=T)
          }
  )
# END ---------------------------------------------------------------------
  
  
# Create and add structure tokens  ----------------------------------------
  onclick("addstr",
          {
            
            if(input$strtypes=="Compartments"){

              b <- gsub(".*\\$MODEL(.*?)\\$.*","\\1",input$ace)
              c<- regmatches(b, gregexpr("COMP.*\\(.*?\\)", b))[[1]]
              d <- gsub(".*\\((.*)\\).*","\\1",c)
              e<- gsub("\\,","",d)
              
              cmtn <- which(e==input$selectcmt)
              ncmts <- length(e)
              
            tokenTypes<- data.frame(tokengroup = input$tokengroupinput,tokenset = c("None",paste0("+",1:input$ncmt,"cmt")))
            lookup <- data.frame(tnumber=rep(c(1,2,3),4),tokenset=rep(c("None","+1cmt","+2cmt","+3cmt"),each=3),
                                 token= c("N/A","N/A","N/A",
                                          as.character(p1cmttext(input$selectcmt,cmtn,ncmts)[1]),
                                          as.character(p1cmttext(input$selectcmt,cmtn,ncmts)[2]),
                                          "(0 10)\n(0 10)",
                                          as.character(p2cmttext(input$selectcmt,cmtn,ncmts)[1]),
                                          as.character(p2cmttext(input$selectcmt,cmtn,ncmts)[2]),
                                          "(0 10)\n(0 10)\n(0 10)\n(0 10)",
                                          as.character(p3cmttext(input$selectcmt,cmtn,ncmts)[1]),
                                          as.character(p3cmttext(input$selectcmt,cmtn,ncmts)[2]),
                                          "(0 10)\n(0 10)\n(0 10)\n(0 10)\n(0 10)\n(0 10)"), stringsAsFactors=FALSE)
            tokenToAdd<-merge(lookup,tokenTypes,by="tokenset")%>%arrange(tokenset,tnumber)%>%select(tokengroup,tokenset,token)
            alltokens <- filter(alltokens,!(tokengroup%in%tokenToAdd$tokengroup & tokenset%in%tokenToAdd$tokenset))
            alltokens<-rbind(alltokens,tokenToAdd)
            
            tokensetlist<- filter(alltokens,tokengroup==input$tokengroupinput)$tokenset%>%as.character()
            updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
            
            
            tokenlist<- filter(alltokens,tokenset==unique(tokensetlist)[1])$token%>%as.character()
            updateAwesomeRadio(session,inputId = "tokeninput",choices =as.character(tokenlist),selected = NULL)
            
            hide("Structure",anim=T)
            }
          }
  )
# END ---------------------------------------------------------------------
  
  
# Show panel of tokensets when double click token group  ------------------
  
  onevent("dblclick","tokengroupinput",{
    
    selectedtokengroup <- input$tokengroupinput
    removeUI(selector=".tokeneditdiv",immediate=T) #.tokeneditdiv class is tokenset panel div
    
    #only display panel if tokensets exist for selected tokengroup
    if (length(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])>0){
    insertUI(selector="#tokenedit",where="beforeEnd",
             ui = {
               column(12, align="center",
               textInput("edittokengroupname","Token Group Name",value = selectedtokengroup,width = "33%"),
               lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])), function(i) {
                column(12,
               textInput(paste0("tokenset",i),"Token Set Name",value =unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i] ,width = "33%"),
               lapply(1:length(alltokens$token[alltokens$tokengroup==input$tokengroupinput & alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]]),function(j){
              column(6,
                      textAreaInput(paste0("edittoken",i,j),paste("Token",j),value = alltokens$token[alltokens$tokengroup==input$tokengroupinput & 
                                                                                                 alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]][j],
                                    resize="vertical"))}),
               class="tokeneditdiv")}),
               column(12,
               actionButton("addtokensetedit","Add Token Set"),
               actionButton("addtokenedit","Add Token"),
               actionButton("savetokenedit",label="Save")),class="tokeneditdiv", id="tokeneditdiv")
             }, session = session, immediate=T)
    }
  })
  
# END ---------------------------------------------------------------------
  
  
# Save tokens after editing  ----------------------------------------------
  
  onclick("savetokenedit",{
          alltokens$tokengroup[alltokens$tokengroup==input$tokengroupinput]<-input$edittokengroupname
          
          for (i in 1:length(unique(alltokens$tokenset[alltokens$tokengroup==input$edittokengroupname]))){
          alltokens$tokenset[alltokens$tokengroup==input$edittokengroupname &
                               alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$edittokengroupname])[i]]<- input[[paste0("tokenset",i)]]
          
          
          for (j in 1:length(alltokens$token[alltokens$tokengroup==input$edittokengroupname & alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$edittokengroupname])[i]])){
            alltokens$token[alltokens$tokengroup==input$edittokengroupname & 
                              alltokens$tokenset==unique(alltokens$tokenset[alltokens$tokengroup==input$edittokengroupname])[i]][j] <- input[[paste0("edittoken",i,j
                                                                                                                                                     )]]
          }
            
          }

          x <- input$ace
          tokengrouptext <- paste0("\\{",input$tokengroupinput,"\\}")

          x<-gsub(tokengrouptext, paste0("\\{",input$edittokengroupname,"\\}"),x)

          
          updateTextAreaInput(session,"ace",value=x)
          
          tokensetlist <- filter(alltokens,tokengroup==input$edittokengroupname)$tokenset%>%as.character()
          updateAwesomeRadio(session,inputId = "tokensetinput",choices =unique(tokensetlist),selected =unique(tokensetlist)[1])
          
          tokenlist <- filter(alltokens,tokengroup==input$edittokengroupname,tokenset==unique(tokensetlist)[1])$token%>%as.character()
          updateAwesomeRadio(session,inputId = "tokeninput",choices =unique(tokenlist),selected =unique(tokenlist)[1])
          
          
          hide("tokeneditdiv",anim=T)
          
  }
  )
  
# END ---------------------------------------------------------------------
  
  

  onclick("refreshmods",{
    
    #delete models folder if it exists
    unlink("models",force = TRUE,recursive = T)
    
    #create new models folder
    dir.create("models")
    dir.create("models/All")
    
    allmods <- allcombos(alltokens)
    write.csv(allmods,"allmods.csv")
    # write.csv(data.frame("Number","Path","OFV","Fitness","S","C","NTHETA","NETA",names(allmods)),"allmodsresults.csv")
    
    allmods[] <- lapply(allmods, as.character) 
    
    ctlstream <- input$ace
    #copy data to subfolder
    a<- strsplit(ctlstream,"\n")[[1]]
    datapath <- strsplit(a[startsWith(a,"$DATA")]," ")[[1]][2]
    file.copy(datapath,paste0("models"))
    
    
    alltokens$token <- as.character(alltokens$token) 
    alltokens$tokenset <- as.character(alltokens$tokenset) 
    
    ctlstream <- gsub("\\$DATA\\s*", "\\$DATA ../../",ctlstream) 
    create.model<-create.model 
    
    cl<-makeCluster(20)  
    registerDoSNOW(cl) 
    
    # for each row, for each column replace {tokengroup} with the cretaed token
    a <- foreach ( i = 1:dim(allmods)[1],.packages = c("dplyr","gsubfn"),.combine = "rbind") %dopar% {
      x <- ctlstream
      
      x<-create.model(x,alltokens,allmods[i,],names(allmods))
      
      dir.create(paste0("models/All/mod",i))
      writeLines(x,paste0("models/All/mod",i,"/mod.ctl"))
      
      maxtheta <- gsub("[THETA\\(\\)]", "", regmatches(x, gregexpr("THETA\\(.*?\\)", x))[[1]])
      maxtheta<-max(as.numeric(c(maxtheta,0)),na.rm=T)
      
      maxeta <- gsub("[ETA\\(\\)]", "", regmatches(x, gregexpr("\\bETA\\(.*?\\)", x))[[1]])
      maxeta<-max(as.numeric(c(maxeta,0)),na.rm=T)
      a <-cbind(data.frame(Number=i,Path=paste0("models/All/mod",i,"/mod.ctl")),OFV="",Fitness="",S="",C="",NTHETA=maxtheta,NETA=maxeta,allmods[i,])

      
      write.csv(a,paste0("models/All/mod",i,"/mod.csv"),row.names = F)
      a
    }

    
    stopCluster(cl)
    write.csv(a,"allmodsresults.csv",row.names = F)
  })
  
  
  
  # Create table to be displayed in UI
  allmods2<- reactive({
    
    
    z<-Sys.time()
    # input$refreshmods
    input$refreshmods2
    allmodsresults <- read.csv("allmodsresults.csv")
    
    if (file.exists("modelruntemp.csv")&file.info("modelruntemp.csv")$size>0){
      modelsran <- read.csv("modelruntemp.csv",header = F)%>%distinct(1,.keep_all=T)
      print(modelsran)
      checkpresent <- file.exists(paste0(modelsran[,1],"/results/PsN_execute_plots.R"))
      resultspresent <- modelsran[checkpresent,]
      print(dim(resultspresent))
      if (dim(resultspresent)[1]>0){
      print(resultspresent)
      results <- lapply(resultspresent[,1],retrieveresults)
      removerefreshed <- modelsran[-checkpresent,]
      print(removerefreshed)
      write.table(removerefreshed,
                  "modelruntemp.csv", sep = ",",row.names = F,col.names = F)
      print(results)
      resultsdf <- data.frame()
      resultsdf <- rbind(resultsdf, do.call(rbind, results))
      resultsdf[,1]<-as.numeric(as.character(resultsdf[,1]))
      mods <- resultsdf[,1]
      print(match(mods,allmodsresults$Number))
      print(allmodsresults)
      match(allmodsresults$Number,mods)
      
      #logic for below - match function finds the indices of allmodsresults in which the number equals the 
      # value in mods. na.omit removes na values. 1,3,5,6 are columns to be updated
      # right side !is.na is safegaurd for models that were run but have since been deleted before results were fetched
      allmodsresults[na.omit(match(mods,allmodsresults$Number)),c(1,3,5,6)] <- resultsdf[!is.na(match(mods,allmodsresults$Number)),]
      write.csv(allmodsresults,"allmodsresults.csv",row.names = F)
      }
      
    }

    
    updateRadioGroupButtons(session,"selecteddir",selected=selecteddir,choices =list.dirs("./models",full.names = F,recursive = F))
    delay(100,runjs('$(".radiobtn").on("dragenter",function () {
        y=$(this).find("input").val()
          console.log("hi")
      })'))
    allmodsresults
    
    
    
    
  })
  
 

  # call PsN execute to run model
  onclick("runmod",{
    
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    for (i in 1:length(input$modeltable_selected)){
    path<-as.character(allmods2()[input$modeltable_selected[i]+1,2])
    relpath <- sub("/mod.ctl","",path)
    unlink(paste0(homepath,'/',relpath,"/results"),recursive = TRUE)
    run.model(paste0(homepath,'/',relpath),basename(relpath))
    }
    
  })
  
  # Delete selected model/s
  onclick("deletemods",{
    
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    for (i in 1:length(input$modeltable_selected)){
      path<-as.character(allmods2()[input$modeltable_selected[i]+1,2])
      write.csv(allmods2()[-(input$modeltable_selected[i]+1),],"allmodsresults.csv",row.names=F)
      relpath <- sub("/mod.ctl","",path)
      unlink(paste0(homepath,'/',relpath),recursive = TRUE)
    }
    
      runjs('$("#refreshmods2").click()')
  })
  
  # open file explorer
  onclick("openlocation",{
  
      path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
      relpath <- sub("/mod.ctl","",path)
      
      shell(gsub("/","\\\\",paste("explorer", relpath)))
    
  })
  
  onclick("addsimilar",{
    
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
    
    phen <- read.csv(paste0(relpath,"/mod.csv"),as.is=T)[-(1:8)]
    allmods <- read.csv("allmods.csv",as.is=T)[-1]
    
    similarmods <- which(apply(allmods, 1, function(x) sum( x ==phen)>=(dim(allmods)[2]-1)))
    
    copyTo <- input$selecteddir
    copy <- paste0("models/All/mod",similarmods)
    for (i in copy){
    if(!is.null(copyTo)){
      print(copyTo)
      print(i)
      file.copy(i, paste0("models/",copyTo), recursive=TRUE)
      
      mod <-read.csv(paste0(i,"/mod.csv"),as.is=T)
      print(mod)
      i <- mod$Number[1]
      a <-mutate(mod,Path=paste0("models/",copyTo,"/mod",i,"/mod.ctl"))
      print(a)
      write.csv(a,paste0("models/",copyTo,"/mod",i,"/mod.csv"),row.names = F)
    }
    }
    
    runjs('$("#refreshmods2").click()')
    
    
  })
  
  # All models modal --------------------------------------------------------
  observeEvent(input$viewmod, {
    showModal(
      table
    )
    })
  
  # Render table
  output$modeltable <- DT::renderDataTable({
    allmods2()
  },
  # filter = list(position = "top"),
  options = list(
    select = list(style = 'os', # set 'os' select style so that ctrl/shift + click in enabled
                  items = 'row'), # items can be cell, row or column
    pageLength = 100,
    dom = 'tpf', 
    columnDefs = list(list(className = 'dt-center', targets = "_all"))),
  rownames= FALSE, 
  extensions = 'Select', 
  selection="none", 
  callback = JS(
    "table.on( 'click.dt', 'tbody td', function (e) {", # react on click event
    "var type = table.select.items();", # get the items setting of Select extension
    "var idx = table[type + 's']({selected: true}).indexes().toArray();", # get the index of selected items
    "var DT_id = table.table().container().parentNode.id;", # get the output id of DT
    "Shiny.onInputChange('modeltable' + '_selected', idx);", # send the index to input$outputid_selected
    "})"
  ), server = F)
  
  table <-  modalDialog(fluidPage(div(column(9,
                                  actionButton("runmod","Run Selected Model"),
                                  actionButton("refreshmods","Recreate Models"),
                                  actionButton("refreshmods2","Refresh"),
                                  actionButton("deletemods","Delete Selected Models"),
                                  actionButton("editmod","Edit Selected Model"),
                                  actionButton("openlocation","Open Model Location"),
                                  actionButton("addsimilar","Add Similar Models"),
                                  DT::dataTableOutput('modeltable'),

                                  radioGroupButtons("selecteddir",NULL,choices = c("All"),selected = "All"),
                                  actionButton("showadddir",icon("plus","fa-1x")),
                                  span(id= "add-directory", style="display: none;",div(style="display: inline-block;vertical-align:top; width: 150px;",textInput("boog",NULL,"boog",width="auto")),actionButton("adddir","Add")),
                                  
                                  ondrop="drop(event)", ondragover="allowDrop(event)")),
                                  div(column(3, tabsetPanel(id = "modelinfo",
                                                        tabPanel("Plots",
                                                                 bsCollapse(bsCollapsePanel("Covariate Model",
                                                                 actionButton("ranparvscov","ETAs vs Covariates"),
                                                                 actionButton("covscatter","Covariate Scatter Plots")),
                                                                 bsCollapsePanel("Goodness of Fit",
                                                                                 actionButton("basicgof","Basic GOF"),
                                                                                 actionButton("indplots","Individual Plots"),
                                                                                 actionButton("dvpredipred","DV vs Pred, IPRED")),
                                                                 multiple = T,open=c("Covariate Model"))
                                                                 ),
                                                        tabPanel("Parameters",
                                                                 bsCollapse(bsCollapsePanel("Structural (THETA)",
                                                                                            tableOutput("thetas")
                                                                                            ),
                                                                            bsCollapsePanel("IIV (OMEGA)",
                                                                                            tableOutput("omegas")),
                                                                            bsCollapsePanel("RUV (SIGMA)",
                                                                                            tableOutput("sigmas")),
                                                                            multiple = T,
                                                                            open=c("Structural (THETA)","IIV (OMEGA)","RUV (SIGMA)"))),
                                                        tabPanel("Covariate Model",
                                                                 tableOutput("regress"),
                                                                 tableOutput("regress.parms")),
                                                        tabPanel("PsN",
                                                                 actionButton("runvpc","VPC")))))),easyClose = T)
  
  onclick("showadddir",{
    toggle("add-directory",anim = T,animType = "fade")
  })
  
  onclick("adddir",{
    dir.create(paste0(getwd(),"/models/",input$boog))
    hide("add-directory",anim = T,animType = "fade")
    
    updateRadioGroupButtons(session,"selecteddir",selected=input$selecteddir,choices =list.dirs("./models",full.names = F,recursive = F))
    delay(100,runjs('$(".radiobtn").on("dragenter",function () {
        y=$(this).find("input").val()
          console.log("hi")
      })'))

    
  })
  
  
  
  
  observeEvent(input$startGA,{invalidate$future=isolate(invalidate$future)+1})
  observeEvent(input$nextGA,{invalidate$future=isolate(invalidate$future)+1})
  
  
  observeEvent({
    invalidate$future
    },{
    print("started")
      delay(10,{
    f <<- future({
      b<-1
      while(length(b)>0){
        a <- shell("tasklist /v",intern=T)
        b<- as.vector(na.omit(str_extract(a,"mod[0-9]* - execute")))
      }}) %plan% multiprocess
    invalidate$nextGA <- isolate(invalidate$nextGA)+1})
  },ignoreInit = T)
#   
#   
  observe({
    invalidate$nextGA
  
    if(!is.null(f)){
      if (!resolved(f)){
        print("not done")
        invalidateLater(10000, session)
      }
    
      if(resolved(f)){
        print("boog")
        runjs('$("#refreshGAmods").click()')
        runjs('$("#nextGA").click()')
      }
    }
    
  })
  
  

  
# Genetic Algorithm Modal -------------------------------------------------
  
  GAmods<- reactive({
    input$refreshGAmods
    
    
    allmods1 <- data.frame ()
    allmodslist <-list()
    mod.directories <- read.csv("GAgens.csv",as.is=T)%>%
      filter(Generation == as.numeric(input$selectedgen))%>%select(Dir)
  

    allmodslist<-lapply(paste0("./",mod.directories$Dir),retrieveresults)
    

    allmods1<- rbind(allmods1, do.call(rbind, allmodslist))
    allmods1<-mutate(allmods1,Number = as.numeric(Number))%>%arrange(Number)

    
    allmods1
    
    
  })
  
  GAmodsAllGens<- reactive({
    input$refreshGAmods
    
    
    allmods1 <- data.frame ()
    allmodslist <-list()
    mod.directories <- read.csv("GAgens.csv",as.is=T)
    mod.directories.max.gen <- filter(mod.directories,Generation==max(Generation))
    
    
    allmodslist<-lapply(paste0("./",mod.directories.max.gen$Dir),retrieveresults)
    
    
    allmods1<- rbind(allmods1, do.call(rbind, allmodslist))
    allmods1<-mutate(allmods1,Number = as.numeric(Number),Generation = mod.directories.max.gen$Generation)%>%
      arrange(Number)%>%
      select(Generation,Fitness)
    GAProgress <- rbind(GAProgress,allmods1)
    print(GAProgress)


  })
  
  output$GAtable <- DT::renderDataTable({
    GAmods()
  },#filter = 'top',
  options = list(
    select = list(style = 'os', # set 'os' select style so that ctrl/shift + click in enabled
                  items = 'row'), # items can be cell, row or column
    pageLength = 100,
    dom = 't', 
    columnDefs = list(list(className = 'dt-center', targets = "_all"))),
  rownames= FALSE, 
  extensions = 'Select', 
  selection="none", 
  callback = JS(
    "table.on( 'click.dt', 'tbody td', function (e) {", # react on click event
    "var type = table.select.items();", # get the items setting of Select extension
    "var idx = table[type + 's']({selected: true}).indexes().toArray();", # get the index of selected items
    "var DT_id = table.table().container().parentNode.id;", # get the output id of DT
    "Shiny.onInputChange('modeltable' + '_selected', idx);", # send the index to input$outputid_selected
    "})"
  ), server = F)
  
  output$GAPlot <- renderPlot({
    GAProgress <<- GAmodsAllGens()
    GAmodsAllGensGrouped <- group_by(GAProgress,Generation)
    summary <- summarize(GAmodsAllGensGrouped, meanfit = mean(as.numeric(Fitness),na.rm=T),minfit = min(Fitness,na.rm=T))
    ggplot(summary,aes(x=Generation,y=minfit),color="blue")+geom_line() +geom_line(aes(y=meanfit),color="red") 
  })
  
  observeEvent(input$openGA, {
    showModal(
      GAmodal
    )
  })
  
  onclick("startGA",{
    GAProgress <- data.frame (Generation = numeric(0), Fitness = numeric(0))
    
    nmods <- length(list.files("./models/All"))
    nindiv = 40
    if (nindiv>nmods) nindiv<- nmods
    initiateGA(nmods,nindiv)
  })
  
  onclick("nextGA",{
    GAgens <- read.csv("GAgens.csv")
    mods <- filter(GAgens, Generation==max(Generation))
    mods <- paste0("./",mods$Dir)
    maxgen<-nextGA(mods,alltokens)
    updateRadioGroupButtons(session,"selectedgen",choices = 1:maxgen)
  })
  
  GAmodal <-  modalDialog(fluidPage(div(column(9,
                                               actionButton("startGA","Start"),
                                               actionButton("nextGA","Next"),
                                               actionButton("refreshGAmods","Refresh"),
                                               actionButton("openlocation2","Open Model Location"),
                                             DT::dataTableOutput('GAtable'))),
                                  div(column(3, tabsetPanel(id = "modelinfo2",
                                                            tabPanel("Plots",
                                                                     bsCollapse(bsCollapsePanel("Covariate Model",
                                                                                                actionButton("ranparvscov2","ETAs vs Covariates"),
                                                                                                actionButton("covscatter2","Covariate Scatter Plots")),
                                                                                bsCollapsePanel("Goodness of Fit",
                                                                                                actionButton("basicgof2","Basic GOF"),
                                                                                                actionButton("indplots2","Individual Plots"),
                                                                                                actionButton("dvpredipred2","DV vs Pred, IPRED")),
                                                                                multiple = T,open=c("Covariate Model")),
                                                                     plotOutput("GAPlot",height = "200px")
                                                            ),
                                                            tabPanel("Parameters",
                                                                     bsCollapse(bsCollapsePanel("Structural (THETA)",
                                                                                                tableOutput("thetas2")
                                                                     ),
                                                                     bsCollapsePanel("IIV (OMEGA)",
                                                                                     tableOutput("omegas2")),
                                                                     bsCollapsePanel("RUV (SIGMA)",
                                                                                     tableOutput("sigmas2")),
                                                                     multiple = T,
                                                                     open=c("Structural (THETA)","IIV (OMEGA)","RUV (SIGMA)"))),
                                                            tabPanel("PsN",
                                                                     actionButton("runvpc2","VPC"))))),
                                  column(12,radioGroupButtons("selectedgen",NULL,choices = c(1,2),selected = NULL))),easyClose = T)
  
  
  
  # Edit control stream modal -----------------------------------------------
  
  editcontrolstream <- modalDialog(fluidPage(actionButton("savecontrol","Save"),textAreaInput("editcontrol",
                                                           label = NULL,
                                                           width="100%",
                                                           height = "200px",
                                                           cols=NULL,
                                                           value="")))
  
  
  onclick("savecontrol",{
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    x<- input$editcontrol
    writeLines(x,paste0(path))
    })
  
  
  
  observeEvent(input$editmod, {
    
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA/"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    file <-paste0(homepath,"/",path)
    if (file.exists(file)){
      updateTextAreaInput(session,"editcontrol",value=paste(readLines(file),collapse = "\n"))
    }
    showModal(
      editcontrolstream
    )})
  
  
  #  ------------------------------------------------------------------------
  
  datafile <- reactive({
    x<- input$ace
    a<- strsplit(x,"\n")[[1]]
    datapath <- strsplit(a[startsWith(a,"$DATA")]," ")[[1]][2]
    
    read.csv(datapath)
  })
  
  output$data <- DT::renderDataTable({
    datafile()
  },filter = 'top',options = list(pageLength = 100,scrollY = "70vh",scrollX="100%",dom = 't', autoWidth = TRUE,columnDefs = list(list(className = 'dt-center', targets = "_all"))),rownames= FALSE)
  
 
  # Diagrams ----------------------------------------------------------------
  
    output$distributiondiagram <- renderGrViz({
      b <- gsub(".*\\$MODEL(.*)\\$.*","\\1",input$ace)
      c<- regmatches(b, gregexpr("COMP.*\\(.*?\\)", b))[[1]]
      d <- gsub(".*\\((.*)\\).*","\\1",c)
      e<- gsub("\\,","",d)
      
      cmtn <- which(e==input$selectcmt)
      ncmts <- length(e)
      if (input$ncmt==1) a<-p1cmt(substr(input$selectcmt,1,7),cmtn,ncmts)   
      if (input$ncmt==2) a<-p2cmt(substr(input$selectcmt,1,7),cmtn,ncmts) 
      if (input$ncmt==3) a<-p3cmt(substr(input$selectcmt,1,7),cmtn,ncmts)  
      a
      
      
    })

  
   # Parameter tables --------------------------------------------------------
  
  output$thetas <- renderTable({
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
      
    a<-readLines(paste0(homepath,"/",relpath,"/results/raw_results_structure"))[-1]
    b<-strsplit(a,"[=|,]")
    c<-b[which(lengths(b)==3)]
    d <- do.call(cbind,c)
    final<- as.data.frame(d,stringsAsFactors = F)
    names(final)<- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2])+1,as.numeric(final$theta[2])+as.numeric(final$theta[3]))
    omegaindices <-  seq(as.numeric(final$omega[2])+1,as.numeric(final$omega[2])+as.numeric(final$omega[3]))
    sigmaindices <-  seq(as.numeric(final$sigma[2])+1,as.numeric(final$sigma[2])+as.numeric(final$sigma[3]))
    
    thetaSE <-  seq(as.numeric(final$setheta[2])+1,as.numeric(final$setheta[2])+as.numeric(final$setheta[3]))
    omegaSE <-  seq(as.numeric(final$seomega[2])+1,as.numeric(final$seomega[2])+as.numeric(final$seomega[3]))
    sigmaSE <-  seq(as.numeric(final$sesigma[2])+1,as.numeric(final$sesigma[2])+as.numeric(final$sesigma[3]))
    
    etaSHR <-  seq(as.numeric(final$shrinkage_eta[2])+1,as.numeric(final$shrinkage_eta[2])+as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath,"/",relpath,"/results/raw_results_mod.csv"))
    
    THETA <- rawResults[thetaindices] 
    THETASE <- abs(rawResults[thetaSE]/rawResults[thetaindices]*100)
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] *100
    ETASHR <- rawResults[etaSHR]
    
    SIGMA <- rawResults[sigmaindices]
    SIGMARSE <- rawResults[sigmaSE]/rawResults[sigmaindices]*100
    
    THETATBL <- as.data.frame(cbind(t(THETA),t(THETASE)))
    names(THETATBL) <- c("Estimate","RSE (%)")
    
    OMEGATBL <- as.data.frame(cbind(t(OMEGA),t(OMEGARSE),t(ETASHR)))
    names(OMEGATBL) <- c("Estimate","RSE (%)", "Shrinkage (%)")
    
    SIGMATBL <- as.data.frame(cbind(t(SIGMA),t(SIGMARSE)))
    names(SIGMA) <- c("Estimate","RSE (%)")
    THETATBL
  },rownames=T)
  
  output$omegas <- renderTable({
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
    
    a<-readLines(paste0(homepath,"/",relpath,"/results/raw_results_structure"))[-1]
    b<-strsplit(a,"[=|,]")
    c<-b[which(lengths(b)==3)]
    d <- do.call(cbind,c)
    final<- as.data.frame(d,stringsAsFactors = F)
    names(final)<- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2])+1,as.numeric(final$theta[2])+as.numeric(final$theta[3]))
    omegaindices <-  seq(as.numeric(final$omega[2])+1,as.numeric(final$omega[2])+as.numeric(final$omega[3]))
    sigmaindices <-  seq(as.numeric(final$sigma[2])+1,as.numeric(final$sigma[2])+as.numeric(final$sigma[3]))
    
    thetaSE <-  seq(as.numeric(final$setheta[2])+1,as.numeric(final$setheta[2])+as.numeric(final$setheta[3]))
    omegaSE <-  seq(as.numeric(final$seomega[2])+1,as.numeric(final$seomega[2])+as.numeric(final$seomega[3]))
    sigmaSE <-  seq(as.numeric(final$sesigma[2])+1,as.numeric(final$sesigma[2])+as.numeric(final$sesigma[3]))
    
    etaSHR <-  seq(as.numeric(final$shrinkage_eta[2])+1,as.numeric(final$shrinkage_eta[2])+as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath,"/",relpath,"/results/raw_results_mod.csv"))
    
    THETA <- rawResults[thetaindices] 
    THETASE <- abs(rawResults[thetaSE]/rawResults[thetaindices]*100)
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] *100
    ETASHR <- rawResults[etaSHR]
    
    SIGMA <- rawResults[sigmaindices]
    SIGMARSE <- rawResults[sigmaSE]/rawResults[sigmaindices]*100
    
    THETATBL <- as.data.frame(cbind(t(THETA),t(THETASE)))
    names(THETATBL) <- c("Estimate","RSE (%)")
    
    OMEGATBL <- as.data.frame(cbind(t(OMEGA),t(OMEGARSE),t(ETASHR)))
    names(OMEGATBL) <- c("Estimate","RSE (%)", "Shrinkage (%)")
    
    SIGMATBL <- as.data.frame(cbind(t(SIGMA),t(SIGMARSE)))
    names(SIGMA) <- c("Estimate","RSE (%)")
    OMEGATBL
  },rownames=T)
  
  output$sigmas <- renderTable({
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
    
    a<-readLines(paste0(homepath,"/",relpath,"/results/raw_results_structure"))[-1]
    b<-strsplit(a,"[=|,]")
    c<-b[which(lengths(b)==3)]
    d <- do.call(cbind,c)
    final<- as.data.frame(d,stringsAsFactors = F)
    names(final)<- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2])+1,as.numeric(final$theta[2])+as.numeric(final$theta[3]))
    omegaindices <-  seq(as.numeric(final$omega[2])+1,as.numeric(final$omega[2])+as.numeric(final$omega[3]))
    sigmaindices <-  seq(as.numeric(final$sigma[2])+1,as.numeric(final$sigma[2])+as.numeric(final$sigma[3]))
    
    thetaSE <-  seq(as.numeric(final$setheta[2])+1,as.numeric(final$setheta[2])+as.numeric(final$setheta[3]))
    omegaSE <-  seq(as.numeric(final$seomega[2])+1,as.numeric(final$seomega[2])+as.numeric(final$seomega[3]))
    sigmaSE <-  seq(as.numeric(final$sesigma[2])+1,as.numeric(final$sesigma[2])+as.numeric(final$sesigma[3]))
    
    etaSHR <-  seq(as.numeric(final$shrinkage_eta[2])+1,as.numeric(final$shrinkage_eta[2])+as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath,"/",relpath,"/results/raw_results_mod.csv"))
    
    THETA <- rawResults[thetaindices] 
    THETASE <- abs(rawResults[thetaSE]/rawResults[thetaindices]*100)
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] *100
    ETASHR <- rawResults[etaSHR]
    
    SIGMA <- rawResults[sigmaindices]
    SIGMARSE <- rawResults[sigmaSE]/rawResults[sigmaindices]*100
    
    THETATBL <- as.data.frame(cbind(t(THETA),t(THETASE)))
    names(THETATBL) <- c("Estimate","RSE (%)")
    
    OMEGATBL <- as.data.frame(cbind(t(OMEGA),t(OMEGARSE),t(ETASHR)))
    names(OMEGATBL) <- c("Estimate","RSE (%)", "Shrinkage (%)")
    
    SIGMATBL <- as.data.frame(cbind(t(SIGMA),t(SIGMARSE)))
    names(SIGMATBL) <- c("Estimate","RSE (%)")
    SIGMATBL
  },rownames=T)
  
  output$regress <- renderTable({
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
    
    
    
    cotab<-read.table(paste0(relpath,"/cotab1"),skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)
    catab<-read.table(paste0(relpath,"/catab1"),skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)
    patab <- read.table(paste0(relpath,"/patab1"),skip = 1,header=T,stringsAsFactors = F)%>%distinct(ID,.keep_all = T)
    
    ETAS <- names(patab)[grepl("ETA",names(patab))]
    CO<- names(cotab)[!(grepl("ID",names(cotab)))]
    CA <- names(catab)[!(grepl("ID",names(catab)))]
    
    etatab <- patab[c(1,which(grepl("ETA",names(patab))))]
    
    alltab <- right_join(cotab,catab)%>%right_join(etatab)

    
    
    k = 0
    corlist <- list()
    for (i in ETAS){
      ETAi <- which(names(alltab)==i)
      for (j in CO){
        k=k+1
        PARAMi <- which(names(alltab)==j)
        cor <- lm (alltab[,ETAi]~alltab[,PARAMi])
        corlist[[k]]<-(c(R="Linear",
                         Model=paste0(j,"on",i),
                         Estimate =round(as.numeric(cor[[1]][2]),3), 
                         STDError = round(as.numeric(coef(summary(cor))[, 2][2] ),3)))
      }
    } 
    
    allcor<-data.frame()
    
    allcor<- rbind(allcor, do.call(rbind, corlist))
    

  },rownames=F)
  
  output$regress.parms <- renderTable({

    
    xpose.runno <- '1'
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
    relpath <- sub("/mod.ctl","",path)
    working.directory <- paste0(homepath,'/',relpath)
    setwd(working.directory)
    xpdb<-xpose.data(xpose.runno)
    setwd(homepath)
    
  

    
    k=0
    all <- list()
    for (i in xpdb@Prefs@Xvardef$parms){
      y <- select_(xpdb@Data,i)
      for (j in xpdb@Prefs@Xvardef$covariates){
        x <- select_(xpdb@Data,j)
        if(!is.factor(x[,1])){
          k=k+1
          a<-rbind(linmod(center(x[,1]),as.numeric(as.character(y[,1])),xlab=names(x),ylab=names(y)),
                   expmod(center(x[,1]),as.numeric(as.character(y[,1])),xlab=names(x),ylab=names(y)),
                   powmod(center2(x[,1]),as.numeric(as.character(y[,1])),xlab=names(x),ylab=names(y)))
          all[[k]] <- a
        }
      }
    }
    
    allcor<-data.frame()
    allcor<- rbind(allcor, do.call(rbind, all))
    
    
    
    
  },rownames=F)
  
  
  # Plot functions ----------------------------------------------------------
  
  #ETAs vs Covariates
  onclick("ranparvscov",
          {
            xpose.runno <- '1'
            homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
            relpath <- sub("/mod.ctl","",path)
            working.directory <- paste0(homepath,'/',relpath)
            setwd(working.directory)
            pdf.filename <- paste0('etavscov.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file=pdf.filename,width=10,height=7,title=pdf.title)
            print(ranpar.vs.cov(xpdb))
            dev.off()
            shell("etavscov.pdf",wait = F)
            setwd(homepath)
          })
  
  #Basic Goodness of fit - (DV vs POPPRED, DV vs IPRED, |IWRES| vs DV, CWRES vs Time)
  onclick("basicgof",
          {
            xpose.runno <- '1'
            homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
            relpath <- sub("/mod.ctl","",path)
            working.directory <- paste0(homepath,'/',relpath)
            setwd(working.directory)
            pdf.filename <- paste0('basicGOF.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file=pdf.filename,width=10,height=7,title=pdf.title)
            print(basic.gof(xpdb))
            dev.off()
            shell("basicGOF.pdf",wait = F)
            setwd(homepath)
          })
  
  onclick("indplots",
          {
            xpose.runno <- '1'
            homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
            relpath <- sub("/mod.ctl","",path)
            working.directory <- paste0(homepath,'/',relpath)
            setwd(working.directory)
            pdf.filename <- paste0('indplots.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file=pdf.filename,width=10,height=7,title=pdf.title)
            print(ind.plots(xpdb))
            dev.off()
            shell("indplots.pdf",wait = F)
            setwd(homepath)
          })
  
  onclick("dvpredipred",
          {
            xpose.runno <- '1'
            homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
            relpath <- sub("/mod.ctl","",path)
            working.directory <- paste0(homepath,'/',relpath)
            setwd(working.directory)
            pdf.filename <- paste0('dvpredipred.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file=pdf.filename,width=10,height=7,title=pdf.title)
            print(dv.vs.pred.ipred(xpdb))
            dev.off()
            shell("dvpredipred.pdf",wait = F)
            setwd(homepath)
          })
  
  onclick("covscatter",
          {
            xpose.runno <- '1'
            homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            path<-as.character(allmods2()[input$modeltable_selected[1]+1,2])
            relpath <- sub("/mod.ctl","",path)
            working.directory <- paste0(homepath,'/',relpath)
            setwd(working.directory)
            pdf.filename <- paste0('covscatter.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file=pdf.filename,width=10,height=7,title=pdf.title)
            print(cov.splom(xpdb))
            dev.off()
            shell("covscatter.pdf",wait = F)
            setwd(homepath)
          })
  
}








# Run the application 
shinyApp(ui = ui, server = server)


