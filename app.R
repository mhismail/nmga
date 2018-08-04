#STABLE-test
graphics.off()
rm(list = ls())

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
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("future")) install.packages("future")
if (!require("pkga")){
  if (!require("devtools")) install.packages("devtools")
  library(devtools)
  install("pkga")
  library(pkga)
}



#Set plot theme
theme_set(theme_bw())

# InitiateSCM <- function (path,control=NULL,alltokens=NULL,allmods=NULL){
#   results_dir <- "SCM"
#   i=1
#   if(dir.exists(results_dir)){
#     while(dir.exists(results_dir)){
#       results_dir<- paste0("SCM(",i,")")
#       if(suppressWarnings(dir.create(results_dir))){
#         break
#       }
#       i = i+1
#     }
#   }else{
#     suppressWarnings(dir.create(results_dir))
#   }
#   
#   copyTo<- paste0(results_dir,"/Base")
#   dir.create(copyTo)
#   CheckThenCreate(path,control,alltokens,allmods)
#   CopyModel(path,copyTo,control,alltokens,allmods) 
#   #   homepath <-getwd()
#   #   relpath <- path
#   # 
#   #   RunModel(paste0(homepath,'/',relpath),basename(relpath))
#   #      
#   #   
#   #   
#   #   write.csv(data.frame(Generation=1,Inidvidual=pop,Dir=popmods), "GAgens.csv",row.names=F)
# } 

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  div("NMGA", class = "title-bar"),
  column(12,
         span(icon("plus", "fa-2x"), id="newproj", class=c("nav-bar-el")),
         span(icon("save", "fa-2x"), id="saveproj", class=c("nav-bar-el")),
         span(icon("folder-open", "fa-2x"), id="dir", class=c("nav-bar-el")),
         span("Directory", id="proj", class=c("nav-bar-el")),
         span(actionButton("viewmod", "View Models", class="nav-bar-button")),
         span(actionButton("openGA", "Initiate Genetic Algorithm", class="nav-bar-button")),
         span(actionButton("GAsettings", icon("gear", "fa-2x"), class="nav-bar-button")),
         span(actionButton("openSCM", "Initiate SCM", class="nav-bar-button")),
         span(actionButton("SCMsettings", icon("gear", "fa-2x"), class="nav-bar-button")),
         div(id = "new_project_div", style="display: none;",
             span(h2("Control Stream", class="nav-bar-button"),
                  actionButton("select_ctl", "Select Control Stream", class="nav-bar-button"), class="nav-bar-el"),
             span(h2("Control Stream", class="nav-bar-button"),
                  actionButton("select_ctl", "Select Control Stream", class="nav-bar-button"), class="nav-bar-el")),
         class="nav-bar"),
  column(6,
         tabsetPanel(id="tabs", tabPanel("Control Stream",
                                         textAreaInput("ace",
                                                       label=NULL,
                                                       width="100%",
                                                       cols=NULL,
                                                       value="Select a directory with a single .ctl file")),
                     tabPanel("Preview",
                              textAreaInput("previewtext",
                                            label = NULL,
                                            width="100%",
                                            cols=NULL,
                                            value="")),
                     tabPanel("Data",
                              DT::dataTableOutput("data")))),
  column(6,
         column(12,
                align="center",
                radioGroupButtons("tokentype",
                                  NULL,
                                  choices = c("Covariate","ETA","EPS","Structure", "Custom"),
                                  selected = NULL)),
         column(12,
                column(4,
                       awesomeRadio(inputId = "tokengroupinput",
                                    "Token Group", 
                                    choices = c("")),
                       actionButton("deletetokengroup",
                                    "Delete",
                                    style="background-color: #b71616; color: white;")
                ),
                column(4,
                       awesomeRadio(inputId = "tokensetinput",
                                    "Token Set", 
                                    choices = c("")),
                       actionButton("preview",
                                    "Preview"),
                       actionButton("deletetokenset",
                                    "Delete",
                                    style="background-color: #b71616; color: white;")
                ),
                column(4,
                       awesomeRadio(inputId = "tokeninput",
                                    "Token", 
                                    choices = c(""))
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
                           checkboxGroupButtons("covtypes",
                                                "Select Covariate relationships",
                                                choices= c("None", "Linear","Power","Exponential","Proportional") ),
                           column(6,
                                  materialSwitch(inputId = "center", 
                                                 label = "Center Covariate (Median)", 
                                                 status = "success"),
                                  
                                  actionButton("addcovs",
                                               "Add Selected Token Sets")),
                           column(6,
                                  "Required {Token Group} in:", 
                                  HTML("<h4 id='covrequiredpk'>$PK</h4>
                                       <h4 id='covrequiredtheta'>$THETA</h4>"),
                                  class="required-holder")
                           ),
                       div(id = "ETA",
                           
                           class = "tokentype",
                           column(6,
                                  checkboxGroupButtons("etatypes",
                                                       "Select IIV relationships",
                                                       choices= c("None","Normal","Logarithmic","Normal (Proportional)") ),
                                  actionButton("addetas",
                                               "Add Selected Token Sets")),
                           column(6,
                                  "Required {Token Group} in:", 
                                  HTML("<h4 id='etarequiredpk'>$PK</h4>
                                       <h4 id='etarequiredomega'>$OMEGA</h4>"),
                                  class="required-holder")
                           
                           ),
                       div(id = "EPS",
                           class = "tokentype",
                           checkboxGroupButtons("epstypes",
                                                "Select Residual Error Models",
                                                choices= c("Additive","Proportional","Additive+Proportional") ),
                           actionButton("addeps",
                                        "Add Selected Token Sets")
                       ),
                       div(id = "Structure",
                           class = "tokentype",
                           radioGroupButtons("strtypes",
                                             "Select Structural Component",
                                             choices= c("Compartments","Tlag","Michaelis-Menten") ),
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
                                           radioGroupButtons("ncmt",
                                                             "Additional Compartments",
                                                             choices= c("1","2","3") ),
                                           actionButton("addstr",
                                                        "Add Selected Token Sets")))
                ),
                       div(id = "Custom",
                           class = "tokentype",
                           actionButton("addcustom",
                                        "Add Blank Token Set")
                           ))), class="token-div-holder"
         ),

  includeScript("www/jquery.highlighttextarea.js"),
  includeScript("www/drag.js"),
  includeCSS("www/style.css"),
  includeCSS("www/jquery.highlighttextarea.min.css"),
  
  
  # Add tooltips ------------------------------------------------------------
  
  bsTooltip(id = "dir", title = "Load Project"),
  bsTooltip(id = "saveproj", title = "Save Project"),
  bsTooltip(id = "proj", title = "Open in directory in File Explorer"),
  bsTooltip(id = "GAsettings", title = "Genetic Algorithm Settings"),
  bsTooltip(id = "SCMsettings", title = "SCM Settings"),
  bsTooltip(id = "deletetokengroup", title = "Delete selected token group"),
  bsTooltip(id = "preview", title = "Preview control stream with selected token set"),
  bsTooltip(id = "deletetokenset", title = "Delete selected token set")
  )

# Server ------------------------------------------------------------------

server <- function(input, output,session) {
  # Create a new project
  onclick("newproj",{
    toggle("new_project_div")
  })
  
  #these varaibles are used to trigger next generation when no more model tasks are found to be running
  invalidate <- reactiveValues(nextGA=1,
                               future=1)
  
  #copy model to new directory when dragged
  observe({
    if(!is.null(input$draggedfile[1])){
      copyTo <- paste0("models/", input$draggedfile[1])
      
      CopyModel(input$draggedfile[2], copyTo, control = input$ace, alltokens = alltokens, allmods = allmods)
    }
  }, suspended = F)
  
  
  # Define globals
  alltokens <- data.frame(tokengroup = character(0), 
                         tokenset = character(0), 
                         token = character(0), 
                         stringsAsFactors = FALSE)
  
  GAProgress <- data.frame (Generation = numeric(0), 
                            Fitness = numeric(0))
  
  f <- NULL
  onclick("dir", {
    directory <- choose.dir("M:\\Users\\mhismail-shared\\Rprogramming\\genetic algorithm\\")
    if (is.na(directory)){
      directory <- getwd()
    }
    
    setwd(directory)
    if(length(list.files(pattern = ".*.ctl"))>0){
      updateTextAreaInput(session,
                          inputId = "ace",
                          value = paste(readLines(list.files(pattern = ".*.ctl")[1]),
                                      collapse = "\n"))
      updateTextAreaInput(session,
                          inputId = "preview",
                          value = paste(readLines(list.files(pattern = ".*.ctl")[1]),
                                      collapse = "\n")) 
    }
    
    
    
    if (file.exists("alltokens.csv")){
      alltokens <- read.csv("alltokens.csv", as.is = T)
    }
    
    if (file.exists("Allmodsresults.csv.gz")){
      zz <- gzfile('Allmodsresults.csv.gz')   
      allmods <<- read.csv(zz, stringsAsFactors = F)[-1:-9] #model phenotypes
    }
    runjs({
      paste("$('#proj').text('Current Directory:", basename(directory), "')")
    })
  })
  
  onclick("proj",{
    shell(gsub("/", "\\\\", paste("explorer", getwd())))
  }
  )
  
  
  
  
  # Saving and loading projects ---------------------------------------------
  # Save project on click ---------------------------------------------------
  onclick("saveproj",{
    write.csv(alltokens, "alltokens.csv", row.names = F)
    writeLines(input$ace, list.files(pattern = ".*.ctl")[1])
  })
  
  # Hide compartment numbers if Compartments isnt selected ------------------
  observe({
    if(input$strtypes != "Compartments"){
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
    choices <- gsub("\\{|:(.*)|\\}", "", 
                   regmatches(x, gregexpr("\\{.*?\\}", x))[[1]])
    choices2 <- gsub("\\{|:(.*)|\\}", "", 
                     regmatches(toString(alltokens$token), gregexpr("\\{.*?\\}", toString(alltokens$token)))[[1]])
    choices <- unique(c(choices, choices2))
    if(length(choices) == 0){choices = c("")}
    updateAwesomeRadio(session,
                       inputId = "tokengroupinput",
                       choices = choices, selected = choices[1])
    
    # Parse control stream, search for "$INPUT", set covariate choices --------
    # Choices will be words following $INPUT, filtering out NONMEM reserved ---
    a <- strsplit(x, "\n")[[1]]
    if(length(strsplit(a[startsWith(a, "$INPUT")], "\\s+")) > 0){
      covlist <- strsplit(a[startsWith(a, "$INPUT")], "\\s+")[[1]][-1]
      covlist <- covlist[!(covlist %in% c("DATE=DROP", "DATE", "RATE", "CMT", "ID", "TIME", "AMT",
                                        "SS", "ADDL", "II", "DV", "MDV", "EVID", "DUR"))]
      updateSelectizeInput(session,
                           "cov",
                           choices = covlist)
    }
    
    # Parse control stream for compartments present
    b <- gsub(".*\\$MODEL(.*?)\\$.*", "\\1", x)
    c<- regmatches(b, gregexpr("COMP.*\\(.*?\\)", b))[[1]]
    d <- gsub(".*\\((.*)\\).*", "\\1", c)
    e<- gsub("\\,", "", d)
    updateSelectizeInput(session, "selectcmt", choices = e)
  }, priority = 100)
  # END ---------------------------------------------------------------------
  
  
  # Parse control stream to see if required tokens are present and --------
  #       change UI accordingly by adding/removie css class ---------------
  
  observe({
    x <- input$ace
    
    #check if required tokens are present
    selectedtokengroup <- paste0("\\{",input$tokengroupinput, "\\}")
    selectedtokengroupnumbered <-  paste0("\\{",input$tokengroupinput, ":[0-9]\\}")
    
    numberedtokengroups <- as.numeric(gsub(paste0("\\{", selectedtokengroup, ":|\\}"), 
                                           "",
                                           regmatches(x, gregexpr(paste0("\\{", selectedtokengroup, ":([0-9])\\}"), x))[[1]]))
    
    for (k in numberedtokengroups){
      x<- gsub(paste0("\\{", selectedtokengroup, ":", k, "\\}"), tokenlist[k], x)
    }
    
    PKblock <- gsub(".*\\$PK(.*?)[(\\$.*)$]", "\\1", x) # Text between $PK and next $
    THETAblock <- gsub(".*\\$THETA(.*?)[(\\$.*)$]", "\\1", x) # Text betwen $THETA and $
    OMEGAblock <- gsub(".*\\$OMEGA(.*?)[(\\$.*)$]", "\\1", x) # Text betwen $OMEGA and $
    
    
    if(grepl(selectedtokengroup, PKblock)|grepl(selectedtokengroupnumbered, PKblock)){
      removeClass("covrequiredpk", "red")
      addClass("covrequiredpk", "green")
    }else{
      removeClass("covrequiredpk", "green")
      addClass("covrequiredpk", "red")}
    
    if(grepl(selectedtokengroup, THETAblock) | grepl(selectedtokengroupnumbered, THETAblock)){
      removeClass("covrequiredtheta", "red")
      addClass("covrequiredtheta", "green")
    }else{
      removeClass("covrequiredtheta", "green")
      addClass("covrequiredtheta", "red")}
    
    if(grepl(selectedtokengroup, PKblock) | grepl(selectedtokengroupnumbered, PKblock)){
      removeClass("etarequiredpk", "red")
      addClass("etarequiredpk", "green")
    }else{
      removeClass("etarequiredpk", "green")
      addClass("etarequiredpk", "red")}
    
    if(grepl(selectedtokengroup, OMEGAblock) | grepl(selectedtokengroupnumbered, OMEGAblock)){
      removeClass("etarequiredomega", "red")
      addClass("etarequiredomega", "green")
    }else{
      removeClass("etarequiredomega", "green")
      addClass("etarequiredomega", "red")}
  }, priority = 100)
  
  # END ---------------------------------------------------------------------
  
  
  
  # filter tokenset list for selected token group ---------------------------
  observe({
    if(!is.null(input$tokengroupinput)){
      selectedTokenGroup <- input$tokengroupinput
      tokensetlist <- filter(alltokens, tokengroup == selectedTokenGroup)$tokenset %>% as.character()
      
      hide("tokeneditdiv", anim = T) #close editing tab on selection of new tokengroup
      updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
    }
  }, priority = -10)
  # END ---------------------------------------------------------------------
  
  
  # check if cov is in tokengroup and set to selected value if match is found ----------------------
  observe({
    selectedTokenGroup <- input$tokengroupinput
    x <- input$ace
    a <- strsplit(x, "\n")[[1]]
    if(length(strsplit(a[startsWith(a, "$INPUT")], "\\s+")) > 0){
      covlist <- strsplit(a[startsWith(a, "$INPUT")], "\\s+")[[1]][-1]
      covlist <- covlist[!(covlist %in% c("DATE=DROP", "DATE", "RATE", "CMT", "ID", "TIME", "AMT",
                                        "SS", "ADDL", "II", "DV", "MDV", "EVID", "DUR"))]
      bestguess <- covlist[str_detect(selectedTokenGroup, covlist)][which(nchar(covlist[str_detect(selectedTokenGroup, covlist)]) == max(nchar(covlist[str_detect(selectedTokenGroup, covlist)])))]
      updateSelectizeInput(session, "cov", selected = bestguess)
    }
  })
  # END ---------------------------------------------------------------------
  
  
  # filter token list for selected token set ---------------------------
  observe({
    if(!is.null(input$tokengroupinput) & !is.null(input$tokensetinput)){
      selectedTokenGroup <- input$tokengroupinput
      selectedTokenSet <- input$tokensetinput
      tokenlist <- filter(alltokens, tokengroup == selectedTokenGroup, tokenset == selectedTokenSet)$token %>% as.character()
      updateAwesomeRadio(session, inputId = "tokeninput", choices = tokenlist, selected = unique(tokenlist)[1])
    }
  }, priority = -11)
  # END ---------------------------------------------------------------------
  
  # highlight selected tokengroup -------------------------------------------
  # highlight package used: https://github.com/garysieling/jquery-highlighttextarea
  onclick("tokengroupinput",
          {
            updateTabsetPanel(session,"tabs", selected = "Control Stream")
            
            runjs("$('#ace').highlightTextarea('destroy')")
            fun<-paste0("$('#ace').highlightTextarea({
                        words: [{
                        color: '#FFFF00',
                        words: ['{", input$tokengroupinput, "}',
                        '{", input$tokengroupinput, ":[0-9]}']
                        }
                        , {
                        color: '#ff1b0f',
                        words: []
                        }]
          });")
          
            
            
            delay(1, runjs(fun)) #needs delay for tokengroup change to register first
            })
  
  # END ---------------------------------------------------------------------
  
  
  # show how token set appears in control stream ----------------------------
  onclick("preview",
          {
            x <- input$ace
            selectedTokenGroup <- input$tokengroupinput
            selectedTokenSet <- input$tokensetinput
            tokengrouptext <- paste0("\\{", input$tokengroupinput, "\\}")
            tokenlist <- filter(alltokens, tokengroup == selectedTokenGroup, tokenset == selectedTokenSet)$token %>% as.character()
            tokenlist[tokenlist == "N/A"] <- ""
            for (i in 1:length(tokenlist)){
              x <- sub(tokengrouptext, tokenlist[i], x)
            }
            x <- gsub("\\{[^//}]*\\}", "", x)
            updateTextAreaInput(session, "previewtext", value = x)
            words <- paste0("'", tokenlist, "',", collapse = "")
            words <- substr(words, 1, nchar(words) - 1)
            
            updateTabsetPanel(session, "tabs", selected = "Preview")
          })
  # END ---------------------------------------------------------------------
  
  
  
  # display selected token type div (ie covariate, ETA,EPS,etc.) ------------
  onclick("tokentype",
          {
            hide(selector = ".tokentype")
            showElement(id = input$tokentype)
          })
  # END ---------------------------------------------------------------------
  
  
  
  # Delete functions -------------------------------------------------------
  onclick("deletetokengroup",
          {
            x <- input$ace
            x <- gsub(paste0("\\{", input$tokengroupinput, "\\}"), "", x)
            alltokens <- filter(alltokens, tokengroup != input$tokengroupinput)
            updateTextAreaInput(session, "ace", value = x)
          })
  
  onclick("deletetokenset",
          {
            alltokens <- filter(alltokens, !(tokengroup == input$tokengroupinput & tokenset == input$tokensetinput))
            tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
            updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
          })
  # END ---------------------------------------------------------------------
  
  
  
  # Add to existing tokengroup ----------------------------------------------
  
  onclick("addtokensetedit", {
    n_tokensets <- length(unique(filter(alltokens, 
                              tokengroup == input$tokengroupinput)$tokenset))
    

    n_tokens <- length(filter(alltokens, 
                              tokengroup == input$tokengroupinput, 
                              tokenset == input$tokensetinput)$token)
    
    
    tokenToAdd <- data.frame(tokengroup = input$tokengroupinput, tokenset = paste("Tokenset", n_tokensets + 1), token = rep("N/A", n_tokens))
    
    # If tokenset exists already, replace it
    alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
    alltokens <- rbind(alltokens, tokenToAdd)
    
    tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
    updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
    
    
    tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
    updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
    
    selectedtokengroup <- input$tokengroupinput
    removeUI(selector = ".tokeneditdiv", immediate = T) #.tokeneditdiv class is tokenset panel div
    
    #only display panel if tokensets exist for selected tokengroup
    if (length(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput]) > 0){
      insertUI(selector = "#tokenedit", where = "beforeEnd",
               ui = {
                 column(12, align = "center",
                        textInput("edittokengroupname", "Token Group Name", value = selectedtokengroup,width = "33%"),
                        lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])), function(i) {
                          column(12,
                                 textInput(paste0("tokenset", i),"Token Set Name", value = unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i], width = "33%"),
                                 lapply(seq_along(alltokens$token[alltokens$tokengroup == input$tokengroupinput & alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i]]),function(j){
                                   column(6,
                                          textAreaInput(paste0("edittoken", i, j), paste("Token", j), value = alltokens$token[alltokens$tokengroup == input$tokengroupinput & 
                                                                                                                           alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i]][j],
                                                        resize = "vertical"))}),
                                 class = "tokeneditdiv")}),
                        column(12,
                               actionButton("addtokensetedit", "Add Token Set"),
                               actionButton("addtokenedit", "Add Token"),
                               actionButton("savetokenedit", label = "Save")), class = "tokeneditdiv", id = "tokeneditdiv")
               }, session = session, immediate = T)
    }
  })
  
  # END---------------------------------------------------------------------
  
  # Add token to existing tokenset ----------------------------------------------
  
  onclick("addtokenedit",{
    
    tokensets <- unique(filter(alltokens,tokengroup == input$tokengroupinput)$tokenset)
    tokenToAdd <-data.frame(tokengroup = input$tokengroupinput, tokenset = tokensets, token = "N/A")
    
    alltokens<-rbind(alltokens, tokenToAdd)
    
    tokensetlist<- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
    updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
    
    
    tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
    updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
    
    selectedtokengroup <- input$tokengroupinput
    removeUI(selector = ".tokeneditdiv", immediate = T) #.tokeneditdiv class is tokenset panel div
    
    #only display panel if tokensets exist for selected tokengroup
    if (length(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput]) > 0){
      insertUI(selector = "#tokenedit", where = "beforeEnd",
               ui = {
                 column(12, align = "center",
                        textInput("edittokengroupname", "Token Group Name", value = selectedtokengroup, width = "33%"),
                        lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])), function(i) {
                          column(12,
                                 textInput(paste0("tokenset",i),"Token Set Name",value = unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i] ,width = "33%"),
                                 lapply(1:length(alltokens$token[alltokens$tokengroup == input$tokengroupinput & alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i]]),function(j){
                                   column(6,
                                          textAreaInput(paste0("edittoken", i, j), paste("Token", j), value = alltokens$token[alltokens$tokengroup == input$tokengroupinput & 
                                                                                                                           alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i]][j],
                                                        resize = "vertical"))}),
                                 class = "tokeneditdiv")}),
                        column(12,
                               actionButton("addtokensetedit", "Add Token Set"),
                               actionButton("addtokenedit", "Add Token"),
                               actionButton("savetokenedit", label = "Save")), class = "tokeneditdiv", id = "tokeneditdiv")
               }, session = session, immediate = T)
    }
  })
  
  # END---------------------------------------------------------------------
  
  
  # Create and add covariate tokens  ----------------------------------------
  
  onclick("addcovs", {
          tokenTypes <- data.frame(tokengroup = input$tokengroupinput, tokenset = input$covtypes, stringsAsFactors = FALSE)
          x <- input$ace
          a <- strsplit(x,"\n")[[1]]
          covlist <- strsplit(a[startsWith(a, "$INPUT")], "\\s+")[[1]][-1]
          
          if(any(covlist == input$cov)){
            values <- as.numeric(datafile()[which(covlist == input$cov)][[1]])
            median <- median(values,na.rm = T) %>% round(2)
            
            if (input$center == T){
              
              lookup <- data.frame(tokenset = rep(c("None", "Linear", "Power", "Exponential", "Proportional"), each = 2),
                                   token = c("N/A",
                                            "N/A",
                                            paste0("+(", input$cov, "-", median, ")*THETA([1])"),
                                            paste0("(-10,.001,10);", input$tokengroupinput),
                                            paste0("*(",input$cov, "/", median, ")**THETA([1])"),
                                            paste0("(-5,.001,10);", input$tokengroupinput),
                                            paste0("*exp((", input$cov, "-", median, ")*THETA([1]))"),
                                            paste0("(-50,.001,50);", input$tokengroupinput),
                                            paste0("*(1+(", input$cov, "-", median, ")*THETA([1]))"),
                                            paste0("(-100,.001,100);", input$tokengroupinput)), stringsAsFactors=FALSE)
              
            }
          }
          
          if (input$center == F){
            
            lookup <- data.frame(tokenset = rep(c("None", "Linear", "Power", "Exponential", "Proportional"), each = 2),
                                 token= c("N/A",
                                          "N/A",
                                          paste0("+(", input$cov, ")*THETA([1])"),
                                          paste0("(-100,.001,100);", input$tokengroupinput),
                                          paste0("*", input$cov, "**THETA([1])"),
                                          paste0("(-100,.001,100);", input$tokengroupinput),
                                          paste0("*exp(", input$cov, ")*THETA([1])"),
                                          paste0("(-100,.001,100);", input$tokengroupinput),
                                          paste0("*(1+", input$cov, "*THETA([1]))"),
                                          paste0("(-100,.001,100);", input$tokengroupinput)), stringsAsFactors = FALSE)
            
          }
          
          
          
          
          
          tokenToAdd <- merge(tokenTypes, lookup, by = "tokenset") %>% select(tokengroup, tokenset, token)
          
          # If tokenset exists already, replace it
          alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
          alltokens <-rbind(alltokens, tokenToAdd)
          
          tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
          updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
          
          tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
          updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
          hide("Covariate", anim = T)
          }
  )
  
  # END ---------------------------------------------------------------------
  
  
  # Create and add eta tokens  ----------------------------------------
  
  onclick("addetas",
          {
            tokenTypes <- data.frame(tokengroup = input$tokengroupinput, tokenset = input$etatypes)
            lookup <- data.frame(tokenset = rep(c("None", "Normal", "Logarithmic", "Normal (Proportional)"), each=2),
                                 token= c("N/A",
                                          "N/A",
                                          paste0("+ETA([1])"),
                                          paste0("(0.01);", input$tokengroupinput),
                                          paste0("*exp(ETA([1]))"),
                                          paste0("(0.01);", input$tokengroupinput),
                                          paste0("*(1+ETA([1]))"),
                                          paste0("(0.01);", input$tokengroupinput)), stringsAsFactors = FALSE)
            tokenToAdd <- merge(tokenTypes, lookup, by = "tokenset") %>% select(tokengroup, tokenset, token)
            alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
            
            alltokens <- rbind(alltokens, tokenToAdd)
            
            tokensetlist<- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
            updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
            
            
            tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
            updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
            hide("ETA", anim = T)
          }
  )
  # END ---------------------------------------------------------------------
  
  
  # Create and add eps tokens  ----------------------------------------
  
  onclick("addeps",
          {
            tokenTypes <- data.frame(tokengroup = input$tokengroupinput, tokenset = input$epstypes)
            lookup <- data.frame(tokenset = rep(c("Additive", "Proportional", "Additive+Proportional"), each=2),
                                 token = c(paste0("Y=F+EPS([1])"),
                                          "(0.01); Additive",
                                          paste0("Y=F*(1+EPS([1]))"),
                                          "(0.01); Proportional",
                                          paste0("Y=F*(1+EPS([1]))+EPS([2])"),
                                          "(0.01); Additive \n (0.01);Proportional"), stringsAsFactors = FALSE)
            tokenToAdd <-merge(tokenTypes, lookup, by = "tokenset") %>% select(tokengroup, tokenset, token)
            alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
            
            alltokens <- rbind(alltokens, tokenToAdd)
            
            tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
            updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
            
            
            tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
            updateAwesomeRadio(session, inputId = "tokeninput",choices = as.character(tokenlist), selected = NULL)
            
            hide("EPS", anim = T)
          }
  )
  # END ---------------------------------------------------------------------
  
  
  # Create and add structure tokens  ----------------------------------------
  onclick("addstr",
          {
            
            if(input$strtypes == "Compartments"){
              
              b <- gsub(".*\\$MODEL(.*?)\\$.*", "\\1", input$ace)
              c <- regmatches(b, gregexpr("COMP.*\\(.*?\\)", b))[[1]]
              d <- gsub(".*\\((.*)\\).*", "\\1", c)
              e <- gsub("\\,", "", d)
              
              cmtn <- which(e == input$selectcmt)
              ncmts <- length(e)
              
              tokenTypes <- data.frame(tokengroup = input$tokengroupinput, tokenset = c("None", paste0("+", 1:input$ncmt, "cmt")))
              lookup <- data.frame(tnumber = rep(c(1, 2, 3), 4), tokenset = rep(c("None", "+1cmt", "+2cmt", "+3cmt"), each = 3),
                                   token = c("N/A", "N/A", "N/A",
                                            as.character(p1cmttext(input$selectcmt, cmtn, ncmts)[1]),
                                            as.character(p1cmttext(input$selectcmt, cmtn, ncmts)[2]),
                                            "(0 10)\n(0 10)",
                                            as.character(p2cmttext(input$selectcmt, cmtn, ncmts)[1]),
                                            as.character(p2cmttext(input$selectcmt, cmtn, ncmts)[2]),
                                            "(0 10)\n(0 10)\n(0 10)\n(0 10)",
                                            as.character(p3cmttext(input$selectcmt, cmtn, ncmts)[1]),
                                            as.character(p3cmttext(input$selectcmt, cmtn, ncmts)[2]),
                                            "(0 10)\n(0 10)\n(0 10)\n(0 10)\n(0 10)\n(0 10)"), stringsAsFactors = FALSE)
              tokenToAdd <- merge(lookup, tokenTypes, by = "tokenset") %>% arrange(tokenset, tnumber) %>% select(tokengroup, tokenset, token)
              alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
              alltokens <- rbind(alltokens, tokenToAdd)
              
              tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
              updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
              
              
              tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
              updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
              
              hide("Structure", anim = T)
            }
          }
  )
  
  onclick("addcustom",{
    tokenToAdd <- data.frame(tokengroup = input$tokengroupinput, tokenset = "New Token", token = "N/A")
    alltokens <- filter(alltokens, !(tokengroup %in% tokenToAdd$tokengroup & tokenset %in% tokenToAdd$tokenset))
    alltokens <-rbind(alltokens, tokenToAdd)
    
    tokensetlist <- filter(alltokens, tokengroup == input$tokengroupinput)$tokenset %>% as.character()
    updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
    
    
    tokenlist <- filter(alltokens, tokenset == unique(tokensetlist)[1])$token %>% as.character()
    updateAwesomeRadio(session, inputId = "tokeninput", choices = as.character(tokenlist), selected = NULL)
    
    hide("Custom", anim = T)
  }
  )
          
  # END ---------------------------------------------------------------------
  
  
  # Show panel of tokensets when double click token group  ------------------
  
  onevent("dblclick", "tokengroupinput", {
    
    selectedtokengroup <- input$tokengroupinput
    removeUI(selector= ".tokeneditdiv", immediate=T) #.tokeneditdiv class is tokenset panel div
    
    #only display panel if tokensets exist for selected tokengroup
    if (length(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput]) > 0){
      insertUI(selector="#tokenedit", where = "beforeEnd",
               ui = {
                 column(12, align = "center",
                        textInput("edittokengroupname", "Token Group Name", value = selectedtokengroup, width = "33%"),
                        lapply(1:length(unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])), function(i) {
                          column(12,
                                 textInput(paste0("tokenset", i), "Token Set Name", value = unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i], width = "33%"),
                                 lapply(1:length(alltokens$token[alltokens$tokengroup == input$tokengroupinput & alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$tokengroupinput])[i]]), function(j){
                                   column(6,
                                          textAreaInput(paste0("edittoken", i, j), paste("Token", j), value = alltokens$token[alltokens$tokengroup == input$tokengroupinput & 
                                                                                                                           alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup==input$tokengroupinput])[i]][j],
                                                        resize = "vertical"))}),
                                 class = "tokeneditdiv")}),
                        column(12,
                               actionButton("addtokensetedit", "Add Token Set"),
                               actionButton("addtokenedit", "Add Token"),
                               actionButton("savetokenedit", label = "Save")), class = "tokeneditdiv", id = "tokeneditdiv")
               }, session = session, immediate = T)
    }
  })
  
  # END ---------------------------------------------------------------------
  
  
  # Save tokens after editing  ----------------------------------------------
  
  onclick("savetokenedit",{
    alltokens$tokengroup[alltokens$tokengroup == input$tokengroupinput] <- input$edittokengroupname
    
    for (i in 1:length(unique(alltokens$tokenset[alltokens$tokengroup == input$edittokengroupname]))){
      alltokens$tokenset[alltokens$tokengroup == input$edittokengroupname &
                           alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$edittokengroupname])[i]] <- input[[paste0("tokenset", i)]]
      
      
      for (j in 1:length(alltokens$token[alltokens$tokengroup == input$edittokengroupname & alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$edittokengroupname])[i]])){
        alltokens$token[alltokens$tokengroup == input$edittokengroupname & 
                          alltokens$tokenset == unique(alltokens$tokenset[alltokens$tokengroup == input$edittokengroupname])[i]][j] <- input[[paste0("edittoken", i, j
                          )]]
      }
      
    }
    
    x <- input$ace
    tokengrouptext <- paste0("\\{", input$tokengroupinput, "\\}")
    
    x <- gsub(tokengrouptext, paste0("\\{", input$edittokengroupname, "\\}"), x)
    
    
    updateTextAreaInput(session, "ace", value = x)
    
    tokensetlist <- filter(alltokens, tokengroup == input$edittokengroupname)$tokenset %>% as.character()
    updateAwesomeRadio(session, inputId = "tokensetinput", choices = unique(tokensetlist), selected = unique(tokensetlist)[1])
    
    tokenlist <- filter(alltokens, tokengroup == input$edittokengroupname, tokenset == unique(tokensetlist)[1])$token %>% as.character()
    updateAwesomeRadio(session, inputId = "tokeninput", choices = unique(tokenlist), selected = unique(tokenlist)[1])
    
    
    hide("tokeneditdiv", anim = T)
    
  }
  )
  
  # END ---------------------------------------------------------------------
  
  
  
  onclick("refreshmods", {
    
    #delete models folder if it exists
    unlink("models/All", force = TRUE, recursive = T)
    
    #create new models folder
    if (!dir.exists("models")){
      dir.create("models")
    }
    dir.create("models/All")
    
    allmods <- AllCombos(alltokens)
    
    allmods[] <- lapply(allmods, as.character) 
    allmods <<- allmods # define allmods globally
    
    ctlstream <- input$ace
    #copy data to subfolder
    a <- strsplit(ctlstream, "\n")[[1]]
    datapath <- strsplit(a[startsWith(a, "$DATA")], "\\s+")[[1]][2]
    file.copy(datapath, paste0("models"))
    
    
    alltokens$token <- as.character(alltokens$token) 
    alltokens$tokenset <- as.character(alltokens$tokenset) 
    
    ctlstream <- gsub("\\$DATA\\s*", "\\$DATA ../../", ctlstream) 
    CreateModel <- CreateModel 
    z <- Sys.time()
    nmods <- dim(allmods)[1]
    a <- cbind(data.frame(Number = 1:nmods, Path = paste0("models/All/mod", 1:nmods, "/mod.ctl"), OFV="", Fitness="", S="", C="", NTHETA="", NETA="", NEPS="", allmods))
    
    
    gz1 <- gzfile("Allmodsresults.csv.gz", "w")
    write.csv(a, gz1, row.names=F)
    close(gz1)
    # write.csv(a,"Allmodsresults.csv",row.names = F)
  })
  
  
  
  # Create table to be displayed in UI
  allmods2<- reactive({
    z <- Sys.time()
    input$refreshmods2     
    input$refreshmods
    
    selecteddir <- input$selecteddir
    
    if (selecteddir == "All"){
      zz = gzfile('Allmodsresults.csv.gz')   
      allmodsresults <- read.csv(zz, stringsAsFactors = F)
      # allmodsresults <- read.csv("Allmodsresults.csv")
      
      #if (models have been run since laste refresh)
      if (file.exists("modelruntemp.csv") & file.info("modelruntemp.csv")$size > 0){ 
        modelsran <- read.csv("modelruntemp.csv", header = F) %>% filter(!duplicated(V1))
        checkpresent <- file.exists(paste0(modelsran[,1], "/results/PsN_execute_plots.R"))
        resultspresent <- modelsran[checkpresent,]
        if (dim(resultspresent)[1] > 0){
          results <- lapply(resultspresent[,1], RetrieveResults)
          removerefreshed <- modelsran[!checkpresent,]
          write.table(removerefreshed,
                      "modelruntemp.csv", sep = ",", row.names = F, col.names = F)
          resultsdf <- data.frame()
          resultsdf <- rbind(resultsdf, do.call(rbind, results))
          resultsdf[,1] <- as.numeric(as.character(resultsdf[,1]))
          mods <- resultsdf[,1]
          match(allmodsresults$Number, mods)
          
          #logic for below - match function finds the indices of allmodsresults in which the number equals the 
          # value in mods. na.omit removes na values. 1,3,5,6 are columns to be updated
          # right side !is.na is safegaurd for models that were run but have since been deleted before results were fetched
          allmodsresults[na.omit(match(mods, allmodsresults$Number)), c(1, 3, 5, 6)] <- resultsdf[!is.na(match(mods, allmodsresults$Number)),]
          
          gz1 <- gzfile("Allmodsresults.csv.gz", "w")
          write.csv(allmodsresults, gz1, row.names = F)
          close(gz1)
          
          # write.csv(allmodsresults,"Allmodsresults.csv",row.names = F)
        }
        
      }
    }else{
      allmods1 <- data.frame ()
      allmodslist <-list()
      selecteddir <- input$selecteddir
      mod.directories <- paste0("./models/", selecteddir, "/", list.files(paste0("./models/", selecteddir)))
      
      nmods <- length(mod.directories)
      allmodslist <- lapply(mod.directories, RetrieveResultsEach)
      
      allmods1 <- rbind(allmods1, do.call(rbind, allmodslist))
      allmodsresults <- mutate(allmods1, Number = as.numeric(Number)) %>% arrange(Number)
    }
    
    
    updateRadioGroupButtons(session, "selecteddir", selected = selecteddir, choices = list.dirs("./models", full.names = F, recursive = F))
    delay(100, runjs('$(".radiobtn").on("dragenter",function () {
                    y=$(this).find("input").val()
                    console.log("hi")
  })'))
    allmodsresults
    })
  
  
  
  # call PsN execute to run model
  onclick("runmod", {
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    for (i in 1:length(input$modeltable_selected)){
      mod <- input$modeltable_selected[i]
      all <- allmods2()
      path <- as.character(all[all$Number == mod, 2]) # Column 2 is model path
      relpath <- sub("/mod.ctl", "", path)
      
      CheckThenCreate(relpath, input$ace, alltokens, allmods)
      unlink(paste0(homepath,'/', relpath, "/results"), recursive = TRUE) 
      RunModel(paste0(homepath, '/', relpath), basename(relpath))
    }
    
  })
  
  # Delete selected model/s
  onclick("deletemods",{
    
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    for (i in 1:length(input$modeltable_selected)){
      all <- allmods2()
      
      path<-as.character(all[all$Number == input$modeltable_selected[i],2])
      if(input$selecteddir == "All"){
        alert("Cannot delete from 'All' directory")
      }
      relpath <- sub("/mod.ctl", "", path)
      unlink(paste0(homepath, '/', relpath), recursive = TRUE)
    }
    
    runjs('$("#refreshmods2").click()')
  })
  
  # open file explorer
  onclick("openlocation", {
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    CheckThenCreate(relpath, input$ace, alltokens, allmods)
    shell(gsub("/", "\\\\", paste("explorer", relpath)))
    
  })
  
  onclick("addsimilar", {
    all <- allmods2() 
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    
    phen <- read.csv(paste0(relpath, "/mod.csv"), as.is = T)[-(1:9)]
    zz <- gzfile('Allmodsresults.csv.gz')   
    allmods <- read.csv(zz, stringsAsFactors = F)[-1:-9]
    # allmods <- read.csv("Allmodsresults.csv",as.is=T)[-1:-9]
    
    similarmods <- which(apply(allmods, 1, function(x) sum( x == phen) >= (dim(allmods)[2] - 1)))
    
    copyTo <- paste0("models/",input$selecteddir)
    copy <- paste0("models/All/mod",similarmods)
    for (i in copy){
      CopyModel(i,copyTo,control=input$ace,alltokens=alltokens,allmods=allmods)
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
  filter = list(position = "top"),
  options = list(
    select = list(style = 'os', # set 'os' select style so that ctrl/shift + click in enabled
                  items = 'row'), # items can be cell, row or column
    pageLength = 100,
    dom = 'Bfrtip', 
    buttons = I('colvis'),
    columnDefs = list(list(className = 'dt-center', targets = "_all"))),
  rownames = FALSE, 
  extensions = c('Select', "Buttons"), 
  selection = "none", 
  callback = JS(
    "table.on( 'click.dt', 'tbody td', function (e) {", # react on click event
    "var type = table.select.items();", # get the items setting of Select extension
    "var idx = table[type + 's']({selected: true}).indexes().toArray();
    var values = [];
    for (var i = 0 ; i<idx.length; i++){
    values.push(table.table().column(0).data()[idx[i]])
    }", # get the index of selected items
    "var DT_id = table.table().container().parentNode.id;", # get the output id of DT
    "Shiny.onInputChange('modeltable' + '_selected', values);", # send the index to input$outputid_selected
    "})"
  ), server = T)
  
  table <- modalDialog(fluidPage(div(column(9,
                                             actionButton("runmod", "Run Selected Model"),
                                             actionButton("refreshmods", "Recreate Models"),
                                             actionButton("refreshmods2", "Refresh"),
                                             actionButton("deletemods", "Delete Selected Models"),
                                             actionButton("editmod", "Edit Selected Model"),
                                             actionButton("openlocation", "Open Model Location"),
                                             actionButton("addsimilar", "Add Similar Models"),
                                             DT::dataTableOutput('modeltable'),
                                             
                                             radioGroupButtons("selecteddir", 
                                                               NULL, 
                                                               choices = c("All"), 
                                                               selected = "All"),
                                             actionButton("showadddir", icon("plus", "fa-1x")),
                                             span(id = "add-directory", 
                                                  style = "display: none;", 
                                                  div(style="display: inline-block;vertical-align:top; width: 150px;", 
                                                      textInput("boog", NULL, "New", width = "auto")), 
                                                  actionButton("adddir", "Add")),
                                            
                                             ondrop="drop(event)", ondragover="allowDrop(event)")),
                                  div(column(3, 
                                             tabsetPanel(id = "modelinfo",
                                                            tabPanel("Plots",
                                                                     bsCollapse(bsCollapsePanel("Covariate Model",
                                                                                                actionButton("ranparvscov","ETAs vs Covariates"),
                                                                                                actionButton("parvscov","Parameters vs Covariates"),
                                                                                                actionButton("parvspar","Parameters vs Parameters"),
                                                                                            
                                                                                                actionButton("covscatter","Covariate Scatter Plots")),
                                                                                bsCollapsePanel("Goodness of Fit",
                                                                                                actionButton("basicgof","Basic GOF"),
                                                                                                actionButton("indplots","Individual Plots"),
                                                                                                actionButton("dvpredipred","DV vs Pred, IPRED")),
                                                                                multiple = T, open = c("Covariate Model"))
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
                                                                     open = c("Structural (THETA)", "IIV (OMEGA)", "RUV (SIGMA)"))),
                                                            tabPanel("Covariate Model",
                                                                     tableOutput("regress.parms")),
                                                            tabPanel("PsN",
                                                                     actionButton("runvpc", "VPC")))))), easyClose = T)
  
  onclick("showadddir", {
    toggle("add-directory", anim = T, animType = "fade")
  })
  
  onclick("adddir",{
    dir.create(paste0(getwd(), "/models/", input$boog))
    hide("add-directory", anim = T, animType = "fade")
    
    updateRadioGroupButtons(session, 
                            "selecteddir", 
                            selected = input$selecteddir,
                            choices = list.dirs("./models", full.names = F, recursive = F))
    delay(100, runjs('$(".radiobtn").on("dragenter",function () {
                    y=$(this).find("input").val()})'))
    })
  

  observeEvent(input$startGA, {invalidate$future = isolate(invalidate$future) + 1})
  observeEvent(input$nextGA, {invalidate$future = isolate(invalidate$future) + 1})
  
  
  observeEvent({
    invalidate$future
  }, {
    print("started")
    delay(10, {
      f <<- future({
        b <- 1
        while(length(b) > 0){
          a <- shell("tasklist /v",intern=T)
          b <- as.vector(na.omit(str_extract(a, "mod[0-9]* - execute")))
        }}) %plan% multiprocess
      invalidate$nextGA <- isolate(invalidate$nextGA) + 1})
  }, ignoreInit = T)
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
        runjs('$("#refreshGAmods").click()')
        runjs('$("#nextGA").click()')
      }
    }
    
  })
  
  
  
  
  # Genetic Algorithm Modal -------------------------------------------------
  
  GAmods<- reactive({
    input$refreshGAmods
    
    
    allmods1 <- data.frame ()
    allmodslist <- list()
    mod.directories <- read.csv("GAgens.csv", as.is = T)%>%
      filter(Generation == as.numeric(input$selectedgen)) %>% select(Dir)
    
    
    allmodslist <- lapply(paste0("./", mod.directories$Dir), RetrieveResultsEach)
    
    
    allmods1 <- rbind(allmods1, do.call(rbind, allmodslist))
    allmods1 <- mutate(allmods1, Number = as.numeric(Number)) %>% arrange(Number)
    
    
    allmods1
    
    
  })
  
  GAmodsAllGens<- reactive({
    input$refreshGAmods
    
    
    allmods1 <- data.frame ()
    allmodslist <- list()
    mod.directories <- read.csv("GAgens.csv", as.is = T)
    mod.directories.max.gen <- filter(mod.directories, Generation == max(Generation))
    
    
    allmodslist <- lapply(paste0("./", mod.directories.max.gen$Dir), RetrieveResultsEach)
    
    
    allmods1 <- rbind(allmods1, do.call(rbind, allmodslist))
    allmods1 <- mutate(allmods1, Number = as.numeric(Number), Generation = mod.directories.max.gen$Generation) %>%
      arrange(Number) %>%
      select(Generation, Fitness)
    GAProgress <- rbind(GAProgress, allmods1)
    
    
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
  rownames = FALSE, 
  extensions = 'Select', 
  selection = "none", 
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
    summary <- summarize(GAmodsAllGensGrouped, meanfit = mean(as.numeric(Fitness),na.rm=T), 
                         minfit = min(Fitness, na.rm=T))
    
    ggplot(summary, aes(x=Generation,y=minfit), color="blue") + 
      ylim(min(summary$minfit), mean(summary$minfit)*1.1) + 
      geom_line() + 
      geom_line(aes(y=meanfit), color="red") 
  })
  
  observeEvent(input$openGA, {
    showModal(
      GAmodal
    )
  })
  
  onclick("startGA", {
    GAProgress <- data.frame (Generation = numeric(0), 
                              Fitness = numeric(0))
    
    nmods <- dim(allmods)[1]
    nindiv <- 20
    if (nindiv > nmods) nindiv <- nmods
    InitiateGA(nmods, nindiv, control = input$ace, alltokens = alltokens, allmods = allmods, seed = 10)
  })
  
  onclick("nextGA",{
    
    maxgen <- NextGA(alltokens, control = input$ace, allmods = allmods)
    updateRadioGroupButtons(session, "selectedgen", choices = 1:maxgen)
  })
  
  GAmodal <-  modalDialog(fluidPage(div(column(9,
                                               actionButton("startGA", "Start"),
                                               actionButton("nextGA", "Next"),
                                               actionButton("refreshGAmods", "Refresh"),
                                               actionButton("openlocation2", "Open Model Location"),
                                               DT::dataTableOutput('GAtable'))),
                                    div(column(3, tabsetPanel(id = "modelinfo2",
                                                              tabPanel("Plots",
                                                                       bsCollapse(bsCollapsePanel("Covariate Model",
                                                                                                  actionButton("ranparvscov2", "ETAs vs Covariates"),
                                                                                                  actionButton("covscatter2", "Covariate Scatter Plots")),
                                                                                  bsCollapsePanel("Goodness of Fit",
                                                                                                  actionButton("basicgof2", "Basic GOF"),
                                                                                                  actionButton("indplots2", "Individual Plots"),
                                                                                                  actionButton("dvpredipred2", "DV vs Pred, IPRED")),
                                                                                  multiple = T, open = c("Covariate Model")),
                                                                       plotOutput("GAPlot", height = "200px")
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
                                                                       open = c("Structural (THETA)", "IIV (OMEGA)", "RUV (SIGMA)"))),
                                                              tabPanel("PsN",
                                                                       actionButton("runvpc2", "VPC"))))),
                                    column(12, radioGroupButtons("selectedgen", NULL, choices = c(1), selected = NULL))), easyClose = T)
  
  # GA settings modal -------------------------------------------------------
  GAsettingsmodal <- modalDialog(fluidPage(
    radioGroupButtons("settingsselect",
                      choices = c("Fitness", "Selection", "Crossover", "Mutation")),
    div(id = "Fitness",
        class = "GA-settings-div",
        numericInput("THETA","Theta Penalty", value = 5, width = "60%"),
        numericInput("ETA", "Eta Penalty", value = 5, width = "60%"),
        numericInput("EPS", "Epsilon Penalty", value = 5, width = "60%"),
        numericInput("covsuccess", "Successful Covariance", value = -10, width = "60%"),
        numericInput("minsuccess", "Successful Minimization", value = -10, width = "60%")
        
        #         numericInput("covwarn","Covariance Warning",value = 0),
        #         numericInput("boundary","Parameter/s Near Boundary",value = 0),
        #         numericInput("rounding","Rounding Errors",value = 0),
        #         numericInput("zerograd","Zero Gradient",value = 0),
        #         numericInput("finalzerograd","Final Zero Gradient",value = 0),
        #         numericInput("hessianreset","Hessian Reset",value = 0),
        #         numericInput("ssingular","S-Matrix Singular",value = 0)
        
    ),
    div(id = "Selection",
        class = "GA-settings-div",
        radioButtons("selectionmethod", "Selection Algorithm", choices = c("2-Way Tournament", "Scaled Roullette", "Rank Roullette"))
    ),
    div(id = "Crossover",
        class = "GA-settings-div",
        numericInput("crossoverrate", "Crossover Rate", value = 0.7, min = 0, max = 1, width = "60%")
    ),
    div(id = "Mutation",
        class = "GA-settings-div",
        numericInput("mutationrate", "Mutation Rate", value = 0.05, min = 0, max = 1, width = "60%")
    )
  ),
  size = "s",
  class = "GAsettingsmodal")
  
  observeEvent(input$GAsettings, {
    showModal(
      GAsettingsmodal
    )
  })
  
  onclick("settingsselect",
          {
            hide(selector = ".GA-settings-div")
            showElement(id = input$settingsselect)
          })
  
  
  # Edit control stream modal -----------------------------------------------
  
  editcontrolstream <- modalDialog(fluidPage(actionButton("savecontrol", "Save"), textAreaInput("editcontrol",
                                                                                              label = NULL,
                                                                                              width = "100%",
                                                                                              height = "200px",
                                                                                              cols = NULL,
                                                                                              value = "")))
  
  
  onclick("savecontrol", {
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    x <- input$editcontrol
    writeLines(x, paste0(path))
  })
  
  
  
  observeEvent(input$editmod, {
    all <- allmods2()
    
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA/"
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    file <- paste0(homepath, "/", path)
    relpath <- sub("/mod.ctl", "", path)
    
    CheckThenCreate(relpath, input$ace, alltokens, allmods)
    if (file.exists(file)){
      updateTextAreaInput(session, "editcontrol", value = paste(readLines(file), collapse = "\n"))
    }
    showModal(
      editcontrolstream
    )})
  
  
  #  ------------------------------------------------------------------------
  
  # SCM Modal --------------------------------------------------------
  observeEvent(input$openSCM, {
    showModal(
      scm_table
    )
  })
  
  # Render table
  output$modeltable3 <- DT::renderDataTable({
    allmods2()
  },
  filter = list(position = "top"),
  options = list(
    select = list(style = 'os', # set 'os' select style so that ctrl/shift + click in enabled
                  items = 'row'), # items can be cell, row or column
    pageLength = 100,
    dom = 'tpf', 
    columnDefs = list(list(className = 'dt-center', targets = "_all"))),
  rownames = FALSE, 
  extensions = 'Select', 
  selection = "none", 
  callback = JS(
    "table.on( 'click.dt', 'tbody td', function (e) {", # react on click event
    "var type = table.select.items();", # get the items setting of Select extension
    "var idx = table[type + 's']({selected: true}).indexes().toArray();
    var values = [];
    for (var i = 0 ; i<idx.length; i++){
    values.push(table.table().column(0).data()[idx[i]])
    }", # get the index of selected items
    "var DT_id = table.table().container().parentNode.id;", # get the output id of DT
    "Shiny.onInputChange('modeltable' + '_selected', values);", # send the index to input$outputid_selected
    "})"
  ), server = T)
  
  
  onclick("startSCM", {
    
    mod <- input$modeltable_selected[1]
    all <- allmods2()
    path <- as.character(all[all$Number == mod, 2])
    relpath <- sub("/mod.ctl", "", path)
    InitiateSCM(path = relpath, control = input$ace, alltokens = alltokens, allmods = allmods)
  })
  
  
  
  
  scm_table <- modalDialog(fluidPage(div(column(9,
                                                 actionButton("openlocation3", "Open Model Location"),
                                                 actionButton("startSCM", "Start SCM"),
                                                 DT::dataTableOutput('modeltable3'),
                                                 
                                                 radioGroupButtons("selecteddir3", NULL, choices = c("All"), selected = "All")
                                                 
  )),
  div(column(3, tabsetPanel(id = "modelinfo3",
                            tabPanel("Plots3",
                                     bsCollapse(bsCollapsePanel("Covariate Model",
                                                                actionButton("ranparvscov3", "ETAs vs Covariates"),
                                                                actionButton("parvscov3", "Parameters vs Covariates"),
                                                                actionButton("parvspar3", "Parameters vs Parameters"),
                                                                
                                                                actionButton("covscatter3", "Covariate Scatter Plots")),
                                                bsCollapsePanel("Goodness of Fit",
                                                                actionButton("basicgof3", "Basic GOF"),
                                                                actionButton("indplots3", "Individual Plots"),
                                                                actionButton("dvpredipred3", "DV vs Pred, IPRED")),
                                                multiple = T, open = c("Covariate Model"))
                            ),
                            tabPanel("Parameters3",
                                     bsCollapse(bsCollapsePanel("Structural (THETA)",
                                                                tableOutput("thetas3")
                                     ),
                                     bsCollapsePanel("IIV (OMEGA)",
                                                     tableOutput("omegas3")),
                                     bsCollapsePanel("RUV (SIGMA)",
                                                     tableOutput("sigmas3")),
                                     multiple = T,
                                     open = c("Structural (THETA)", "IIV (OMEGA)", "RUV (SIGMA)"))),
                            tabPanel("Covariate Model",
                                     tableOutput("regress.parms3")),
                            tabPanel("PsN3",
                                     actionButton("runvpc3", "VPC")))))), easyClose = T)
  
  # ---------------------------------------------------------------------
  
  datafile <- reactive({
    x <- input$ace
    a <- strsplit(x, "\n")[[1]]
    datapath <- strsplit(a[startsWith(a, "$DATA")], "\\s+")[[1]][2]
    
    read.csv(datapath)
  })
  
  output$data <- DT::renderDataTable({
    datafile()
  }, filter = 'top', options = list(pageLength = 100,
                                  dom = 't', 
                                  autoWidth = TRUE,
                                  columnDefs = list(list(className = 'dt-center', 
                                                         targets = "_all"))), 
  rownames = FALSE)
  
  
  
  
  # Parameter tables --------------------------------------------------------
  
  output$thetas <- renderTable({
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    
    a <- readLines(paste0(homepath, "/", relpath, "/results/raw_results_structure"))[-1]
    
    b <- strsplit(a, "[=|,]")
    
    c <- b[which(lengths(b) == 3)]
    
    
    d <- do.call(cbind, c)
    
    final <- as.data.frame(d, stringsAsFactors = F)
    
    
    names(final) <- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2]) + 1,
                        as.numeric(final$theta[2]) + as.numeric(final$theta[3]))
    omegaindices <-  seq(as.numeric(final$omega[2]) + 1,
                         as.numeric(final$omega[2]) + as.numeric(final$omega[3]))
    sigmaindices <-  seq(as.numeric(final$sigma[2]) + 1,
                         as.numeric(final$sigma[2]) + as.numeric(final$sigma[3]))
    
    thetaSE <-  seq(as.numeric(final$setheta[2]) + 1, 
                    as.numeric(final$setheta[2]) + as.numeric(final$setheta[3]))
    omegaSE <-  seq(as.numeric(final$seomega[2]) + 1, 
                    as.numeric(final$seomega[2]) + as.numeric(final$seomega[3]))
    sigmaSE <-  seq(as.numeric(final$sesigma[2]) + 1,
                    as.numeric(final$sesigma[2]) + as.numeric(final$sesigma[3]))
    
    etaSHR <-  seq(as.numeric(final$shrinkage_eta[2]) + 1,
                   as.numeric(final$shrinkage_eta[2]) + as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath, "/", relpath, "/results/raw_results_mod.csv"))
    
    THETA <- rawResults[thetaindices] 
    
    THETASE <- abs(rawResults[thetaSE]/rawResults[thetaindices] * 100)
    
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] * 100
    ETASHR <- rawResults[etaSHR]
    
    SIGMA <- rawResults[sigmaindices]
    SIGMARSE <- rawResults[sigmaSE]/rawResults[sigmaindices] * 100
    
    THETATBL <- as.data.frame(cbind(t(THETA), t(THETASE)))
    
    names(THETATBL) <- c("Estimate", "RSE (%)")
    
    THETATBL
  }, rownames = T)
  
  output$omegas <- renderTable({
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    
    a <- readLines(paste0(homepath, "/", relpath, "/results/raw_results_structure"))[-1]
    b <- strsplit(a, "[=|,]")
    c <- b[which(lengths(b) == 3)]
    d <- do.call(cbind, c)
    final <- as.data.frame(d, stringsAsFactors = F)
    names(final) <- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2]) + 1,
                        as.numeric(final$theta[2]) + as.numeric(final$theta[3]))
    omegaindices <- seq(as.numeric(final$omega[2]) + 1,
                        as.numeric(final$omega[2]) + as.numeric(final$omega[3]))
    sigmaindices <- seq(as.numeric(final$sigma[2]) + 1,
                        as.numeric(final$sigma[2]) + as.numeric(final$sigma[3]))
    
    thetaSE <- seq(as.numeric(final$setheta[2]) + 1,
                    as.numeric(final$setheta[2]) + as.numeric(final$setheta[3]))
    omegaSE <- seq(as.numeric(final$seomega[2]) + 1,
                    as.numeric(final$seomega[2]) + as.numeric(final$seomega[3]))
    sigmaSE <- seq(as.numeric(final$sesigma[2]) + 1,
                    as.numeric(final$sesigma[2]) + as.numeric(final$sesigma[3]))
    
    etaSHR <- seq(as.numeric(final$shrinkage_eta[2]) + 1,
                   as.numeric(final$shrinkage_eta[2]) + as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath, "/", relpath, "/results/raw_results_mod.csv"))
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] * 100
    ETASHR <- rawResults[etaSHR]
    
    
    OMEGATBL <- as.data.frame(cbind(t(OMEGA), t(OMEGARSE), t(ETASHR)))
    names(OMEGATBL) <- c("Estimate", "RSE (%)", "Shrinkage (%)")
    
    OMEGATBL
  }, rownames = T)
  
  output$sigmas <- renderTable({
    homepath <-getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    
    a <- readLines(paste0(homepath, "/", relpath, "/results/raw_results_structure"))[-1]
    b <- strsplit(a, "[=|,]")
    c <- b[which(lengths(b) == 3)]
    d <- do.call(cbind, c)
    final<- as.data.frame(d, stringsAsFactors = F)
    names(final) <- final[1,]
    
    thetaindices <- seq(as.numeric(final$theta[2]) + 1,
                        as.numeric(final$theta[2]) + as.numeric(final$theta[3]))
    omegaindices <- seq(as.numeric(final$omega[2]) + 1,
                        as.numeric(final$omega[2]) + as.numeric(final$omega[3]))
    sigmaindices <- seq(as.numeric(final$sigma[2]) + 1,
                        as.numeric(final$sigma[2]) + as.numeric(final$sigma[3]))
    
    thetaSE <- seq(as.numeric(final$setheta[2]) + 1, 
                   as.numeric(final$setheta[2]) + as.numeric(final$setheta[3]))
    omegaSE <- seq(as.numeric(final$seomega[2]) + 1,
                   as.numeric(final$seomega[2]) + as.numeric(final$seomega[3]))
    sigmaSE <- seq(as.numeric(final$sesigma[2]) + 1,
                   as.numeric(final$sesigma[2]) + as.numeric(final$sesigma[3]))
    
    etaSHR <- seq(as.numeric(final$shrinkage_eta[2]) + 1,
                   as.numeric(final$shrinkage_eta[2]) + as.numeric(final$shrinkage_eta[3]))
    
    rawResults <- read.csv(paste0(homepath, "/", relpath, "/results/raw_results_mod.csv"))
    
    THETA <- rawResults[thetaindices] 
    THETASE <- abs(rawResults[thetaSE]/rawResults[thetaindices] * 100)
    
    OMEGA <- rawResults[omegaindices] 
    OMEGARSE <- rawResults[omegaSE]/rawResults[omegaindices] * 100
    ETASHR <- rawResults[etaSHR]
    
    SIGMA <- rawResults[sigmaindices]
    SIGMARSE <- rawResults[sigmaSE]/rawResults[sigmaindices] * 100
    
    THETATBL <- as.data.frame(cbind(t(THETA), t(THETASE)))
    names(THETATBL) <- c("Estimate", "RSE (%)")
    
    OMEGATBL <- as.data.frame(cbind(t(OMEGA), t(OMEGARSE), t(ETASHR)))
    names(OMEGATBL) <- c("Estimate","RSE (%)", "Shrinkage (%)")
    
    SIGMATBL <- as.data.frame(cbind(t(SIGMA), t(SIGMARSE)))
    names(SIGMATBL) <- c("Estimate", "RSE (%)")
    SIGMATBL
  },rownames = T)
  
  
  
  
  output$regress.parms <- renderTable({
    
    
    xpose.runno <- '1'
    homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
    all <- allmods2()
    
    path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
    relpath <- sub("/mod.ctl", "", path)
    working.directory <- paste0(homepath, '/', relpath)
    setwd(working.directory)
    xpdb<-xpose.data(xpose.runno)
    setwd(homepath)
    
    k=0
    all <- list()
    for (i in xpdb@Prefs@Xvardef$parms){
      y <- distinct(xpdb@Data,ID,.keep_all = T) %>% select_(i)
      for (j in xpdb@Prefs@Xvardef$covariates){
        x <- distinct(xpdb@Data,ID, .keep_all=T) %>% select_(j)
        if(!is.factor(x[, 1])){
          k <- k+1
          
          a <- rbind(linmod(center(x[, 1]), as.numeric(as.character(y[, 1])), xlab = names(x), ylab = names(y)),
                   expmod(center(x[, 1]), as.numeric(as.character(y[, 1])), xlab = names(x), ylab = names(y)),
                   powmod(center2(x[, 1]), as.numeric(as.character(y[, 1])), xlab = names(x), ylab = names(y)))
          all[[k]] <- a
        }
      }
    }
    
    allcor <- data.frame()
    allcor <- rbind(allcor, do.call(rbind, all))
    
    
    
    
  }, rownames = F)
  
  
  # Plot functions ----------------------------------------------------------
  
  #ETAs vs Covariates
  onclick("ranparvscov",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"\
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('etavscov.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(ranpar.vs.cov(xpdb))
            dev.off()
            shell("etavscov.pdf", wait = F)
            setwd(homepath)
          })
  
  onclick("parvscov",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"\
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('parvscov.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(parm.vs.cov(xpdb))
            dev.off()
            shell("parvscov.pdf", wait = F)
            setwd(homepath)
          })
  
  onclick("parvspar",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"\
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1],2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('parvspar.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(parm.vs.parm(xpdb))
            dev.off()
            shell("parvspar.pdf", wait = F)
            setwd(homepath)
          })
  
  #Basic Goodness of fit - (DV vs POPPRED, DV vs IPRED, |IWRES| vs DV, CWRES vs Time)
  onclick("basicgof",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('basicGOF.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(basic.gof(xpdb))
            dev.off()
            shell("basicGOF.pdf", wait = F)
            setwd(homepath)
          })
  
  onclick("indplots",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('indplots.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb<-xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(ind.plots(xpdb, logy = T))
            dev.off()
            shell("indplots.pdf", wait = F)
            setwd(homepath)
          })
  
  onclick("dvpredipred",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('dvpredipred.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(dv.vs.pred.ipred(xpdb))
            dev.off()
            shell("dvpredipred.pdf", wait = F)
            setwd(homepath)
          })
  
  onclick("covscatter",
          {
            xpose.runno <- '1'
            homepath <- getwd()#"M:/Users/mhismail-shared/Rprogramming/genetic algorithm/GA"
            all <- allmods2()
            
            path <- as.character(all[all$Number == input$modeltable_selected[1], 2])
            relpath <- sub("/mod.ctl", "", path)
            working.directory <- paste0(homepath, '/', relpath)
            setwd(working.directory)
            pdf.filename <- paste0('covscatter.pdf')
            pdf.title <- 'execute diagnostic plots run 1'
            xpdb <- xpose.data(xpose.runno)
            pdf(file = pdf.filename, width = 10, height = 7, title = pdf.title)
            print(cov.splom(xpdb))
            dev.off()
            shell("covscatter.pdf", wait = F)
            setwd(homepath)
          })
  }

#  ------------------------------------------------------------------------


# Run the application 
shinyApp(ui = ui, server = server)


