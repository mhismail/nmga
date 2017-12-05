a<- paste(readLines("../test/1.ctl"),collapse="\n")
b<- gsub("\\$DATA\\s*", "\\$DATA ../",a)
b
