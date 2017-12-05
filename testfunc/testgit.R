library(git2r)

setwd("../../gittest")
path <-getwd()
repo<-init(".")



for (i in 2:10000){ 
  writeLines(paste0("Hello world!",i), file.path(path, paste0("example.txt")))
  add(repo, "example.txt")
  commit_1 <- commit(repo, paste0("commit message",i))
  branch_create(commit_1,paste0("mod",i))
  print(i)
}

branches(repo)
?branch_create
