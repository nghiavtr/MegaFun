## Take the workdir and core arguments
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
workdir = args[1]
design.matrix = args[2]

# workdir= paste0(curdir,"/output/eqv")
# design.matrix="output/X_matrix/X_matrix.RData"
#core = 1 #default

source("/path/to/MegaFun/Rsource.R")

load(design.matrix)
options(stringsAsFactors=FALSE)
setwd(workdir)

if(!dir.exists("Ycount")) dir.create("Ycount")
setwd(paste(workdir,"/Ycount",sep=""))
flist = list.files(workdir,pattern="eqClass.txt",recursive=TRUE,full.names = TRUE)
tx_length = flist[1]
#initialization for parallel computing
#library(foreach)
#library(doParallel)
#registerDoParallel(cores=core)

#res=foreach(id = 1:length(flist),.combine=c) %dopar% {
for (id in 1:length(flist)){
  y = NULL
  y = crpcount(flist[id])
  samplename = paste(y$samplename,".RData",sep="")

  save(y,file=samplename)
#  return(flist[id])
}
res=NULL

flist = list.files(paste(workdir,"/Ycount",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
npat = sapply(CRP,nrow)   # number of occupancy patterns per clusterloc2 = which(npat>1) 
loc2 = which(npat>1) 

if(length(loc2) == 0){
  Y = NULL
  save(Y,file='../Ycount.RData')
  stop("No eqv in Y")
}

CCRP1 = CCRP[loc2]
Y=NULL
for(id in 1:length(flist)){
 cat("Merging results from sample ",flist[id],' ...\n')
 load(flist[id])
 if(id==1){
  Y = y[[1]][loc2]
 }
 if(id>1){
  y1 = y[[1]][loc2]
  for(i in 1:length(Y)){
   Y[[i]] = cbind(Y[[i]],sample1=y1[[i]][,'sample1'])
   }
  }
 }



samplename1 = NULL
for(id in 1:length(flist))
 {
  s.1 = strsplit(flist[id],"/")[[1]]
  s = s.1[length(s.1)]
  s = gsub(".RData","",s)
  samplename1 = c(samplename1,s)
 }

 setwd(workdir)

 for(i in 1:length(Y))
 {
 y2 = Y[[i]]
 xloc = which(colnames(y2) != "sample1")
 y2.1 = y2[,-xloc]
 y3 = cbind(CCRP1[[i]],y2.1)
 Y[[i]] = as.matrix(y3)
 }

save(Y,samplename1,tx_length,file='Ycount.RData')

cat("\n...Done...\n")
