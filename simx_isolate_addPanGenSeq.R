### 2) add pan-genome sequences
rm(list = ls())

core = 8
crpPath="temp/"
xmatrix="X_matrix/X_matrix.RData"



args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="core") core=as.integer(res[2])  
  if (res[1]=="crppath") crpPath=as.character(res[2])
  if (res[1]=="xmatrix") xmatrix=as.character(res[2])
  if (res[1]=="pangenome"){panfn=as.character(res[2])}
}

system(paste0("mkdir ",crpPath),ignore.stdout = TRUE, ignore.stderr = TRUE)


suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(doParallel)))

cl <- makePSOCKcluster(core)   #
registerDoParallel(cl)

estfun = function(mat){
  CRP.y = mat$crpcount
  ## cluster info
  npat = sapply(CRP.y,nrow)   # number of occupancy patterns per cluster
  loc1 = which(npat==1)  # clusters with 1 pattern

  est.all = NULL
  if(length(loc1) > 0){
    TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
     names(TC1)=names(CRP.y[loc1])
    ## single tx
    est.all =  c(TC1)
  }
  return(est.all)
}

suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(data.table)))

set.seed(2023)

#load("postNameList.RData")

pan_genome = readDNAStringSet(panfn)
names(pan_genome) = sapply(names(pan_genome), function(s) {
  s = unlist(strsplit(s, "\\."))
  s = s[length(s)]
  s = unlist(strsplit(s, "\\|"))
  paste0(s[1], "|", s[2])
})

load(xmatrix)

nr = sapply(CRP, nrow)
names(nr)=NULL
pick = which(nr > 1)
CRP = CRP[pick]


### nghiavtr/07Nov2024 - start
# mapping between CRP and tx

#map from crp to tx
c2t=sapply(names(CRP), strsplit, split=' ')#crp to tx  
cLen=lengths(c2t)
t2c=rep(names(cLen),cLen) #tx to crp
tVec=unlist(c2t)  
names(t2c) = tVec
crpID=seq(1:length(CRP)) #get CRP indicies
names(crpID)=names(CRP)

x=tVec
names(x)=NULL
x1 = sapply(x, function(s) {
  s = unlist(strsplit(s, "\\."))
  s = s[length(s)]
  s = unlist(strsplit(s, "\\|"))
  paste0(s[1], "|", s[2])
})

p2t=tapply(names(x1),x1,c) #pan_genome names to tx, it is 1-1 mapping
x1_cname=t2c[names(x1)]
p2c=tapply(x1_cname,x1,c) #pan_genome names to CRP names, it is 1-1 mapping
x1_cID=crpID[x1_cname]
p2cID=tapply(x1_cID,x1,c) #pan_genome names to crpID, it is 1-1 mapping
cID2p=tapply(names(p2cID),p2cID,c) #crpID to pan_genome names , it is 1-n

### nghiavtr/07Nov2024 - end

###############
### now add pan genome data if necessary, then build X-matrix

res = foreach(i = 1:length(CRP), .packages = c("Biostrings","data.table")) %dopar% {
  temp_path = paste0(crpPath,"/", i)
  if(!dir.exists(temp_path)) dir.create(temp_path)
  
  panGenelist=cID2p[[as.character(i)]]
  
  iso_exist_file=paste0(temp_path, "/temp.list")
  if (file.exists(iso_exist_file)){
    iso_exist = read.table(iso_exist_file,header=FALSE)
    pick=which(! (panGenelist %in% iso_exist[,1]) )
  }else{
    pick=c(1:length(panGenelist))
  }

  pan_seqs = NULL
  #pick = which(len_seqs == 0) #if can not found anything from isolate genome, then use the data from pan genome
  if(length(pick) > 0){
    #pick = panGenelist[pick]
    #pick = paste0(pick$species, "|", pick$uniref90)
    pan_seqs = pan_genome[panGenelist[pick]]
    names(pan_seqs) = paste0(names(pan_seqs), "@1")
    if(length(pan_seqs) != length(pick)){ #nghiavtr: not sure why we need this
      next
    }
  }

  #now add the sequence from pan genome if avail
  system(paste0("mkdir ", temp_path, "/isolate_seqs"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (length(pan_seqs) > 0) writeXStringSet(pan_seqs, file = paste0(temp_path, "/isolate_seqs/isolate_seqs.fasta"),append=TRUE)

  return(i)
} #of foreach


cat("\n add pan-genome sequence: done!")
