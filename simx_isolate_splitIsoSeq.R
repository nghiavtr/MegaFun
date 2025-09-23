### 1) Split isolate sequences
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
  if (res[1]=="iso"){iso_Dir=as.character(res[2])}
  if (res[1]=="pangenome"){panfn=as.character(res[2])}
  if (res[1]=="smap"){smapfn=as.character(res[2])}
  if (res[1]=="core") core=as.integer(res[2])  
  if (res[1]=="crppath") crpPath=as.character(res[2])
  if (res[1]=="xmatrix") xmatrix=as.character(res[2])
}

if(!dir.exists(crpPath)) system(paste0("mkdir ",crpPath),ignore.stdout = TRUE, ignore.stderr = TRUE)

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
### extract fasta sequences from isolate genomes and distribute them into individual CRPs. This requires to scan all isolate genomes once.
### this substantially improve the speed in comparison with older versions

smap=fread(smapfn)

flist=list.files(iso_Dir,"*.fasta")
i=1
postNameList=list()
#for (i in 1:length(flist)){
res = foreach(i = 1:length(flist), .packages = c("Biostrings","data.table")) %dopar% {

  sname=strsplit(basename(flist[i]),"\\.")[[1]][1]
  cat("\nProcessing species ",i, ":", sname)
  fa_in=readDNAStringSet(paste0(iso_Dir,"/", flist[i]))

  fa_name=names(fa_in)
  
  ### nghiavtr/10May2025 - fix the bugs of species name 
  #use correct name - NOTE: the species name in fa_name (ncbi name) and sname by the flist name (by chocophlan name) can be different
  choco_sname=sname
  if (!choco_sname %in% smap$choco_sname) choco_sname=smap$choco_sname[smap$ncbi_sname==sname]
  fa_UniRef90name = sapply(fa_name, function(s) {
    s = unlist(strsplit(s, "\\."))
    s = s[length(s)]
    s = unlist(strsplit(s, "\\|"))
    #paste0(s[1], "|", s[2])
    paste0(choco_sname, "|", s[2])
    #return(s[2])
  })
  names(fa_UniRef90name)=NULL

  #get protein name
  #fa_species_name=paste0(sname,"|",fa_name)
  fa_species_name=fa_UniRef90name

  keep=which(fa_species_name %in% names(p2cID))
  fa_in=fa_in[keep]
  fa_species_name=fa_species_name[keep]

  #get crpID
  fa_crpID=p2cID[fa_species_name]
  #revise names
  x=tapply(c(1:length(fa_species_name)),fa_species_name,c)
  fa_species_isoname=fa_species_name
  for (j in 1:length(x)){
    x1=x[[j]]
    fa_species_isoname[x1]=paste0(fa_species_isoname[x1],"@",1:length(x1))
  }
  names(fa_in)=fa_species_isoname
  fa_out=tapply(fa_in,fa_crpID,c)
  fa_out_species_name=tapply(fa_species_name,fa_crpID,unique)

  for (j in 1:length(fa_out)){
    temp_path = paste0(crpPath,"/", names(fa_out)[j])
    if(!dir.exists(temp_path)) dir.create(temp_path)
    system(paste0("mkdir ", temp_path, "/isolate_seqs"), ignore.stdout = TRUE, ignore.stderr = TRUE)

    writeXStringSet(fa_out[[j]], file = paste0(temp_path,"/isolate_seqs", "/isolate_seqs.fasta"),append=TRUE)
    write.table(fa_out_species_name[[j]], file = paste0(temp_path, "/temp.list"), row.names = F, col.names = F, quote = F,append=TRUE)
  }

  #postNameList[[sname]]=sort(postName)
  return(sname)
}
res=NULL

###### 
cat("\n split isolate sequence: done!")
