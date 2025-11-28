## 27 Nov 2019 / Nghia:
# add standard error for the estimates
## 04 June 2019 / Nghia:
# add parameter: isoform.method=average/total to report the expression of the individual members of a paralog i) average (default) or ii) total from the paralog set
## 01 Apr 2019 / Wenjiang:
# add "merge.paralogs" parameter to turn on/off the paralog merging in XAEM. The default is off, which will generate the same set of isoforms between different projects. To turn it on, just add "merge.paralogs=TRUE"
# Example of command: Rscript buildCRP.R in=eqClass.txt isoform.out=X_matrix.RData merge.paralogs=TRUE

## Take the workdir and core arguments
rm(list = ls())

design.matrix="X_matrix/X_matrix.RData"
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
    for (i in 1:length(args)){
        res=unlist(strsplit(args[i],"="))
        if (res[1]=="in"){pan_eqv_dir=as.character(res[2])}
        if (res[1]=="Xmatrix"){design.matrix=as.character(res[2])}
    }
} else {
    print("Error: input folder is missing!")
    q(save="no")
}

#curdir = getwd()
#workdir= paste0(curdir,"/", sample_Dir, "/eqv")
curdir = getwd()
workdir= paste0(curdir,"/", pan_eqv_dir)

core = 8 #default
merge.paralogs = TRUE ## default is to combine paralogs in the updated X to obtain the best performance
fout="XAEM_isoform_expression.RData"
foutr="XAEM_paralog_expression.RData"
#design.matrix="X_matrix/X_matrix.RData"
isoform.method="average" #  "average" or "total"
remove.ycount=TRUE

saveSubset=TRUE #save singleton and paralogs 
noBiasResult=FALSE #export the results without bias correction
foutr_noBias="XAEM_noBiasCor.RData"

source("/path/to/MegaFun/Rsource.R")
load(design.matrix)


options(stringsAsFactors=FALSE)
setwd(workdir)

#load input data
load("Ycount.RData")
#set parallel
library(foreach)
library(doParallel)
ncores = detectCores()
nc = min(ncores,core)     # use 8 or 16 as needed!!
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)

##### start from here
X.y= Y


fun = function(crpdat, maxiter.X=5, modify=TRUE){
  #xloc = grep('N', colnames(crpdat))
  xloc = which(colnames(crpdat) != "sample1")
  X0 = matrix(crpdat[,xloc], ncol=length(xloc))
  Ymat = crpdat[,-xloc]
  est = AEM(X0, Ymat, maxiter.X=maxiter.X, modify=modify)
  return(est)
}
X.y= Y

EST <- foreach(i= 1:length(X.y)) %dopar% fun(X.y[[i]], maxiter.X=5)

names(EST) = names(X.y)

#save(EST,file="Updated_X.Rdata")

####
#### Do not run add paralog step 
####
# if(!merge.paralogs)
# {
# x.all = list()
# X.y=Y
# for(i in 1:length(X.y))
# {
#  x1=x2=NULL
#  x1 = EST[[i]]$X
#  x.y =  X.y[[i]]
#  #xloc = grep('N', colnames(x.y))
#  xloc = which(colnames(x.y) != "sample1")
#  colnames(x1) = colnames(x.y)[xloc]
#  x.all[[i]] = x1
# }
# }
#
####
#### run add paralog step with X from EST result
####
if(merge.paralogs)
{
x.all = list()
X.y=Y
for(i in 1:length(X.y))
{
 x1=x2=NULL
 x1 = EST[[i]]$X
 x.y =  X.y[[i]]
 #xloc = grep('N', colnames(x.y))
 xloc = which(colnames(x.y) != "sample1")
 colnames(x1) = colnames(x.y)[xloc]
 x.all[[i]] = try(ccrpfun(x1),silent = TRUE)
}
run.err=NULL
for(i in 1:length(X.y))
 if(class(x.all[[i]])[1]== "try-error")
  run.err=c(run.err,i)

 for(i in run.err)
 {
  x1=x2=NULL
  x1 = EST[[i]]$X
  x.y =  X.y[[i]]
  #xloc = grep('N', colnames(x.y))
  xloc = which(colnames(x.y) != "sample1")
  colnames(x1) = colnames(x.y)[xloc]
  x.all[[i]] = ccrpfun(x1,clim=5) #clim shoule be smaller, such as 10
 }
}
#save(x.all, file="One_more_Collapsing_X.RData")

cat("\n...Estimation using AEM algorithm...\n")
X.y=Y
beta.all = list()# use the new X to calculate new beta
err.all = list()# keep standard error
for(i in 1:length(X.y))
{
 x2=NULL
 x.y =  X.y[[i]]
 #xloc = grep('N', colnames(x.y))
 xloc = which(colnames(x.y) != "sample1")
 y = x.y[,-xloc]
 x2 = x.all[[i]]
 beta1 = foreach(j=1:ncol(y)) %dopar% EM(x2,y[,j], maxiter = 50)
 beta2 = Reduce(rbind,beta1)
 rownames(beta2) = NULL
 beta.all = c(beta.all,list(beta2))
  err = colSums((x2 %*% t(beta2) - y)^2)
  err.all = c(err.all, list(err))
#  if (merge.paralogs){ #compute standard error only if using merge.paralogs=TRUE
   # s2=getSE(x2,beta2,y)
#    se.all = c(se.all,list(s2))   
#  }
}
if (saveSubset) save(beta.all, err.all,file="Beta_final_paralog.Rdata")

##### process singletons
## singletons
estfun = function(mat){
  CRP.y = mat$crpcount
  #sum(!(names(CRP.y)==names(CRP))) 
  ## cluster info
  npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
  table(npat)
  loc1 = which(npat==1)  # clusters with 1 pattern

  TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
   names(TC1)=names(CRP.y[loc1])
  ## single tx
  est.all =  c(TC1)
  return(est.all)
}

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
result_est=NULL
for(id in 1:length(flist)){ # call crpcount()
  load(flist[id])
  est1 = estfun(y)# estimation step
  result_est = cbind(result_est,est1)
  #cat("sample ",id,'\n')
}

#samplename1 = gsub(pattern = "pq", replacement = "N", x = samplename1)
colnames(result_est)=samplename1

if (saveSubset) save(result_est,file="Est_result_Singletons.RData")

seq1 = Reduce(cbind,beta.all);seq1=t(seq1)
XAEM_count = rbind(result_est,seq1)
save(XAEM_count, file = "XAEM_count.RData")



