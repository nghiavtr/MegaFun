### 3) create y count matrix
rm(list = ls())

core = 8
bamPath="" #must have
bamPattern="_sort.bam"
outPath="sam_temp/"
crpPath="temp/"
xmatrix="X_matrix/X_matrix.RData"



args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="in") bamPath=as.character(res[2])
  if (res[1]=="pattern") bamPattern=as.character(res[2])
  if (res[1]=="core") core=as.integer(res[2])  
  if (res[1]=="crppath") crpPath=as.character(res[2])
  if (res[1]=="outpath") outPath=as.character(res[2])
  if (res[1]=="xmatrix") xmatrix=as.character(res[2])
}



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


load(xmatrix)

nr = sapply(CRP, nrow)
pick = which(nr > 1)
CRP = CRP[pick]

#system("mkdir temp")

bamList=list.files(bamPath,paste0(bamPattern,"$"))

temp_data_path = outPath
system(paste0("mkdir ",temp_data_path,"/", "eqv"),ignore.stdout = TRUE, ignore.stderr = TRUE)
i=2500
res = foreach(i = 1:length(CRP), .packages = c("Biostrings","data.table")) %dopar% {
# for(i in 1:30){
  print(i)
  tryCatch({

  temp_path = paste0(crpPath,"/", i)
  eqv_path = paste0(temp_data_path,"/", "eqv/",i)
  system(paste0("mkdir ", eqv_path), ignore.stdout = TRUE, ignore.stderr = TRUE)

  x = CRP[[i]]

  eqc_path = paste0(temp_path, "/X_matrix/eqClass.txt")
  xmat_path = paste0(temp_path, "/X_matrix/X_matrix.RData")
  
  for(k in 1:length(bamList)){

    bamFn=paste0(bamPath,"/",bamList[k])
    sampleName=basename(bamFn)
    temp_bam_path = paste0(temp_data_path,"/", sampleName)
    myCrpID=i
    temp_crp_path = paste0(temp_bam_path,"/", myCrpID)

    eqv_sample_path = paste0(eqv_path,"/", sampleName)
    system(paste0("mkdir ", eqv_sample_path), ignore.stdout = TRUE, ignore.stderr = TRUE)


    if (!dir.exists(temp_crp_path)){ #if not existing, create a null fasta file
      system(paste0("mkdir ",temp_crp_path))
      writeXStringSet(DNAStringSet(NULL), file = paste0(temp_crp_path, "/pan_temp.fasta"))
    }
  
    # system("bowtie2 -f -x output/bowtie2_index/isolate_seqs -U pan_temp.fasta -p 1 --local -a | samtools sort -n -o output/stool_simulated_data_s1.bam")
    system(paste0("minimap2 --secondary=yes -t 1 -Y -ax sr ", temp_path, "/isolate_seqs/isolate_seqs.fasta ", temp_crp_path, "/pan_temp.fasta | samtools sort -n -o ", temp_crp_path, "/pan_temp.bam"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("/path/to/MegaFun/bin/bam2eqv_isolate ", temp_crp_path, "/pan_temp.bam -o ", eqv_sample_path, "/eqClass.txt -nm 5 -t 1"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("rm ",temp_crp_path, "/pan_temp.bam")) #remove the generated bam file
  }

  #eqv_path = paste0(getwd(), "/", temp_path, "/eqv") #created in advance
  
#  workdir=paste0(getwd(), "/",eqv_path)
#  design.matrix=xmat_path
  system(paste0("Rscript /path/to/MegaFun/create_count_matrix_isolate.R ", getwd(), "/",eqv_path, " ", xmat_path), ignore.stdout = TRUE, ignore.stderr = TRUE)

  # singleton
  flist = list.files(paste0(eqv_path, "/Ycount"),pattern="RData",recursive=TRUE,full.names = TRUE)
  singleton_est=NULL
  for(id in 1:length(flist)){ # call crpcount()
    load(flist[id])
    est1 = estfun(y)# estimation step
    if(!is.null(est1)){
      singleton_est = cbind(singleton_est,est1)
    }
    #cat("sample ",id,'\n')
  }

  Y_final = NULL
  if(file.exists(paste0(eqv_path, "/Ycount.RData"))){
    load(paste0(eqv_path, "/Ycount.RData"))
    Y_final = Y
  }

  system(paste0("rm ", eqv_path, " -f -r"), ignore.stdout = TRUE, ignore.stderr = TRUE)

  list(Y_final, singleton_est)
  }, error = function(e){
    list(NULL, NULL)
  })
}

new_Y = lapply(res, function(r) r[[1]])
singleton = lapply(res, function(r) r[[2]])
# save(new_Y, new_CRP, new_CCRP, file = "output/new_XY.RData")
save(new_Y, file = "output/isolate_XY.RData")
save(singleton, file = "output/isolate_singleton.RData")
#system("rm temp -f -r")
