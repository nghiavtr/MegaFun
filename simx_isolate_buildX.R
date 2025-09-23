### 3) build CRP at isolate level
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
}

system(paste0("mkdir ",crpPath),ignore.stdout = TRUE, ignore.stderr = TRUE)


suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(doParallel)))

cl <- makePSOCKcluster(core)
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

suppressMessages(suppressWarnings(library(data.table)))

set.seed(2023)

load(xmatrix)

nr = sapply(CRP, nrow)
names(nr)=NULL
pick = which(nr > 1)
CRP = CRP[pick]

### nghiavtr/07Nov2024 - end

###############
### extract fasta sequences from isolate genomes and distribute them into individual CRPs. This requires to scan all isolate genomes once.
### this substantially improve the speed in comparison with older versions

### now generate X-matrix
res = foreach(i = 1:length(CRP), .packages = c("Biostrings","data.table")) %dopar% {
#  print(i)
  tryCatch({

  temp_path = paste0(crpPath,"/", i)
  if(!dir.exists(temp_path)) dir.create(temp_path)


  system(paste0("mkdir ", temp_path, "/X_matrix"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("/path/to/MegaFun/bin/simx ", temp_path, "/isolate_seqs -o ", temp_path, "/X_matrix -t 1 -d 0.1 -m 50"), ignore.stdout = FALSE, ignore.stderr = FALSE)
  # system("bowtie2-build output/isolate_seqs/isolate_seqs.fasta output/bowtie2_index/isolate_seqs", ignore.stdout = T)
  # system("bowtie2 -f -x output/bowtie2_index/isolate_seqs -U output/X_matrix/isolate_seqs_sim.fasta -p 1 --very-sensitive -a | samtools sort -n -o output/X_matrix/sample_01.bam")
  system(paste0("minimap2 --secondary=yes -t 1 -Y -ax sr ", temp_path, "/isolate_seqs/isolate_seqs.fasta ", temp_path, "/X_matrix/isolate_seqs_sim.fasta.gz | samtools sort -n -o ", temp_path, "/X_matrix/sample_01.bam"), ignore.stdout = FALSE, ignore.stderr = FALSE)

  system(paste0("/path/to/MegaFun/bin/gentc_isolate ", temp_path, "/X_matrix/sample_01.bam -o ", temp_path, "/X_matrix/eqClass.txt -nm 5 -t 1"), ignore.stdout = FALSE, ignore.stderr = FALSE)
  eqc_path = paste0(temp_path, "/X_matrix/eqClass.txt")
  xmat_path = paste0(temp_path, "/X_matrix/X_matrix.RData")
  system(paste0("Rscript /path/to/MegaFun/buildCRP_isolate.R ", eqc_path, " ", xmat_path), ignore.stdout = FALSE, ignore.stderr = FALSE)

  #system(paste0("rm ",temp_path, "/X_matrix/sample_01.bam")) #remove bam file
  #system(paste0("rm ",temp_path, "/X_matrix/isolate_seqs_sim.fasta.gz")) #remove generated fasta file
  

  list(NULL, NULL)

  }, error = function(e){
    list(NULL, NULL)
  })# of trycatch
} #of foreach

cat("\n Build X-matrix at isolate level: done!")


