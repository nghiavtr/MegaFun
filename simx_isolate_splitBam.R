### 2) extract reads from each bam file for individual CRPs
rm(list = ls())

core = 8
bamPath="" #must have
bamPattern="_sort.bam"
outPath="sam_temp/"
#crpPath="temp/"
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
  if (res[1]=="outpath") outPath=as.character(res[2])
  if (res[1]=="xmatrix") xmatrix=as.character(res[2])
  if (res[1]=="pangenome"){panfn=as.character(res[2])}
}


if(!dir.exists(outPath)) system(paste0("mkdir ",outPath),ignore.stdout = TRUE, ignore.stderr = TRUE)

suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(doParallel)))
cl <- makePSOCKcluster(core)   #
registerDoParallel(cl)


suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(data.table)))

set.seed(2023)

pan_genome = readDNAStringSet(panfn)
pan_genome_raw_names=names(pan_genome)
names(pan_genome) = sapply(pan_genome_raw_names, function(s) {
  s = unlist(strsplit(s, "\\."))
  s = s[length(s)]
  s = unlist(strsplit(s, "\\|"))
  paste0(s[1], "|", s[2])
})

#nghiavtr/09Feb2025: resolve the issues of matching between gene names from pan_genome and CRP names
#create mapping from raw to standard names
pan_genome_raw_names_map=names(pan_genome)
names(pan_genome_raw_names_map)=pan_genome_raw_names

load("X_matrix/X_matrix.RData")

nr = sapply(CRP, nrow)
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

### nghiavtr/07Nov2024 - end

bamList=list.files(bamPath,paste0(bamPattern,"$"))


system(paste0("mkdir ",outPath),ignore.stdout = TRUE, ignore.stderr = TRUE)

blockSize=1e6 # the size for reading one batch at once to process data
#scan each bam file then extract reads of the bam files to CRP, we will do it only once
#for(k in 1:5){
#res = foreach(k = 1:5, .packages = c("Biostrings","data.table")) %dopar% {
res = foreach(k = 1:length(bamList), .packages = c("Biostrings","data.table")) %dopar% {

### nghiavtr/07Nov2024 - start
#read individual sequencing reads and put them into the corresponding CRP
  #bamFn=paste0("simulated_data/stool_simulated_data_s", k, "_sort.bam")
  #bamFn="simulated_data/stool_simulated_data_s1_sort.bam"
  bamFn=paste0(bamPath,"/",bamList[k])

  cat("\n Sample in processing: ",bamFn)
  sampleName=basename(bamFn)
  temp_bam_path = paste0(outPath,"/", sampleName)
  if (!dir.exists(temp_bam_path)) system(paste0("mkdir ",temp_bam_path))

  samFn=paste0(temp_bam_path,"/",sampleName,".sam")  
  #cmd=paste0("samtools view -h ",bamFn, " > ", samFn) #convert all to sam
  cmd=paste0("samtools view -F 4 ",bamFn, " > ", samFn) #export only mapped reads before convert to sam
  system(cmd, ignore.stderr = TRUE)

  con <- file(samFn, "r", blocking = TRUE) #blocking = FALSE will make a random access to the file  
  x=readLines(con, n=1)
  while (substr(x,1,1)=="@") x=readLines(con, n=1) #skip the header

  count=0
  while(length(x)>0){ #if readLines() does not reach to the end of file
    cat("\n Number of processed reads: ",count)
    count=count+length(x)

    x1=strsplit(x,"\t")

#    t1=proc.time()
    x2=sapply(x1, function(y){y[3]},USE.NAMES = FALSE)
    x3=pan_genome_raw_names_map[x2]    
    x_crpID=p2cID[x3]
 #   proc.time()-t1  

    p=which(!is.na(x_crpID))

    if (length(p)>0){
      x_crpID=x_crpID[p]
      x1=x1[p]

#    t1=proc.time()
      x2=lapply(x1, function(y){return(c(y[10],y[1]))})
      x2=do.call(rbind,x2)
      x_fa=DNAStringSet(x2[,1])
      names(x_fa)=x2[,2]
#    proc.time()-t1


      x_fa_list=tapply(x_fa,x_crpID,c)

      for (i in 1:length(x_fa_list)){
        myCrpID=names(x_fa_list)[i]
        temp_crp_path = paste0(temp_bam_path,"/", myCrpID)
        if (!dir.exists(temp_crp_path)) system(paste0("mkdir ",temp_crp_path))
        writeXStringSet(x_fa_list[[i]], file = paste0(temp_crp_path, "/pan_temp.fasta"),append=TRUE)
      }
    }


    x=readLines(con, n=blockSize)
  }

  close(con)

  system(paste0("rm ",samFn))
  
### nghiavtr/07Nov2024 - end
 list(NA)
}


