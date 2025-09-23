rm(list = ls())
# Rscript get_isolate_genome.R in=isolate_genome_summary.txt.filtered.csv out=data
#extract aligned gene sequence from isolate for MegaFun
#need to install vsearch before running

#nghiavtr/30June2025: exclude duplicates, speed up with parallel processing

isolate_df_fn="isolate_genome_summary.txt.filtered.csv"
outDir="data"
ncores=NULL

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="in") isolate_df_fn=as.character(res[2])
    if (res[1]=="out") outDir=as.character(res[2])
    if (res[1]=="t") ncores=as.integer(res[2])
    
  }



library(BSgenome)
library(parallel)
library(Rsamtools)
library(data.table)
library(Biostrings) #stringDist
library(foreach)
library(doParallel)

if (!is.null(ncores)) nc = ncores else nc = detectCores()
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)


if (!dir.exists(outDir)) dir.create(outDir)
dir.create(paste0(outDir,"/isolate_genome_MegaFun"))


isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)

#exclude the ones which were missing or could not be downloaded
isolate_df=isolate_df[which(isolate_df$missed==0),]


myspecices_list=unique(isolate_df$species_id)
#species=myspecices_list[1]
# extract alignment results from vsearch
# assign ID = UNKNOWN for a isolate gene without matching in pangenome level
for (species in myspecices_list){
  sub_isolate_df=isolate_df[isolate_df$species_id==species,]
  
  fn_sps=paste0(outDir,"/isolate_genome_MegaFun/s__", sub_isolate_df$species_choco_id[1],".fasta")
  if (file.exists(fn_sps)) next()

  dir_sps=paste0(outDir,"/isolate_genome_MegaFun/", species)
  if (!dir.exists(dir_sps)) dir.create(dir_sps)
  system(paste0("rm ",dir_sps,"/*"))


  sps_fasta=NULL
  y=NULL
  cat("\n processing ",species)

  k=0
#  for(i in 1:nrow(sub_isolate_df)){
  res=foreach(i = 1:nrow(sub_isolate_df))%dopar%{
    library(data.table)
    library(Biostrings) #stringDist

    cat("\n iso ",i)
    #species = isolate_df$species_id[i]
    strain = sub_isolate_df$infraspecific_name[i]
    biosample = sub_isolate_df$biosample[i]

    alignFn=paste0(outDir,"/isolate_genome/ALIGN_CDS_PAN_", species, "@", strain, "@", biosample, ".txt")
    
    if (file.exists(alignFn)){

      cmd = paste0("cp ",alignFn," temp",i,".txt")
      #cmd = paste0("cp ",outDir,"/isolate_genome/ALIGN_CDS_PAN_", species, "@", strain, "@", biosample, ".txt temp",i,".txt")
      system(cmd)
      cmd1=paste0('grep "Query >" temp',i,'.txt -A 2 > temp',i)
      system(cmd1)
      df = readLines(paste0("temp",i))
      n = (length(df) + 1)/4

      anno = lapply(0:(n-1), function(j) {
        seq_id = gsub("Query >", "", df[j*4 + 1])
        uniref_id = unlist(strsplit(df[j*4 + 3], "\\|"))[c(3,4)]
        pan_len = unlist(strsplit(df[j*4 + 3], "\\|"))[5]
        c(seq_id, uniref_id, pan_len)
      })
      anno = do.call(rbind, anno)
      anno = as.data.frame(anno)
      colnames(anno) = c("seq_id", "uniref90_id", "uniref50_id", "pan_len")
      anno$pan_len = as.numeric(anno$pan_len)

      #clean
      system(paste0("rm temp",i))
      system(paste0("rm temp",i,".txt"))

      cds_path = paste0(outDir,"/isolate_genome/CDS_", species, "@", strain, "@", biosample, ".fasta")
      if (file.exists(cds_path)){
        fasta = readDNAStringSet(cds_path)

        pick = match(anno$seq_id, names(fasta))
        anno$isolate_len = width(fasta)[pick]
        anno$diff_len = abs(anno$pan_len - anno$isolate_len)

        # remove the match if the difference between pan
        # and isolate lens is bigger than 3bp
        pick = which(anno$diff_len > 3)
        if(length(pick) > 0){
          anno = anno[-pick, ]
        }

        pick = match(names(fasta), anno$seq_id)
        pick_na = which(is.na(pick))

        names(fasta) = paste0(names(fasta), "|", anno$uniref90_id[pick], "|", anno$uniref50[pick])
        names(fasta)[pick_na] = gsub("\\|NA", "", names(fasta)[pick_na])
        names(fasta)[pick_na] = paste0(names(fasta)[pick_na], "|UNKNOWN|UNKNOWN")
        #writeXStringSet(fasta, file = paste0(outDir,"/isolate_genome_MegaFun/CDS_UNIREF_", species, "@", strain, "@", biosample, ".fasta")) # missing append=true

        x=as.character(fasta)
        p1=which(!duplicated(x))
        fasta=fasta[p1]

        out_k=paste0(dir_sps,"/",species,i)
        writeXStringSet(fasta, file = out_k, append=FALSE)

      }
    }
    return(NULL)
  }

  #merging using bi-tree approach
  allFiles=list.files(dir_sps,full.names = TRUE)
  r=0
  while (length(allFiles) > 30){ #no less than 30
    r=r+1
    cat("\n round ",r, " ",length(allFiles))
    res=foreach(i = 1:trunc(length(allFiles)/2))%dopar%{
      library(data.table)
      library(Biostrings) #stringDist

    #for (i in 1:trunc(length(allFiles)/2)){
    #  cat(i," ")
      i2=length(allFiles)-i+1
      f1 = readDNAStringSet(allFiles[i])
      f2 = readDNAStringSet(allFiles[i2])
      x1=as.character(f1)
      x2=as.character(f2)
      p=x2%in%x1
      p=which(!p)
      if (length(p)>0){
        f2=f2[p]
        f1=c(f1,f2)
      }
      writeXStringSet(f1, file = allFiles[i], append=FALSE)
      system(paste0("rm ",allFiles[i2]))
      return(NULL)
    }
    allFiles=list.files(dir_sps,full.names = TRUE)
    #length(allFiles)

  }
  # post process to avoid memory issues

  cat("\n slow down to avoid memory issues")
  allFiles=list.files(dir_sps,full.names = TRUE)
  cat("\n remaining: ",length(allFiles))

  if (length(allFiles) > 1){ #at least 2 files or more
    for (i in 1:(length(allFiles)-1)){
      cat(" ",i)
      f1 = readDNAStringSet(allFiles[i])
      x1=as.character(f1)    

      for (j in (i+1):length(allFiles)){
          f2 = readDNAStringSet(allFiles[j])
          x2=as.character(f2)
          p=x2%in%x1
          p=which(!p)
          if (length(p)>0){ #update f2
            f2=f2[p]
            writeXStringSet(f2, file = allFiles[j], append=FALSE)
          }
      }
    }
  }

  allFiles=list.files(dir_sps,full.names = TRUE)
  for (i in 1:length(allFiles)){
    cmd=paste0("cat ",allFiles[i], " >> ",fn_sps) #combine all remaining files together
    system(cmd)
  }
  system(paste0("rm -r ",dir_sps))
  #writeXStringSet(sps_fasta, file = fn_sps, append=TRUE)
  #writeXStringSet(sps_fasta, file = fn_sps, append=FALSE)
}
