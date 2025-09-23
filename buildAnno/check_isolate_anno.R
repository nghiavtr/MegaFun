rm(list = ls())
# Rscript check_isolate_anno.R in=isolate_genome_summary.txt out=data chocoPath=/nfs/NGS/MegaFun_proj/full_chocophlan.v201901_v31 

#check if the download files are correct - redownload or clean

isolate_df_fn="isolate_genome_summary.txt"
chocophlan_path="/nfs/NGS/MegaFun_proj/full_chocophlan.v201901_v31"
outDir="data"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="in") isolate_df_fn=as.character(res[2])
    if (res[1]=="out") outDir=as.character(res[2])
    if (res[1]=="chocoPath") chocophlan_path=as.character(res[2])

    
  }


library(data.table)


isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)


# download genome and annotation for each isolate
for(i in 1:nrow(isolate_df)){
  if (i%%1000==0) print(i)
  link = isolate_df$ftp_path[i]
  id = unlist(strsplit(link, "/"))
  id = id[length(id)]
  fa_path = paste0(link, "/", id, "_genomic.fna.gz")
  gtf_path = paste0(link, "/", id, "_genomic.gtf.gz")
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]
  
  fn1= paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample , ".fna.gz")
  fn2= paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample , ".gtf.gz")
  cmd1=paste0("wget -O ",fn1, " ",fa_path, " ")
  cmd2=paste0("wget -O ",fn2, " ",gtf_path, " ")

  #redownload
  if (file.exists(fn1) & file.info(fn1)$size==0){
    cat("\n redownload: ",i, " ",fn1)
    system(cmd1)
  }
  if (file.exists(fn2) & file.info(fn2)$size==0){
    cat("\n redownload: ",i, " ",fn2)
    system(cmd2)
  }
  
}



# delete them if they are still error
for(i in 1:nrow(isolate_df)){
  if (i%%1000==0) print(i)
  link = isolate_df$ftp_path[i]
  id = unlist(strsplit(link, "/"))
  id = id[length(id)]
  fa_path = paste0(link, "/", id, "_genomic.fna.gz")
  gtf_path = paste0(link, "/", id, "_genomic.gtf.gz")
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]
  
  fn1= paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample , ".fna.gz")
  fn2= paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample , ".gtf.gz")
#  cmd1=paste0("wget -O ",fn1, " ",fa_path, " ")
#  cmd2=paste0("wget -O ",fn2, " ",gtf_path, " ")

  cmd1=paste0("rm ",fn1 )
  cmd2=paste0("rm ",fn2 )

  #delete  
  #cat("\n ",i, " ", fn1, " ", fn2)

  if (file.exists(fn1) & file.info(fn1)$size==0){
    cat("\n delete: ",i, " ",fn1)
    system(cmd1)
    if (file.exists(fn2)){
      cat("\n delete: ",i, " ",fn2)
      system(cmd2)
    }
  }

  if (file.exists(fn2) & file.info(fn2)$size==0){
    cat("\n delete: ",i, " ",fn2)
    system(cmd2)
    if (file.exists(fn1)){
      cat("\n delete: ",i, " ",fn1)
      system(cmd1)
    }
  }


  
}






