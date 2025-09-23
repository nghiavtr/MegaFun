rm(list = ls())
# Rscript get_pan_isolate_anno.R in=isolate_genome_summary.txt out=data chocoPath=/nfs/NGS/MegaFun_proj/full_chocophlan.v201901_v31 

#extract pangenome sequence and download isolate sequence

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

if (!dir.exists(outDir)) dir.create(outDir)
dir.create(paste0(outDir,"/isolate_genome"))
dir.create(paste0(outDir,"/pangenome"))


isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)

### extract pangenome sequence
# need to download Chocophlan and change the path

isolate_df_unique=isolate_df[!duplicated(isolate_df$species_id),]

for(i in 1:nrow(isolate_df_unique)){
  species = isolate_df_unique$species_id[i]
  species_choco_filename = isolate_df_unique$species_choco_filename[i]
  # cp pangenomes
  cmd = paste0("cp ",chocophlan_path,"/", species_choco_filename, " ", outDir,"/pangenome/PAN_", species, ".ffn.gz")
  system(cmd)
  system(paste0("gzip -d ",outDir,"/pangenome/PAN_", species, ".ffn.gz")) #gunzip and change name to the standard
}


############# download isolate sequence

# download genome and annotation for each isolate
for(i in 1:nrow(isolate_df)){
  if (i%%100==0) print(i)
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

  #continue to download if it was interupted - so dont have to re-run
  if (!file.exists(fn1) | file.info(fn1)$size==0){
    cat("\n download: ",i, " ",fn1)
    system(cmd1)
  }
  if (!file.exists(fn2) | file.info(fn2)$size==0){
    cat("\n download: ",i, " ",fn2)
    system(cmd2)
  }
  
}




