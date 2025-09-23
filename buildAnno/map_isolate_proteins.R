rm(list = ls())
# Rscript map_isolate_proteins.R in=isolate_genome_summary.txt.filtered.csv out=data

#vsearch to align isolate to pan to get protein id
#need to install vsearch before running
#very slow


isolate_df_fn="isolate_genome_summary.txt.filtered.csv"
outDir="data"
ncores=0

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


library(data.table)


if (!dir.exists(outDir)) dir.create(outDir)

isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)

#exclude the ones which were missing or could not be downloaded
isolate_df=isolate_df[which(isolate_df$missed==0),]



# align isolate to pan to get protein id

for(i in 1:nrow(isolate_df)){
  if (i %% 100 == 0) cat ("\n Processing ", i)
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]


  cds_path = paste0(outDir,"/isolate_genome/CDS_", species, "@", strain, "@", biosample, ".fasta")
  pan_path = paste0(outDir,"/pangenome/PAN_", species, ".ffn")
  outFile=paste0(outDir,"/isolate_genome/ALIGN_CDS_PAN_", species, "@", strain, "@", biosample, ".txt")

  mycondition=!file.exists(outFile)
  if (file.exists(outFile) & file.info(outFile)$size==0) mycondition=TRUE
  
  if (mycondition){
    if (ncores>0){
      cmd = paste0("vsearch --usearch_global ", cds_path, " --db ", pan_path, " --threads ",ncores," --id 0.90 --alnout ",
                   outFile)
    }else{
    # dont use --threads, then it will take all resourcea
      cmd = paste0("vsearch --usearch_global ", cds_path, " --db ", pan_path, " --id 0.90 --alnout ",
                   outFile)
    }
    system(cmd)
  }
}

