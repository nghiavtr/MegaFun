rm(list = ls())
# Rscript map_isolate_proteins_chunk.R in=isolate_genome_summary.txt.filtered.csv out=data chunkSize=1000 chunk=1

#vsearch to align isolate to pan to get protein id
#need to install vsearch before running
#very slow


isolate_df_fn="isolate_genome_summary.txt.filtered.csv"
outDir="data"
ncores=0
chunk=1
chunkSize=1000

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="in") isolate_df_fn=as.character(res[2])
    if (res[1]=="out") outDir=as.character(res[2])   
    if (res[1]=="chunk") chunkInd=as.integer(res[2])
    if (res[1]=="chunkSize") chunkSize=as.integer(res[2])
    if (res[1]=="t") ncores=as.integer(res[2])
    
  }


library(data.table)


if (!dir.exists(outDir)) dir.create(outDir)

isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)

#exclude the ones which were missing or could not be downloaded
isolate_df=isolate_df[which(isolate_df$missed==0),]

#######
fileNum=nrow(isolate_df)
#chunkInd=1
#chunkSize=1000

chunkStart=(chunkInd-1)*chunkSize+1
chunkEnd=chunkInd*chunkSize

if (chunkEnd > fileNum) chunkEnd=fileNum



# align isolate to pan to get protein id

for(i in 1:nrow(isolate_df))
  if (i >= chunkStart & i <= chunkEnd)
{
  if (i %% 100 == 0) cat ("\n Processing ", i)
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]


  cds_path = paste0(outDir,"/isolate_genome/CDS_", species, "@", strain, "@", biosample, ".fasta")
  pan_path = paste0(outDir,"/pangenome/PAN_", species, ".ffn")
  outFile=paste0(outDir,"/isolate_genome/ALIGN_CDS_PAN_", species, "@", strain, "@", biosample, ".txt")

  if (!file.exists(outFile)){
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

