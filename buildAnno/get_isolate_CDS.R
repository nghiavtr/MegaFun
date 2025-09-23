rm(list = ls())
# Rscript get_isolate_CDS.R in=isolate_genome_summary.txt out=data

#get CDS of isolate



isolate_df_fn="isolate_genome_summary.txt"
outDir="data"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="in") isolate_df_fn=as.character(res[2])
    if (res[1]=="out") outDir=as.character(res[2])
    
  }



library(BSgenome)
library(parallel)
library(Rsamtools)
library(data.table)
library(Biostrings) #stringDist

if (!dir.exists(outDir)) dir.create(outDir)


isolate_df=fread(paste0(isolate_df_fn),sep="\t")
isolate_df=as.data.frame(isolate_df)


############# process to extract isolate sequences

#system(paste0("gzip -d ",outDir,"/isolate_genome/*.gz"))


# remove missing data
isolate_df$missed=0
for(i in 1:nrow(isolate_df)){
#  print(i)
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]

  fa_path = paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample, ".fna.gz")
  gtf_path = paste0(outDir,"/isolate_genome/", species, "@", strain,"@", biosample, ".gtf.gz")

  if (!file.exists(fa_path) | !file.exists(gtf_path)){
    isolate_df$missed[i]=1
  }
}

fwrite(isolate_df,file=paste0(isolate_df_fn,".filtered.csv"),sep="\t")

isolate_df=isolate_df[which(isolate_df$missed==0),]

# extract gene sequence
for(i in 1:nrow(isolate_df)){
  if (i %% 100 == 0) cat ("\n Processing ", i)
  species = isolate_df$species_id[i]
  strain = isolate_df$infraspecific_name[i]
  biosample = isolate_df$biosample[i]
  outFile=paste0(outDir,"/isolate_genome/CDS_", species, "@", strain, "@", biosample, ".fasta")
  
  mycondition=!file.exists(outFile)
  if (file.exists(outFile) & file.info(outFile)$size==0) mycondition=TRUE

  if (mycondition){

    fa_path = paste0(outDir,"/isolate_genome/", species, "@", strain, "@", biosample, ".fna.gz")
    gtf_path = paste0(outDir,"/isolate_genome/", species, "@", strain,"@", biosample, ".gtf.gz")

    fasta <- readDNAStringSet(fa_path)

    gtf <- rtracklayer::import(gtf_path)

    r <- gtf[which(elementMetadata(gtf)[, "type"] == "gene"), ]

    names(fasta) = sapply(names(fasta), function(s) {
      unlist(strsplit(s, " "))[1]
    })

    gene_seq <- getSeq(fasta, r)

    names(gene_seq) <- paste0(elementMetadata(r)[, "gene_id"], "@", species)

    #outFile=paste0(outDir,"/isolate_genome/CDS_", species, "@", strain, "@", biosample, ".fasta")
    writeXStringSet(gene_seq, file = outFile)
  }
}


