# extract output of metaphlan to get the list of species
# Rscript extractSpeciesFromMetaphlan.R in=$metaphlan_out
# Rscript extractSpeciesFromMetaphlan.R in=$metaphlan_out p=0.5 out=$metaphlan_out
rm(list = ls())

metaphlan_outDir="metaphlan_out"
cutoff=0.5 # as the default threshold of Humann (0.5/100 = 0.5)
outdir=NULL

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="in") metaphlan_outDir=as.character(res[2])
    if (res[1]=="out") outdir=as.character(res[2])
    if (res[1]=="p") cutoff=as.double(res[2])
  }

if (is.null(outdir)) outdir=metaphlan_outDir

library(data.table)


flist = list.files(metaphlan_outDir,pattern="profile.txt",recursive=FALSE,full.names = TRUE)
#flist

tax_metaphlan_list=NULL
for (fn in flist){
  x=fread(fn)
  x=as.data.frame(x)
  p=grep("s__",x[,1])
  x1=x[p,]
  s=strsplit(x1[,1],"g__")
  s=do.call(rbind,s)[,2]
  s=paste0("g__",s)

  s1=strsplit(s,"\\|t__")
  s1=do.call(rbind,s1)[,1]
  x1$species=s1
  
  s2=gsub("\\|","\\.",s1)
  x1$speciesName=s2
  
  s3=sapply(s2, function(s) {
      temp <- unlist(strsplit(s, "__"))[3]
      gsub("_", " ", temp)
    })
  
  x1$speciesName0=s3
  
  
  x1$filename=basename(fn)

  tax_metaphlan_list=rbind(tax_metaphlan_list,x1)
}


p=which(tax_metaphlan_list$relative_abundance>cutoff) #cutoff at 2.5%
tax_metaphlan_filtered= tax_metaphlan_list[p,]

exportSpecies=unique(tax_metaphlan_filtered$speciesName)
#length(exportSpecies)
save(tax_metaphlan_list, tax_metaphlan_filtered, exportSpecies, file = paste0(outdir,"/tax_metaphlan_out.RData"))

exportSpecies_table=data.frame(id=exportSpecies)
write.table(exportSpecies_table, file=paste0(outdir,"/metaphlan_identified_species.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
