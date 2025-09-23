
#Rscript cut_genename_pan.R panIn=pangenome panOut=pangenome_shortNames

panIn="pangenome"
panOut="pangenome_shortNames"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="panIn") panIn=as.character(res[2])
  if (res[1]=="panOut") panOut=as.character(res[2])
}

if(!dir.exists(panOut)) dir.create(panOut)

suppressMessages(suppressWarnings(library(Biostrings)))

flist=list.files(panIn,"*.ffn", full.names = TRUE)


panfn="PAN_g__Escherichia.s__Escherichia_coli.ffn"

for (panfn in flist){

  pan_genome = readDNAStringSet(panfn)
  pan_genome_raw_names=names(pan_genome)

  pan_genome_short_names = sapply(pan_genome_raw_names, function(s) {
    s = unlist(strsplit(s, "\\."))
    x = unlist(strsplit(s, "\\|"))
    l=length(s)
    s1=s[c(l-1,l)]    
    s2=paste0("PAN.",s1[1], ".", s1[2])
    s3=paste0(x[1],"|",s2)
    return(s3)
  })
  names(pan_genome_short_names)=NULL
  names(pan_genome)=pan_genome_short_names

  #export to files
  fout=paste0(panOut, "/",basename(panfn))
  writeXStringSet(pan_genome, file = fout)

}


