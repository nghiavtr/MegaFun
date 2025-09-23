### create the map between the short name and the long name of the pangenomes
panFn="genome.fa"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="pangenome") panFn=as.character(res[2])
}

suppressMessages(suppressWarnings(library(Biostrings)))
pan_genome = readDNAStringSet(panFn)
pan_name=names(pan_genome)
pan_shortname = sapply(pan_name, function(s) {
  s = unlist(strsplit(s, "\\."))
  s = s[length(s)]
  s = unlist(strsplit(s, "\\|"))
  paste0(s[1], "|", s[2])
})

#create mapping from short-long names
pan_long2short=pan_shortname
names(pan_long2short)=pan_name
pan_short2long=pan_name
names(pan_short2long)=pan_shortname

save(pan_long2short,pan_short2long,file="panNameMap.RData")
