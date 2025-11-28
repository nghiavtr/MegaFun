eqClassFn="X_matrix/eqClass.txt"
outFn=NULL

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="eqClass") eqClassFn=as.character(res[2])
  if (res[1]=="out") outFn=as.character(res[2])
}

if (is.null(outFn)) {
  outFn=paste0(eqClassFn,"_PAN")
}

library(data.table)

rawmat  = fread(eqClassFn,sep='\t',header=TRUE)
rawmat  = as.data.frame(rawmat)

#write.table(paste0(eqClassFn,"_backup"),sep='\t',header=TRUE, quote=FALSE, row.names=FALSE)


pan_name=rawmat$Transcript
pan_shortname = sapply(pan_name, function(s) {
  s1 = unlist(strsplit(s, "\\."))
  s1[1]=gsub("PAN","",s1[1])
  paste0(s1[1], s1[3])
})
names(pan_shortname)=NULL

rawmat$Transcript=pan_shortname
#write.table(eqClassFn,sep='\t',header=TRUE, quote=FALSE, row.names=FALSE)
write.table(rawmat,file=outFn,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
