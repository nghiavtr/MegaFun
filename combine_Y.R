load("output/isolate_XY.RData")
load("output/isolate_singleton.RData")
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
    for (i in 1:length(args)){
        res=unlist(strsplit(args[i],"="))
        if (res[1]=="in"){pan_eqv_dir=as.character(res[2])}
        if (res[1]=="pangenome"){panfn=as.character(res[2])}
    }
} else {
    print("Error: input folder is missing!")
    q(save="no")
}
#sample_Dir="new_sampleName"

pick = sapply(1:length(new_Y), function(i) is.null(new_Y[[i]]) & is.null(singleton[[i]]))
pick = which(pick == TRUE)

final_Y = new_Y


if(length(pick) > 0){
  load(paste0(pan_eqv_dir, "/Ycount.RData"))
  old_Y = Y[pick]
  final_Y = c(final_Y, old_Y)
}

#final_Y = do.call(c, final_Y)

len=lengths(final_Y)
Y1=final_Y[which(len==1)]
p=which(len>1)
for (i in 1:length(p)){
  j=p[i]
  Y1=c(Y1,final_Y[j][[1]])
}
final_Y=Y1

#process singleton
isolate_singleton = do.call(c, singleton)
isolate_singleton = do.call(rbind, singleton)

#get previous singletons estimated from pangenome
load(paste0(pan_eqv_dir, "/Est_result_Singletons.RData"))
#load("simulated_data/eqv/Est_result_Singletons.RData")

suppressMessages(suppressWarnings(library(Biostrings)))
pan_genome = readDNAStringSet(panfn)
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


#nghiavtr-03/12/2025 - begin: fix the mismatched names
pan_name2 = sapply(pan_name, function(s) {
  s1 = unlist(strsplit(s, "\\."))
  s1[1]=gsub("PAN","",s1[1])
  paste0(s1[1], s1[3])
})
names(pan_name2)=NULL
pan_long2short=pan_shortname
names(pan_long2short)=pan_name2
pan_short2long=pan_name2
names(pan_short2long)=pan_shortname
#nghiavtr-03/12/2025 - end


convertPanName <-function(mat){
  xname=rownames(mat)
  xname1=xname
  
  y1=strsplit(xname," ")
  l1=lengths(y1)
  p1=which(l1==1)
  p2=which(l1!=1)

  if (length(p1)>0){
    xname1[p1]=pan_long2short[xname[p1]]    
  }

  if (length(p2)>0){
    y2=sapply(y1[p2],function(x){
      x2=pan_long2short[x]
      x3=paste(x2,collapse=" ")
    })
    names(y2)=NULL
    xname1[p2]=y2
  }

  return(xname1)
}

rownames(result_est)=convertPanName(result_est)

final_singleton = rbind(result_est, isolate_singleton)
save(final_singleton, file = "output/final_singleton.RData")
Y = final_Y
save(Y, file = "output/final_Y.RData")
