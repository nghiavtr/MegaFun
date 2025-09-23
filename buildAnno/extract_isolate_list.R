rm(list = ls())
# Rscript extract_isolate_list.R speciesList=/nfs/NGS/MegaFun_proj/Nghia_test/MegaFun_anno_prep/metaphlan_out/metaphlan_identified_species.txt out=./ ncbiSummary=ncbi_assembly_summary.txt.gz ncbiName=ncbi_name_db.RData  chocoMap=choco_list_full_chocophlan.v201901_v31.RData

#get the list of isolates for download
#output: isolate_genome_summary.txt


ncbi_db_fn="ncbi_assembly_summary.txt.gz" #download from ncbi
species_list_fn="metaphlan_out/metaphlan_identified_species.txt" #output of MegaFun_anno_pipeline.sh

chocoMapFn="choco_list_full_chocophlan.v201901_v31.RData" #output of extract_ncbi_ID_chocophlan_full_chocophlan.v201901_v31.R

ncbi_name_db_fn="ncbi_name_db.RData"


outDir="./"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="speciesList") species_list_fn=as.character(res[2])
    if (res[1]=="out") outDir=as.character(res[2])
    if (res[1]=="ncbiSummary") ncbi_db_fn=as.character(res[2])
    if (res[1]=="ncbiName") ncbi_name_db_fn=as.character(res[2])
    if (res[1]=="chocoMap") chocoMapFn=as.character(res[2]) #chocophlan_human2
    
  }



library(BSgenome)
library(parallel)
library(Rsamtools)
library(data.table)
#install.packages("stringr")
library(stringr)
library(Biostrings) #stringDist

if (!dir.exists(outDir)) dir.create(outDir)

load(ncbi_name_db_fn)
dim(ncbi_name_db)

ncbi_db=fread(ncbi_db_fn, header = TRUE, quote="" )
x=ncbi_db$organism_name
x1=gsub(" ","_",x)
ncbi_db$species=x1
ncbi_db=as.data.frame(ncbi_db)

#get mapping info of chocophlan - prepared in advance using extract_ncbi_ID_chocophlan.R
load(chocoMapFn)
dim(choco_list)



species_list=fread(species_list_fn, header=FALSE)
if (ncol(species_list)==1){ #output of MetaPhlan
  species_list=as.data.frame(species_list)
  colnames(species_list)="id"
}

if (ncol(species_list)==2){ #manually list with species name and taxID, no header is allowed
  species_list=as.data.frame(species_list)
  colnames(species_list)=c("name","ncbi_ID") #species_name tax_id
}
#head(species_list)

#remove duplicated - nghiavtr/23June2025
p=which(duplicated(species_list[,1]))
if (length(p)>0) species_list=species_list[-p,drop=FALSE,]

sp_list=as.data.frame(species_list)

if(ncol(sp_list)==1){
  sp_list$name <- sapply(sp_list[,1], function(s) {
    temp <- unlist(strsplit(s, "__"))[3]
    gsub("_", " ", temp)
  })

  #check matching between species list (output of metaphlan - very messed up) and chocophlan db
  #first, check the name from ncbi
  m=match(sp_list$name,ncbi_name_db$Name)
  sp_list$ncbi_ID=ncbi_name_db$ID[m]
}



#check if there are unmatched names, if so, check more from chocophlan
p=which(is.na(sp_list$ncbi_ID))
#p
if (length(p)>0){
  m1=match(sp_list$name,choco_list$species_name)
  sp_list$ncbi_ID[p]=choco_list$ID[m1][p]
  # p=which(is.na(sp_list$ncbi_ID))
  # if (length(p)>0){
  #   m2=match(sp_list$name,choco_list_v201901$species_name)
  #   sp_list$ncbi_ID[p]=choco_list_v201901$ID[m2][p]
  # }
}

#sum(is.na(sp_list$ncbi_ID)) #no unmatched data anymore

chocoID=sapply(choco_list[,1],function(x){
  x1=strsplit(x,"\\.")[[1]]
  paste0(x1[1],".",x1[2])
  })
names(chocoID)=NULL

choco_list$chocoID=chocoID

#check if sp_list are overlapping with chocophlan - keep them
p=sp_list$ncbi_ID %in% choco_list$ID | sp_list$name %in% choco_list$chocoID
sp_list=sp_list[p,]
sp_list$ncbi_name=ncbi_name_db$Name[match(sp_list$ncbi_ID, ncbi_name_db$ID)]
sp_list$choco_name=choco_list$species[match(sp_list$ncbi_ID, choco_list$ID )]
sp_list$choco_filename=choco_list[match(sp_list$ncbi_ID, choco_list$ID ),1]

p1=which(is.na(sp_list$choco_name))
if (length(p1) > 0){
  sp_list$choco_name=choco_list$species[match(sp_list$name, choco_list$chocoID )]
  sp_list$choco_filename=choco_list[match(sp_list$name, choco_list$chocoID ),1]
}

chocoID=sapply(sp_list$choco_filename,function(x){
  x1=strsplit(x,"\\.")[[1]]
  paste0(x1[1],".",x1[2])
  })
names(chocoID)=NULL
sp_list$id=chocoID
#head(sp_list)

#remove if ncbi_ID is NA - nghiavtr/23June2025
sp_list=sp_list[!is.na(sp_list$ncbi_ID),]

### now collect isolates from ncbi_db based on ncbi_ID

# get all isolate genomes for 20 species
#p=ncbi_db$taxid %in% sp_list$ncbi_ID
p=ncbi_db$species_taxid %in% sp_list$ncbi_ID #species_taxid can be different from taxid
isolate_df <- ncbi_db[p, ]


isolate_df$infraspecific_name <- sapply(isolate_df$infraspecific_name, function(s) {
  strain <- gsub("strain=", "", s)
  strain <- gsub("-", "", strain)
  strain <- gsub(" ", "", strain)
  strain <- gsub("_", "", strain)
  strain <- gsub("/", "", strain)
  strain <- gsub("\\", "", strain, fixed = TRUE)
})

#pick = match(isolate_df$organism_name, sp_list$name)
#pick = match(isolate_df$taxid, sp_list$ncbi_ID)
pick = match(isolate_df$species_taxid, sp_list$ncbi_ID) #species_taxid can be different from taxid
isolate_df$species_id = sp_list$id[pick]

#isolate_df$species_choco_id=sp_list$choco_name[match(isolate_df$taxid,sp_list$ncbi_ID)]
#isolate_df$species_choco_filename=sp_list$choco_filename[match(isolate_df$taxid,sp_list$ncbi_ID)]
#species_taxid can be different from taxid
isolate_df$species_choco_id=sp_list$choco_name[match(isolate_df$species_taxid,sp_list$ncbi_ID)]
isolate_df$species_choco_filename=sp_list$choco_filename[match(isolate_df$species_taxid,sp_list$ncbi_ID)]


#save data that need to download to download
write.table(isolate_df, file = paste0(outDir,"/isolate_genome_summary.txt"), quote = F, sep = "\t")

sp_list$ncbi_sname=gsub(" ","_",sp_list$ncbi_name)
sp_list$ncbi_sname=paste0("s__",sp_list$ncbi_sname)
#do.call(rbind,strsplit(sp_list$id,"\\."))[,2]
sp_list$choco_sname=paste0("s__",sp_list$choco_name)
write.table(sp_list, file = paste0(outDir,"/species_list_summary.txt"), quote = F, sep = "\t")


