### construct mapping files
rm(list = ls())
### 1) download database from here: 
# wget https://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
# tar -xzf full_chocophlan.v201901_v31.tar.gz
### 2) create the list of files in the database
# ls full_chocophlan.v201901_v31 > full_chocophlan.v201901_v31.txt
### 3) download taxdump from ncbi
# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# tar -xzf taxdump.tar.gz

### run
# Rscript contructTaxChocoDB.R chocoList=full_chocophlan.v201901_v31.txt ncbiNames=taxdump/citations.dmp ncbiCites=taxdump/citations.dmp mapOut=choco_list_full_chocophlan.v201901_v31.RData ncbidbOut=ncbi_name_db.RData

choco_list_fn="full_chocophlan.v201901_v31.txt"
ncbi_cites_fn="taxdump/citations.dmp"
ncbi_names_fn="taxdump/names.dmp"
map_fn="choco_list_full_chocophlan.v201901_v31.RData"
ncbi_db_fn="ncbi_name_db.RData"


args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

if (length(args)>0)
  for (i in 1:length(args)){
    res=unlist(strsplit(args[i],"="))
    if (res[1]=="chocoList") choco_list_fn=as.character(res[2])
    if (res[1]=="ncbiNames") ncbi_names_fn=as.integer(res[2])
    if (res[1]=="ncbiCites") ncbi_cites_fn=as.character(res[2])
    if (res[1]=="mapOut") map_fn=as.character(res[2])
    if (res[1]=="ncbidbOut") ncbi_db_fn=as.character(res[2])
    
  }


library(data.table)
#install.packages("stringr")
library(stringr)
library(Biostrings) #stringDist

choco_list=fread(choco_list_fn, header = FALSE)
choco_list=as.data.frame(choco_list)
s=strsplit(choco_list[,1],"\\.s__")
s=do.call(rbind,s)[,2]
#s=gsub(".centroids.ffn.m8.gz","",s)
s=strsplit(s,"\\.centroids")
s=do.call(rbind,s)[,1]
choco_list$species=s
choco_list$species_name=gsub("_"," ",choco_list$species)

x=choco_list$species_name
x=str_replace_all(x, "[^[:alnum:]]", "") #very convenient from stringr
choco_list$species_name_clean=x

#get data from 
ncbi_cites=fread(ncbi_cites_fn)
ncbi_names=fread(ncbi_names_fn)

x=ncbi_names$V2
x=gsub("\t","",x)
x=gsub('"',"",x)
ncbi_names$V2=x
x=ncbi_names$V4
x=gsub("\t","",x)
x=gsub('"',"",x)
ncbi_names$V4=x
ncbi_names=ncbi_names[,c(1,2,4)]
head(ncbi_names)
colnames(ncbi_names)=c("ID","Name","Type")

x=ncbi_names$Name
x=str_replace_all(x, "[^[:alnum:]]", "")
ncbi_names$species_name_clean=x
head(ncbi_names)
x=NULL
##### now, do matching

head(choco_list)
o=choco_list$species_name_clean %in% ncbi_names$species_name_clean
table(o) #170 unmatched
head(choco_list[!o,])
m=match(choco_list$species_name_clean,ncbi_names$species_name_clean)
p=which(!is.na(m))
#create ncbi ID
choco_list$ID=NA
choco_list$ID[p]=ncbi_names$ID[m[p]]
head(choco_list)

### consider 210 unmatched cases
choco_list_um=choco_list[!o,]
choco_list_um$species_name_clean2=choco_list_um$species_name_clean

#Case1: Bacteroides = Phocaeicola
#pp=grep("Xanthomonadaceae",choco_list_um$species_name_clean)
choco_list_um$species_name_clean2=gsub("Xanthomonadaceae","Lysobacteraceae",choco_list_um$species_name_clean2)
pp=which(choco_list_um$species_name_clean=="enodsymbiontofEscarpiaspicata")
choco_list_um$species_name_clean2[pp]="endosymbiontofEscarpiaspicata"
pp=which(choco_list_um$species_name_clean=="Lactobacilusnuruki")
choco_list_um$species_name_clean2[pp]="Lactobacillusnuruki"
choco_list_um$species_name_clean2=gsub("Bacteroides","Phocaeicola",choco_list_um$species_name_clean2)


o1=choco_list_um$species_name_clean2 %in% ncbi_names$species_name_clean
table(o1)

table(o1) #11 cases resolved
m=match(choco_list_um$species_name_clean2,ncbi_names$species_name_clean)
p=which(!is.na(m))
choco_list_um$ID[p]=ncbi_names$ID[m[p]]
choco_list_um1=choco_list_um[!o1,] #205  remaining
dim(choco_list_um1)


#Case 3: check from citation file, if not use the best match for remaining cases, for example there are authorname et al. in the names


### slow, but ok with 162 cases
i=50
choco_list_um1$ID=NA
for ( i in 1:nrow(choco_list_um1)){
  cat(" ",i)
#  choco_list_um1[i,]
  x=NULL
  p1=grep(choco_list_um1$species_name[i], ncbi_cites$V11)
  cite_ID=ncbi_cites$V13[p1]
#  cite_ID
  if (length(cite_ID) > 0){
    cite_ID=unique(strsplit(cite_ID," ")[[1]])
    cite_ID=as.integer(cite_ID)
    x=ncbi_names[which(ncbi_names$ID %in% cite_ID),]
  }
#  x
  if (length(cite_ID) == 0 || nrow(x)==0){
    nameComp=strsplit(choco_list_um1$species_name[i]," ")[[1]]
    nameComp=ifelse(nameComp=="sp","sp.",nameComp) #use sp. to get exact match
    ## search with the longest components as possible
    p=NULL
    if (length(nameComp)>1){
      for (j in (length(nameComp)):1){
  #      cat(" ",j)
        subName=paste(nameComp[1:j], collapse ="")
        p=grep(subName,ncbi_names$species_name)
        if (length(p)>0){
          break()
        }
      }
    }
    
    x=ncbi_names[p,]
  }


  keep=which(x$Type %in% c("scientific name"))
  if (length(keep)>0)  x=x[keep,] else{
    keep=which(x$Type %in% c("synonym"))
    if (length(keep)>0)  x=x[keep,] else{
      keep=which(x$Type %in% c("equivalent name"))
      if (length(keep)>0)  x=x[keep,]
    }
  }
  
  
  #compute distance between strings, put the choco name at the bottom
  d=pairwiseAlignment(x$Name, choco_list_um1$species_name[i])#from Biostrings
  k=which.max(score(d))
  if (length(k)>0)  choco_list_um1$ID[i]=x$ID[k] #keep the ID
}

#update
choco_list_um1$ID
choco_list_um$ID[!o1]=choco_list_um1$ID
choco_list$ID[!o]=choco_list_um$ID
sum(is.na(choco_list$ID)) #only 1 left

head(ncbi_names)
head(choco_list)
ncbi_names_standard=ncbi_names[ncbi_names$Type=="scientific name",]
dim(ncbi_names_standard)
m=match(choco_list$ID, ncbi_names_standard$ID)
choco_list$ncbi_scientifcName=ncbi_names_standard$Name[m]
#save(choco_list,file="choco_list_full_chocophlan.v201901_v31.RData")
save(choco_list,file=map_fn)


ncbi_name_db=as.data.frame(ncbi_names_standard)
ncbi_name_db=ncbi_name_db[,-c(3,4)]
rownames(ncbi_name_db)=NULL
#save(ncbi_name_db, file="ncbi_name_db.RData")
save(ncbi_name_db, file=ncbi_db_fn)


