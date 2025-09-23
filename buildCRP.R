# 26 Feb 2023 / Nghia: 
# - break down huge matrices to smaller ones by gradually increasing H_thres which started from 0.005
## 15 Otc 2019 / Nghia:
# - fix a bug when separe a CRP into more than 1 CRP
## 30 Apr 2019 / Nghia:
# - compute CRP from CRPCOUNT: fix the bug when colSums(mycrpCount)==0 
## 03 Apr 2019 / Nghia:
# - add H_thres for filtering connections with low counts between two transcripts
# Example of command: Rscript in=eqClass.txt out=X_matrix.RData H=0.025 s=500
## 02 Nov 2018 / Nghia:
# - grow fully TC, 
# - improve speed

### Build transcript cluster and CRP from simulated data
# Before running this script: 
# 1) generate simulated data of the whole transcriptome using genPolyesterSimulation.R
# 2) run GenTC for the simulated data
# Input: eqClass.txt from results of GenTC


### default settings
#eqClassFn="eqClass.txt"
rm(list = ls())
# eqClassFn="X_matrix/pseudo_alignment/eqClass.txt"
# fout='X_matrix/pseudo_alignment/X_matrix.RData'
eqClassFn="X_matrix/eqClass.txt"
fout='X_matrix/X_matrix.RData'
#H_thres=0.0
H_thres=0.025
H_thres_min=0.01
targetsize=500


source("/path/to/MegaFun/Rsource.R")

#core = 8 #default
library(data.table)
library(foreach)
library(doParallel)
ncores = detectCores()
#nc = min(ncores,core)     # use 8 or 16 as needed!!
nc=16
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)


#get input from eqClass.txt in the output directory of GenTC
#lines = readLines(eqClassFn)
#lines = gsub("\'", "", lines)
#rawmat  = read.table(textConnection(lines), header=TRUE, as.is=TRUE, sep='\t')
rawmat  = fread(eqClassFn,sep='\t',header=TRUE)
rawmat  = as.data.frame(rawmat)
# dim(rawmat)

# rawmat = rawmat[rawmat$Count >= 15, ]
# dim(rawmat)


# nghia/24Feb2023: remove noisy eqclasses that read counts contribute too litle to all transcripts (< 0.001)
txTrueCount=tapply(rawmat$Weight,rawmat$Transcript, sum)
rawmat$txTrueCount=txTrueCount[rawmat$Transcript]
rawmat$prop=rawmat$Weight/rawmat$txTrueCount
rawmat$propEq=rawmat$Count/rawmat$txTrueCount
minProp=tapply(rawmat$prop,rawmat$eqClass,max)
minPropEq=tapply(rawmat$propEq,rawmat$eqClass,max)
pick=which(minProp > H_thres_min | minPropEq > H_thres_min) #excluded if reads of individual transcripts in the eqclass contributes < H_thres
rawmat=rawmat[rawmat$eqClass %in% as.integer(names(pick)),]
dim(rawmat)

#reindex
myeqID=unique(rawmat$eqClass)
myeqID=sort(myeqID)
matchID=match(rawmat$eqClass,myeqID)
rawmat$eqClass=matchID
######


tx2eqc = tapply(rawmat$eqClass,rawmat$Transcript,c)  ## map: tx --> eqClass 
a = sapply(tx2eqc, length)
#table(a)
eqc2tx = tapply(rawmat$Transcript, rawmat$eqClass,c) ## map: eqClass --> tx

# neighbors of each tx
fn = function(i) {eqc=as.character(tx2eqc[[i]]);
unique(unlist(eqc2tx[eqc]))
}


#system.time(NB <- foreach(i=1:length(tx2eqc)) %dopar% fn(i) )  ## about 25sec 
#names(NB) = names(tx2eqc)

NB=lapply(tx2eqc,function(x) unique(unlist(eqc2tx[x])))

### build transcript clusters - new codes

#convert name of tx to index
txi=seq(length(NB))
NBi=NB
names(NBi)=txi
#map from tx to crp
t2c_map=unlist(NBi,use.names=FALSE)
t2c_mapi=match(t2c_map,names(NB))
t2c_mapi_group=rep(txi, lengths(NBi))
NBi2=tapply(t2c_mapi,t2c_mapi_group,c)
NBi=NBi2

f1=sapply(NBi,function(x) min(x))
f2=sapply(NBi,function(x) min(f1[x]))

growTimes=1
isOK=FALSE
repeat{  
  f2=sapply(NBi,function(x) min(f1[x]))  
  difNum=sum(f2!=f1)  
  growTimes=growTimes+1
  f1=f2
  if (difNum==0) isOK=TRUE
  if (isOK) break();
}

f1=names(NB)[f1]
names(f1)=names(NB)

NB2=tapply(names(f1),f1,c)
NB2=lapply(NB2,function(x) sort(x))
#create TC3
TC3=NB
TC3[names(f1)]=NB2[f1]

############ 
# nghia/26Feb2023: break down huge matrices to smaller ones by gradually increasing H_thres
 getTC3<-function(matIn){
     #reindex
     myeqID=unique(matIn$eqClass)
     myeqID=sort(myeqID)
     matchID=match(matIn$eqClass,myeqID)
     matIn$eqClass=matchID
     ######
     tx2eqc = tapply(matIn$eqClass,matIn$Transcript,c)  ## map: tx --> eqClass 
     a = sapply(tx2eqc, length)
     #table(a)
     eqc2tx = tapply(matIn$Transcript, matIn$eqClass,c) ## map: eqClass --> tx

     NB=lapply(tx2eqc,function(x) unique(unlist(eqc2tx[x])))

     #convert name of tx to index
     txi=seq(length(NB))
     NBi=NB
     names(NBi)=txi
     #map from tx to crp
     t2c_map=unlist(NBi,use.names=FALSE)
     t2c_mapi=match(t2c_map,names(NB))
     t2c_mapi_group=rep(txi, lengths(NBi))
     NBi2=tapply(t2c_mapi,t2c_mapi_group,c)
     NBi=NBi2

     f1=sapply(NBi,function(x) min(x))
     f2=sapply(NBi,function(x) min(f1[x]))

     growTimes=1
     isOK=FALSE
     repeat{  
       f2=sapply(NBi,function(x) min(f1[x]))  
       difNum=sum(f2!=f1)  
       growTimes=growTimes+1
       f1=f2
       if (difNum==0) isOK=TRUE
       if (isOK) break();
     }

     f1=names(NB)[f1]
     names(f1)=names(NB)

     NB2=tapply(names(f1),f1,c)
     NB2=lapply(NB2,function(x) sort(x))
     #create TC3
     TC3=NB
     TC3[names(f1)]=NB2[f1]

     return(list(TC3=TC3,NB=NB,eqc2tx=eqc2tx,tx2eqc=tx2eqc,rawmat=matIn))
   } #end of the function

 len_TC3=lengths(TC3)
 bigTC3=which(len_TC3>targetsize)

 if (length(bigTC3) > 0){
     rmEqSet=NULL
     H_thres1=H_thres_min
     rawmat1=rawmat
     rawmat1$rawEq=rawmat1$eqClass
     repeat{
       H_thres1=H_thres1+0.005 #increase 0.005 each time
       #if (length(bigTC3) == 0 | H_thres1 > H_thres) break()
       if (length(bigTC3) == 0 ) break() # hard filter by matrix size, dont care about H_thres to make sure it will not be out of memory

       cat("\n",H_thres1)
       mytx=names(bigTC3)
       myeq=rawmat1$rawEq[rawmat1$Transcript %in% mytx]
       myeq=unique(myeq)

       rawmat2=rawmat1[rawmat1$rawEq %in% myeq,]
       minProp=tapply(rawmat2$prop,rawmat2$rawEq,max)
       minPropEq=tapply(rawmat2$propEq,rawmat2$rawEq,max)
       pick=which(minProp <= H_thres1 & minPropEq <= H_thres1) #excluded if reads of individual transcripts in the eqclass contributes < H_thres

       rmEqSet=c(rmEqSet,as.integer(names(pick)))

       pick=which(minProp > H_thres1 | minPropEq > H_thres1) #excluded if reads of individual transcripts in the eqclass contributes < H_thres
       rawmat3=rawmat2[rawmat2$rawEq %in% as.integer(names(pick)),]

       res=getTC3(rawmat3)
       rawmat1=res$rawmat   
       len_TC3=lengths(res$TC3)
       bigTC3=which(len_TC3>targetsize)
       cat("\n",table(len_TC3))
     }#end of repeat

   #update new raw data
   if (length(rmEqSet)>0){
     pick=which(!rawmat$eqClass %in% rmEqSet)
 #    ### process the rmEqSet
 #    for (eq in rmEqSet){
 #      x=rawmat_bk[rawmat_bk$eqClass==eq,]
 #    }
     #update rawmat
     rawmat=rawmat[pick,]
   }
   res=getTC3(rawmat)  
   NB=res$NB
   TC3=res$TC3
   eqc2tx=res$eqc2tx
   tx2eqc=res$tx2eqc
   rawmat=res$rawmat


   }#end of if

############


OTC <- sapply(TC3, paste, collapse=' ') # pasted version
names(OTC) = names(NB)
otcmap = as.list(OTC)
names(otcmap) = names(NB)

#output: list of clusters and tx->cluster map 
clust = names(table(OTC))  # clusters
OTC = strsplit(clust,split=' ')
names(OTC) = clust


#get CRP count
#system.time(
CRPCOUNT <- foreach(i=1:length(OTC)) %dopar%{
#CRPCOUNT=list()
#for (i in 1:length(OTC)){
  if (i %% 100 == 0) cat(" ",i)
  myOTC=unlist(OTC[i])
  names(myOTC)=NULL
  #get binary codes of eqc
  myeqc=unique(unlist(sapply(myOTC,function(x) tx2eqc[x])))
  #eqc2tx[myeqc]
  bcode=lapply(eqc2tx[myeqc],function(x) as.integer(!is.na(match(myOTC,x))))
  #bcode

  #get corresponding count - new codes
  pick=which(rawmat$eqClass %in% myeqc)
  countname=paste(rawmat$Transcript[pick],rawmat$eqClass[pick],sep="__")
  count=rawmat$Weight[pick] 
  myname=expand.grid(myOTC,myeqc)
  myname=paste(myname[,1],myname[,2],sep="__")
  myCount=rep(0,length(myname))
  matchID=match(countname,myname)
  myCount[matchID]=count

  #create TC count matrix
  mycrpCount=matrix(unlist(myCount),nrow=length(bcode),ncol=length(myOTC),byrow=TRUE)
  colnames(mycrpCount)=myOTC
  rownames(mycrpCount)=sapply(bcode, function(x) paste(x,collapse=""))
  
  if (nrow(mycrpCount)>1){ # if there are more than 1 row
    mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
    #check if duplicated rownames
    repID=table(rownames(mycrpCount))
    repID=repID[which(repID>1)]
    if (length(repID)>0){
      rmID=which(rownames(mycrpCount) %in% names(repID))
      sumcrpCount=NULL
      for (j in 1:length(repID)){
        sumcrpCount=rbind(sumcrpCount,colSums(mycrpCount[rownames(mycrpCount) %in% names(repID)[j],]))
      }
      rownames(sumcrpCount)=names(repID)
      mycrpCount=mycrpCount[-rmID,]
      mycrpCount=rbind(mycrpCount,sumcrpCount)
      mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
      mycrpCount=as.matrix(mycrpCount)
    }
  }
#  mycrp=t(t(mycrpCount)/colSums(mycrpCount))
 mycrpCount
# CRPCOUNT[[i]]=mycrpCount
  
#  #extract CRP matrix
}

) #system

#### now compute CRP, use H_thres (default=0) to filter out too low proportion sharing between two transcripts

CRP=list()
for (k in 1:length(CRPCOUNT)){
  mycrpCount=CRPCOUNT[[k]]

  #fix the bug when colSums(mycrpCount)==0 
  txSum=colSums(mycrpCount)
  mycrpCount=mycrpCount[,which(txSum>0),drop=FALSE]
  if(ncol(mycrpCount)==0) next(); 
  
  #x=mycrpCount
  y=t(t(mycrpCount)/colSums(mycrpCount))
  z=apply(y,1,max)
  pick=z>H_thres
  if(sum(pick) > 0){
    x1=mycrpCount[pick,drop=FALSE,]
  } else {
    next
  }
  #decode eq in x1
  myclust=c(1:ncol(x1))
  for (i in 1:nrow(x1)){
    x2=which(x1[i,] >0)
    x3=which(myclust %in% myclust[x2])
    myclust[x3]=min(myclust[x3])
  }
  #generate new CRP
  clustID=as.integer(names(table(myclust)))
  #newx=list()
  for (i in 1:length(clustID)){  
    pick=which(myclust == clustID[i])
    X=mycrpCount[,pick,drop=FALSE]
    X=X[rowSums(X)>0,drop=FALSE,]
    rownames(X)=NULL
    bcode=apply(X,1,function(x) paste(as.integer(x>0),collapse=""))
    ubcode=unique(bcode)
    cbcode=match(bcode,ubcode) #bcode clustering   
    #redistribute values in row
    for (cID in unique(cbcode)){
      pick=which(cbcode==cID)
      if (length(pick)>1){
        for (j in 1:ncol(X)){
          X[pick,j]=sum(X[pick,j])
        }
      }
    }
    pick=which(!duplicated(cbcode))
    X=X[pick,drop=FALSE,]
    rownames(X)=bcode[pick]
    X_names=paste(colnames(X),collapse=" ")
    #normalise to get crp
    X=t(t(X)/colSums(X))    
    #add up results
    newx=list()
    newx[[X_names]]=X
    CRP=c(CRP,newx)
  }
}

CCRP <- lapply(CRP, ccrpfun)

#get txlength
txlength=as.integer(rawmat$RefLength)
names(txlength)=as.character(rawmat$Transcript)
pick=!duplicated(names(txlength))
txlength=txlength[pick]

#export CRP to file
save(CRP,CCRP,txlength,file=fout)



