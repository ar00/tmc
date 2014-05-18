
# date: 16/05/2014
# aim: randomly sampled diagnostic plots (res vs fitted and QQnorm)
# from low to high M0/M1 RMSPE ratios. 
# - red corresponds to ratio in [1,2]
# - green corresponds to ratio in [3,5]
# - blue corresponds to ratio in [6,15]
# return:
# - ../results/model.diagnostic.pdf 


library(caret) # for resolution of linear dependencies in X1
library(GenomicRanges)
library(gplots)
library(plotrix)
library(impute)


load("../data/data.rda")
load("../data/summarizedExperiment.rda")
load("../data/CISBP-RNA.rda")
load("../data/annotation.rda")


source("fit_functions.R")


subset_and_fit_target <- function(target) {

 # RBPs inferred to bind either mRNA UTRs
 v.rbp=c()

 tmp=colnames(rbp.db)[rbp.db[match(
  annotation$HGNC.symbol[match(targetIDs[i],
  annotation$ID)],rownames(rbp.db)),]==1]

 v.rbp=intersect(annotation$ID[is.na(match(annotation$HGNC.symbol,
  tmp))==FALSE],rownames(rbp.prot))

 if(length(v.rbp)>0){ v=c("rna",v.rbp) }
 if(length(v.rbp)==0){ v="rna" }

 # skip genes for which no RBPs are inferred
 if(length(v.rbp)==0){ return(rep("NA",times=5)) }

 # check if the RNA levels are available
 if(all(is.na(t(rna.exprm[is.na(match(rownames(rna.exprm),
   targetIDs[i]))==FALSE,])))){ return(rep("NA",times=5)) }

 # define M1 explanatory matrix
 X1=matrix(0,nrow=dim(exprm)[2],ncol=length(v),
  dimnames=list(colnames(rbp.prot),v))
 # enter explanatory variable in common to M0 and M1 (gene's RNA level)
 X1[,1]=t(rna.exprm[is.na(match(rownames(rna.exprm),targetIDs[i]))==FALSE,])
 # enter additional explanatory variables in M1 (RBP protein levels)
 X1[,2:dim(X1)[2]]=t(rbp.prot[is.na(match(rownames(rbp.prot),v.rbp))==FALSE,])

 # define M0 explanatory vector
 X0=as.numeric(t(rna.exprm[match(targetIDs[i],
  rownames(rna.exprm)),]))

 # define response variable
 Y=as.numeric(t(prot.exprm[match(targetIDs[i],
  rownames(prot.exprm)),]))

 # imputation of N/A values
 Z=t(as.matrix(cbind(Y,X1)))
 Z.knn = t(impute.knn(Z)$data)
 Y=Z.knn[,1]
 X1=Z.knn[,-1]

 m1 = fit(Y,X1)
 m0 = fit(Y,X0)

 return(list(m0,m1,X0,X1,Y))

}


# obtain RBP->target interactions (binary value 1/0)
rbp.db=read.delim("../data/rbp.target.binary.txt",
 header=TRUE,as.is=TRUE,sep="\t")


# obtain NCI-60 protein data
exprm = as.matrix(assays(sset)[[2]])
ind = apply(exprm, 1, function(x) all(is.na(x)))
# remove N/A rows inserted in the construction of summarizedExperiment object
exprm = exprm[!ind,]
# replace 0 values with N/A according to the specification by Prof. Kuster
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(exprm[i,j]==0){ exprm[i,j]=NA }
 }
}
# remove N/A rows from NCI-60 protein data
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]
# rename protein data
prot.exprm=exprm
# obtain RBP protein data
rbp.prot=exprm[is.na(match(rownames(exprm),
 annotation$ID[is.na(match(annotation$HGNC.symbol,rbp.pwm))==FALSE]))==FALSE,]


# obtain RNA data
exprm = as.matrix(assays(sset)[[1]])
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]
rna.exprm=exprm


# select proteins corresponding to inferred RBP targets
tmp=annotation$ID[is.na(match(annotation$HGNC.symbol,rownames(rbp.db)))==FALSE]
targetIDs = rownames(prot.exprm)[is.na(match(rownames(prot.exprm),tmp))==FALSE]


res=data.frame() # relative prediction error
linearComb=list()
remove=list()
for(i in 1 : length(targetIDs)){

 models=subset_and_fit_target(targetIDs[i])
 if(is.list(models)==T){
  m0.rmse = sqrt(mean(residuals(models[[1]])^2)) 
  m1.rmse = sqrt(mean(residuals(models[[2]])^2))
  rmse.ratio = m0.rmse/m1.rmse

  # get RMSE ratios betwen M0 and M1
  res = rbind(res,cbind(targetIDs[i],rmse.ratio))

  # detect linear combinations among predictors
  linearComb[[targetIDs[i]]]=findLinearCombos(models[[4]])
  remove[[targetIDs[i]]]=colnames(models[[4]])[findLinearCombos(models[[4]])$remove]
 }
 show(i)
}

##############################################
# linear combinations in the predictor matrix
##############################################

pdf("../results/linear.comb.pdf")
par(mfrow=c(1,2))
# enumerate variables to remove according to the resolution of 
# linear combinations in the predictor matrix
tmp=as.numeric(lapply(remove,function(x)(length(x))))
hist(tmp,breaks=25,xlab="N. removals against linear combinations",prob=T)

# enumerate linear combinations in the predictor matrix
tmp=as.numeric(lapply(linearComb,function(x)(length(x$linearCombos))))
hist(tmp,breaks=25,xlab="N. linear combinations",prob=T)
dev.off()


###################################################
# diagnostic plots for randomly sampled gene models
###################################################

# sample gene models from low to high RMSE ratio
res.rmse.ratio.1 = res[as.numeric(as.vector(res[,2]))>=1 & 
 as.numeric(as.vector(res[,2]))<2,]
res.rmse.ratio.2 = res[as.numeric(as.vector(res[,2]))>=3 & 
 as.numeric(as.vector(res[,2]))<5,]
res.rmse.ratio.3 = res[as.numeric(as.vector(res[,2]))>=6 &
 as.numeric(as.vector(res[,2]))<15,]


pdf("../results/model.diagnostic.pdf")
par(mfrow=c(3,2))

# first block of samples
s=sample(size=27,x=1:dim(res.rmse.ratio.1)[1],replace=FALSE)
for(i in 1 : length(s)){
 m = subset_and_fit_target(res.rmse.ratio.1[i,1])[[2]]
 if(!all(m$residuals==0)){
  title=paste(res.rmse.ratio.1[i,1],
  annotation$HGNC.symbol[match(res.rmse.ratio.1[i,1],annotation$ID)],sep="-")
  plot(predict(m),m$residuals,main=title,cex=.75,pch=19,col=rgb(1,0,0),
   cex.lab=1.75,cex.axis=1.25)
  mtext("1<= RMSE ratio <2",side=3,col=rgb(1,0,0),line=-1.5,cex=.75)
  qqnorm(rstandard(m),main=title,cex=.75,pch=19,cex.lab=1.75,cex.axis=1.25,col=rgb(1,0,0))
  mtext("1<= RMSE ratio <2",side=3,col=rgb(1,0,0),line=-1.5,cex=.75)
 }
}

# second block of samples
s=sample(size=27,x=1:dim(res.rmse.ratio.2)[1],replace=FALSE)
for(i in 1 : length(s)){
 m = subset_and_fit_target(res.rmse.ratio.2[i,1])[[2]]
 if(!all(m$residuals==0)){
  title=paste(res.rmse.ratio.2[i,1],
  annotation$HGNC.symbol[match(res.rmse.ratio.2[i,1],annotation$ID)],sep="-")
  plot(predict(m),m$residuals,main=title,cex=.75,pch=19,col=rgb(0,1,0),
   cex.lab=1.75,cex.axis=1.25)
  mtext("3<= RMSE ratio <5",side=3,col=rgb(0,1,0),line=-1.5,cex=.75)
  qqnorm(rstandard(m),main=title,cex=.75,pch=19,cex.lab=1.75,cex.axis=1.25,col=rgb(0,1,0))
  mtext("3<= RMSE ratio <5",side=3,col=rgb(0,1,0),line=-1.5,cex=.75)
 }
}

# third block of samples
s=sample(size=27,x=1:dim(res.rmse.ratio.3)[1],replace=FALSE)
for(i in 1 : length(s)){
 m = subset_and_fit_target(res.rmse.ratio.3[i,1])[[2]]
 if(!all(m$residuals==0)){
  title=paste(res.rmse.ratio.3[i,1],
  annotation$HGNC.symbol[match(res.rmse.ratio.3[i,1],annotation$ID)],sep="-")
  plot(predict(m),m$residuals,main=title,cex=.75,pch=19,col=rgb(0,0,1),
   cex.lab=1.75,cex.axis=1.25)
  mtext("6<= RMSE ratio <15",side=3,col=rgb(0,0,1),line=-1.5,cex=.75)
  qqnorm(rstandard(m),main=title,cex=.75,pch=19,cex.lab=1.75,cex.axis=1.25,col=rgb(0,0,1))
  mtext("6<= RMSE ratio <15",side=3,col=rgb(0,0,1),line=-1.5,cex=.75)
 }
}

dev.off()



