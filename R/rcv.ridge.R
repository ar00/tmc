
# date: 16/05/2014 
# aim: 
#  To compare extra-samples prediction error achieved by M0, M1 by 
#  repeated (R=10) K-fold (K=5) cross-validation (RCV).
# 
#  Note that we add elastic-net regularization path for linear regression
#  (we used cv.glmnet() within the glmnet package with nfolds=3 and alpha=0).
#
#  Missing data are imputed, separately in training and test subsets
#  The function for RCV has the option of trimming predctor matrix if
#  linear dependencies among predictors are recorded.
# return: 
# - ../results/rcv.summary.ridge.penalty.txt


library(GenomicRanges)
library(gplots)
library(cvTools)
library(plotrix)
library(impute)


load("../../results/data.rda")
load("../../results/summarizedExperiment.rda")
load("../../annotation/CISBP-RNA.rda")
load("../../results/annotation.rda")


source("fit_functions.R")


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
  if(exprm[i,j]==0){ exprm[i,j] = NA }
 }
}
# remove full N/A rows from NCI-60 protein data
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

# summary of prediction error estimates achieved by M0 and M1
RCV <- data.frame(t(rep(NA,14)))
RCVcolnames <- c("ID","HGNC.symbol","m0.rmspe","m1.rmspe",
 "m0.coef.n","m1.coef.n","cvm","cvsd","lambda.min",
 "rmspe.ratio","m0.rmspe.avg","m1.rmspe.avg",
 "rmspe.ratio.avg","n.iter")
RCV <- RCV[-1,]

# select proteins corresponding to inferred RBP targets
tmp=annotation$ID[is.na(match(annotation$HGNC.symbol,rownames(rbp.db)))==FALSE]
targetIDs = rownames(prot.exprm)[is.na(match(rownames(prot.exprm),tmp))==FALSE]

for(i in 1 : length(targetIDs)){

 show(i)
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
 if(length(v.rbp)==0){ next }

 # check if the RNA levels are available
 if(all(is.na(t(rna.exprm[is.na(match(rownames(rna.exprm),
   targetIDs[i]))==FALSE,])))){ next }

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

 #############################
 # cross-validation
 #############################

 K=5; R=10
 
 cv <- cross_validation_with_penalty_tunel(K=K,R=R,X0=X0,X1=X1,Y=Y)
 if(is.list(cv)==TRUE) {
  RCV=rbind(RCV,cbind(as.character(targetIDs[i]),
   as.character(annotation$HGNC.symbol[match(targetIDs[i],annotation$ID)]),
   cv$m0.rmspe,cv$m1.rmspe,
   cv$m0.coef.n,cv$m1.coef.n,cv$cvm,cv$cvsd,cv$lambda.min, 
   cv$rmspe.ratio,
   cv$m0.rmspe.avg,cv$m1.rmspe.avg,
   cv$rmspe.ratio.avg,
   cv$n.iter))
 }
}
write.table(RCV,file="../results/rcv.summary.ridge.penalty.txt",
 quote=F,row.names=F,col.names=RCVcolnames,sep="\t")
    


