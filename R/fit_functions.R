

library(caret)
library(cvTools)
library(glmnet)


  fit <- function (R,X){
   if(is.matrix(X)==FALSE){
    df=as.data.frame(cbind(R,X))
    colnames(df)=c("Y","Predictor")
    m=lm(Y  ~  Predictor, data=df)
   }
   if(is.matrix(X)==TRUE){
    df=as.data.frame(cbind(R,X))
    colnames(df)=c("Y",colnames(X))
    m=lm(Y  ~ ., data=df)
   }
   return(m)
  }

  fit_ridge <- function (R,X){
   # fit a linear model via penalized maximum likelihood (ridge penalty)
   glmnet.obj=glmnet(X,R,alpha=0)
   return(glmnet.obj)
  }

  fit_predict <- function (Rtrain,Rpredict,Xtrain,Xpredict){
   if(is.matrix(Xtrain)==FALSE){
    df=as.data.frame(cbind(Rtrain,Xtrain))
    colnames(df)=c("Y","Predictor")
    m=lm(Y  ~  Predictor, data=df)
    new.df <- data.frame(Predictor=Xpredict)
    Y.predict=predict.lm(m,new.df, interval="confidence",se.fit=T)
   }
   if(is.matrix(Xtrain)==TRUE){
    df=as.data.frame(cbind(Rtrain,Xtrain))
    colnames(df)=c("Y",colnames(Xtrain))
    m=lm(Y  ~ ., data=df)
    new.df <- data.frame(Xpredict)
    Y.predict=predict.lm(m,new.df, interval="confidence",se.fit=T)
   }
   return(Y.predict)
 }

 fit_ridge_predict <-function(Rtrain,Xtrain,Xpredict){
  # this functions predicts fitted values from a fitted "glmnet" object
  glmnet.obj=glmnet(Xtrain, Rtrain, alpha=0)
  new.df <- as.matrix(Xpredict)
  Y.predict=predict(glmnet.obj,new.df, interval="confidence",se.fit=T)
  return(Y.predict)
 }


# return estimate and confidence intervals of individual model ptredictors
conf_interval <- function (m,alpha){
 estimate = summary(m)$coefficients[,1]
 estimate.se = summary(m)$coefficients[,2]
 critical.value = qt(1-(alpha/2),summary(m1)$fstatistic[3])
 conf.interval = cbind(estimate,estimate-(critical.value*estimate.se),
  estimate+(critical.value*estimate.se))
 return(conf.interval)
}


# summarize estimates of prediction error by repeated K-fold cross-validation
cross_validation <- function (K,R,X0,X1,Y,trim){

 # K is folds number
 # R is repetitions number
 # X0 is the vector of observation of M0 predictor
 # X1 is the matrix of observations of M1 predictors
 # Y is the vector of observations of the gene to model
 # trim takes logical values about the choice of trimming predictors
 #  upon detection of possible linear combinations in predictors matrix

 rmspe.ratio=c();rmspe.ratio.se=c();
 m0.rmspe=c(); m0.rmspe.se=c(); m1.rmspe=c(); m1.rmspe.se=c();
 m0.cor=c();m1.cor=c();cor.ratio=c();cor.df=c();

 for(r in 1 : R){
  folds=cvFolds(length(Y), K = K, type = "random")
  tmp=cbind(folds$which,folds$subsets[,1])
  for(k in 1 : K){
   idx.predict=tmp[tmp[,1]==k,2]
   idx.train=tmp[tmp[,1]!=k,2]
   X0.train=X0[idx.train]
   if(is.matrix(X1)==T){X1.train=X1[idx.train,]}
   if(is.matrix(X1)==F){X1.train=X1[idx.train]}
   Y.train=Y[idx.train]
   X0.predict=X0[idx.predict]
   if(is.matrix(X1)==T){X1.predict=X1[idx.predict,]}
   if(is.matrix(X1)==F){X1.predict=X1[idx.predict]}
   Y.predict=Y[idx.predict]

   # imputation of N/A values in training sub-matrix
   Z1.train=t(as.matrix(cbind(Y.train,X1.train)))
   Z1.train.knn = t(impute.knn(Z1.train)$data)
   Y.train=Z1.train.knn[,1]
   X1.train=Z1.train.knn[,-1]

   # imputation of N/A values in prediction sub-matrix
   Z1.predict=t(as.matrix(cbind(Y.predict,X1.predict)))
   Z1.predict.knn = t(impute.knn(Z1.predict)$data)
   Y.predict=Z1.predict.knn[,1]
   X1.predict=Z1.predict.knn[,-1]

   if(trim==T){
    # remove predictor variables from training subsets
    if(!is.null(findLinearCombos(X1.train)$remove)){
     X1.train = X1.train[,setdiff(1:dim(X1.train)[2],
      findLinearCombos(X1.train)$remove)]
    }
    # remove predictor variables from testing subset
    if(!is.null(findLinearCombos(X1.train)$remove)){
     X1.predict = X1.predict[,setdiff(1:dim(X1.predict)[2],
      findLinearCombos(X1.predict)$remove)]
    }
  }

   # fit model M1
   m1=fit(Y.train,X1.train)
   
   m1.predict = fit_predict(Y.train,Y.predict,X1.train,X1.predict)
   if(!is.null(warnings())){ 
    write.table(paste(i,k,r,warnings(),sep=" "),
    file="rcv.knn.error.out",append=T,quote=F,row.names=F,col.names=F)
   }
   # M1 root mean squared prediction error (RMSPE)  
   rmspe.M1 = NA
   rmspe = signif(as.numeric(rmspe(Y.predict, m1.predict$fit[,1],
    includeSE=T)[1]),digits=3)
   if(is.na(rmspe)==FALSE){ rmspe.M1 = rmspe }
   m1.rmspe = c(m1.rmspe,rmspe.M1)

   # standard error of M1 RMSPE 
   rmspe.se.M1 = NA
   rmspe.se = signif(as.numeric(rmspe(Y.predict,
    m1.predict$fit[,1],includeSE=T)[2]),digits=3)
   if(is.na(rmspe.se)==FALSE){ rmspe.se.M1=rmspe.se }
   m1.rmspe.se = c(m1.rmspe.se,rmspe.se.M1)

   # fit model M0
   m0=fit(Y.train,X0.train) 

   m0.predict = fit_predict(Y.train,Y.predict,X0.train,X0.predict)
   # M0 root mean squared prediction error (RMSPE) 
   rmspe.M0 = signif(as.numeric(rmspe(Y.predict, m0.predict$fit[,1],
    includeSE=T)[1]),digits=3)
   m0.rmspe = c(m0.rmspe,rmspe.M0)

   # standard error of M0 RMSPE 
   rmspe.se.M0 = signif(as.numeric(rmspe(Y.predict,
    m0.predict$fit[,1],includeSE=T)[2]),digits=3)
    m0.rmspe.se = c(m0.rmspe.se,rmspe.se.M0)

   # RMSPE ratio
   rmspe.ratio = c(rmspe.ratio,signif(rmspe.M0/rmspe.M1,digits=3))

   # standard error of RMSPE ratio
   rmspe.ratio.se=c(rmspe.ratio.se,signif((1/rmspe.M1)*sqrt(rmspe.se.M0^2+
    ((1/rmspe.M1^2)*(rmspe.se.M1^2))),digits=3))

   # ratio of cor(predicted,observed) between M1 and M0
   cor.M0 = as.numeric(cor.test(m0.predict$fit[,1],Y.predict)$estimate)
   cor.M1 = as.numeric(cor.test(m1.predict$fit[,1],Y.predict)$estimate)  
   m0.cor=c(m0.cor, signif(cor.M0,digits=3))
   m1.cor=c(m1.cor, signif(cor.M1,digits=3))
   cor.ratio = c(cor.ratio,signif(cor.M1/cor.M0,digits=3))
   cor.df = c(cor.df,as.numeric(cor.test(m0.predict$fit[,1],
    Y.predict)$parameter))
  }
 }

 if(is.null(m1.rmspe)==FALSE & is.null(m0.rmspe)==FALSE){ 

  # mean of RMSPE over K folds and R repetitions
  m1.rmspe.avg = signif(mean(m1.rmspe[is.na(m1.rmspe)==FALSE]),digits=3)
  m0.rmspe.avg = signif(mean(m0.rmspe[is.na(m0.rmspe)==FALSE]),digits=3)

  # standard error of RMSPE mean over K folds and R repetitions
  m1.rmspe.avg.se = signif((1/length(m1.rmspe.se))*sqrt(sum(m1.rmspe.se^2)),
   digits=3)
  m0.rmspe.avg.se = signif((1/length(m0.rmspe.se))*sqrt(sum(m0.rmspe.se^2)),
   digits=3)

  # mean of RMSPE ratio over K folds and R repetitions
  rmspe.ratio.avg = signif(mean(rmspe.ratio[is.na(rmspe.ratio)==FALSE]),
   digits=3)
  # standard error of mean of RMSPE ratio over K folds and R repetitions
  rmspe.ratio.avg.se = NA
  rmspe.ratio.avg.se = signif((1/length(rmspe.ratio.se))*sqrt(sum(rmspe.ratio.se^2)),digits=3)

  # mean of cor(predicted, observed) ratios between M1 and M0
  cor.ratio.avg = signif(mean(cor.ratio[is.na(cor.ratio)==FALSE]),
   digits=3)

  out = list(m0.rmspe = m0.rmspe,m1.rmspe = m1.rmspe,
   m0.rmspe.se = m0.rmspe.se,m1.rmspe.se = m1.rmspe.se,
   m0.coef.n=length(summary(m0)$coefficients), 
   m1.coef.n=length(summary(m1)$coefficients),
   m0.cor=m0.cor,m1.cor=m1.cor,rmspe.ratio = rmspe.ratio,
   rmspe.ratio.se = rmspe.ratio.se,
   cor.ratio=cor.ratio,cor.df=cor.df,
   m0.rmspe.avg = m0.rmspe.avg,m1.rmspe.avg = m1.rmspe.avg,
   m0.rmspe.avg.se = m0.rmspe.avg.se,m1.rmspe.avg.se = m1.rmspe.avg.se,
   rmspe.ratio.avg = rmspe.ratio.avg,rmspe.ratio.avg.se = rmspe.ratio.avg.se,
   cor.ratio.avg=cor.ratio.avg,n.iter=length(m1.rmspe))
 }

 if(is.null(m1.rmspe)==TRUE | is.null(m0.rmspe)==TRUE){
  out=rep("NA",times=20)
 }
 return(out)
}


# this function should be encapsulated within each iteration of the other CV
cross_validation_with_penalty_tunel <- function (K,R,X0,X1,Y){
 # X,Y is X,Ytrain in the cross-validation for prediction error wo penalty
  
 rmspe.ratio=c();rmspe.ratio.se=c();
 m0.rmspe=c(); m0.rmspe.se=c(); m1.rmspe=c(); m1.rmspe.se=c();
 cvm=c(); cvsd=c(); nzero=c(); lambda.min=c();

 for(r in 1 : R){
  folds=cvFolds(length(Y), K = K, type = "random")
  tmp=cbind(folds$which,folds$subsets[,1])
  for(k in 1 : K){
   idx.predict=tmp[tmp[,1]==k,2]
   idx.train=tmp[tmp[,1]!=k,2]
   X0.train=X0[idx.train]
   if(is.matrix(X1)==T){X1.train=X1[idx.train,]}
   if(is.matrix(X1)==F){X1.train=X1[idx.train]}
   Y.train=Y[idx.train]
   X0.predict=X0[idx.predict]
   if(is.matrix(X1)==T){X1.predict=X1[idx.predict,]}
   if(is.matrix(X1)==F){X1.predict=X1[idx.predict]}
   Y.predict=Y[idx.predict]

   # imputation of N/A values in training sub-matrix
   Z1.train=t(as.matrix(cbind(Y.train,X1.train)))
   Z1.train.knn = t(impute.knn(Z1.train)$data)
   Y.train=Z1.train.knn[,1]
   X1.train=Z1.train.knn[,-1]

   # imputation of N/A values in prediction sub-matrix
   Z1.predict=t(as.matrix(cbind(Y.predict,X1.predict)))
   Z1.predict.knn = t(impute.knn(Z1.predict)$data)
   Y.predict=Z1.predict.knn[,1]
   X1.predict=Z1.predict.knn[,-1]

   # do k-fold cross-validation and return:
   # selected lambda (value of lambda that gives minimum cvm)
   # mean cross-validated error at selected lambda
   # estimate of standadrd error of cvm at selected lambda
   # number of non-zero coefficients at selected lambda 
  
   cvob=cv.glmnet(X1.train, Y.train, nfolds=3, alpha=0)
   #coef(cvob, s="lambda.min")
   lambda.min = c(lambda.min,cvob$lambda.min)
   cvm = c(cvm,cvob$cvm[match(cvob$lambda.min,cvob$lambda)])
   cvsd = c(cvsd,cvob$cvsd[match(cvob$lambda.min,cvob$lambda)])
   nzero = c(nzero,cvob$nzero[match(cvob$lambda.min,cvob$lambda)])

   # predict values from fitted glmnet object
   m1.predict = predict(cvob,newx=X1.predict, s="lambda.min", se.fit=T)

   if(!is.null(warnings())){
    write.table(paste(i,k,r,warnings(),sep=" "),
    file="rcv.knn.error.glmnet.out",append=T,quote=F,row.names=F,col.names=F)
   }

   # RMSPE for M1 (Root Mean Squared Prediction Error)
   rmspe.M1 = NA
   residuals=m1.predict-Y.predict
   rmspe = signif(as.numeric(sqrt(mean(residuals^2))),digits=3)
   if(is.na(rmspe)==FALSE){ rmspe.M1 = rmspe }
   m1.rmspe = c(m1.rmspe,rmspe.M1)
   
   # fit model M0
   m0=fit(Y.train,X0.train)

   m0.predict = fit_predict(Y.train,Y.predict,X0.train,X0.predict)
   # RMSPE for M0
   rmspe.M0 = signif(as.numeric(rmspe(Y.predict, m0.predict$fit[,1],
    includeSE=T)[1]),digits=3)
   m0.rmspe = c(m0.rmspe,rmspe.M0)

   # RMSPE ratio
   rmspe.ratio = c(rmspe.ratio,signif(rmspe.M0/rmspe.M1,digits=3))

  }
 }

 if(is.null(m1.rmspe)==FALSE & is.null(m0.rmspe)==FALSE){
  # mean of RMSPE over K folds and R repetitions
  m1.rmspe.avg = signif(mean(m1.rmspe[is.na(m1.rmspe)==FALSE]),digits=3)
  m0.rmspe.avg = signif(mean(m0.rmspe[is.na(m0.rmspe)==FALSE]),digits=3)

  # mean of RMSPE ratio over K folds and R repetitions
  rmspe.ratio.avg = signif(mean(rmspe.ratio[is.na(rmspe.ratio)==FALSE]),
   digits=3)

  out = list(m0.rmspe = m0.rmspe,m1.rmspe = m1.rmspe,
   m0.coef.n=length(summary(m0)$coefficients),m1.coef.n=nzero,
   cvm=cvm,cvsd=cvsd,lambda.min=lambda.min,
   rmspe.ratio = rmspe.ratio,
   m0.rmspe.avg = m0.rmspe.avg,m1.rmspe.avg = m1.rmspe.avg,
   rmspe.ratio.avg = rmspe.ratio.avg,n.iter=length(m1.rmspe))
 }

 if(is.null(m1.rmspe)==TRUE | is.null(m0.rmspe)==TRUE){
  out=rep("NA",times=12)
 }
 return(out)
 
}










