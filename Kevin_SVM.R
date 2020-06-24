library(ggplot2)
library(lattice)
library(e1071)
library(cvAUC)
library(ROCR)
library(caret)
library(pROC)
library(tidyverse)
library(MASS)
rocplot =function (pred , truth , ...){
  # pROC_obj <- roc(pred, truth, ci=TRUE, ci.alpha=0.95, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, print.auc=TRUE)
  # sens.ci <- ci.see(pROC_obj)
  # plot(sens.ci, type="shape", col='lightblue')
  # plot(sens.ci, type="bars")}
  predob = prediction (pred , truth)
  perf = performance (predob , "tpr", "fpr")
  #return(perf)
  #plot(perf, col="grey82", lty=1)
  #plot(perf ,... , spread.scale = 2, colorize=TRUE)
  plot(perf , ..., lwd=2, avg="threshold",  main="ROC - SVM", col="purple", lty = 1)
  xx <- c(0, 0.1, 1.0)
  yy <- c(0,0.1,1.0)
  lines(xx,yy, col="red", lty = 2)}
  

data1 = read.csv('/media/admin/6FD0-88CC/Lauren/WSS_data/Aneurisk_MDrive_MCAcases_Rtable_2WSS_LSA_Threshold.csv', header=T)
data1 = data1[,-1]
Status = data1$ruptureStatus 
levels(Status) <- c(0,1)
#IA_location = data1$IA_Location
#levels(IA_location) <- c(0,1)
Type = data1$IA_Type
levels(Type) <- c(0,1)
rand <- sample(nrow(data1))
data1[rand, ]

#data splitting to keep the ratio of unblanced between U and R
set.seed(1)
iter=100
Error = matrix(0, ncol=4, nrow=iter) #Save confusion matrix
TPR_FPR = matrix(0, ncol=2, nrow=iter) # matrix for true positive and false posiive rates
xout <- rep(as.numeric(0), 47)
yout <- rep(as.numeric(0), 47)
fitted <- vector(mode="list", iter)
labels <- vector(mode="list", iter)
AUCROC <- 0
auc_ROCR <- vector(mode='list', iter)
ci_ROCR<- vector(mode="list", iter)
CIROC <- matrix(0,1,2)
model <- glm(ruptureStatus~., data =data1, family = binomial) %>%stepAIC(trace = FALSE)
summary(model)
for(j in 1:iter){
  #rand <- sample(nrow(data1))
  #data1 <- data1[rand,]
  #split <- sample(1:nrow(data1), size = 0.2*nrow(data1))
  #Train = data1[-split,]
  #Testing data
  #Test =data1[split,]
  
  #Splitting data (0.8 and 0.2)
  index = createDataPartition(Status, p=0.8,list=FALSE, times=1)
   Train =data1[index,]
   Test =data1[-index,]
  
  #Splitting data into Training (50/50 10U and 10R aneurysms) and Testing(Remainder of Data)
  # datasplit <- split(data1, Status)
  # dU <- data.frame(datasplit[2])
  # colnames(dU) = colnames(data1)
  # dR <- data.frame(datasplit[1])
  # colnames(dR) = colnames(data1)
  # dR <- dR[sample(1:nrow(dR)),]
  # dU <- dU[sample(1:nrow(dU)),]
  # Train = data.frame(dU[1:10,])
  # Train <-rbind(Train, dR[1:10,])
  # 
  # Test = data.frame(dU[11:nrow(dU),])
  # Test <- rbind(Test, dR[11:nrow(dR),])
  


  tc <- tune.control(cross = 10)
  tune.out=tune(svm, ruptureStatus~ ., data=Train, kernel ='radial',ranges=list(cost=c(0.1,0.5,1:6),
                                                                               gamma=c(0.5,1,2,3,4) ),tunecontrol = tc)
  summary.tune = summary (tune.out)
  train_control <- trainControl(method="cv", number=10)
  #SVM
  svmfit=svm(ruptureStatus ~ ., data = Train, kernel ='radial',gamma=summary.tune$best.parameters$gamma,cost=summary.tune$best.parameters$cost)
  Error[j,]= c(table(true=Test[,"ruptureStatus"], pred=predict(tune.out$best.model,newdata =Test)))
  TPR_FPR[j,1] = Error[j,4]/(Error[j,4]+Error[j,2])
  TPR_FPR[j,2] = Error[j,3]/(Error[j,3]+Error[j,1])
  svmauc = svm(ruptureStatus~., data=Train,  kernel="radial", gamma=summary.tune$best.parameter$gamma, cost=summary.tune$best.parameters$cost, probability=TRUE)
  svmaucprob <- predict(svmauc, type="prob",  newdata=Test, probability=TRUE)
  svmaucroc <- prediction(attr(svmaucprob, "probabilities")[,2], Test["ruptureStatus"])
  fitted[[j]] =attributes (predict (svmfit ,Test, decision.values=TRUE))$decision.values
  labels[[j]] = Test$ruptureStatus
  predob = prediction (fitted[[j]] , labels[[j]])
  auc_ROCR <- performance(predob, measure = "auc")
  ci_ROCR <- ci.cvAUC(fitted[[j]], labels[[j]], confidence=0.8)$ci
  CIROC <- CIROC[] + ci_ROCR
  AUCROC <- AUCROC + auc_ROCR@y.values[[1]]
  print(predict(tune.out$best.model,newdata =Test))
  print(auc_ROCR@y.values)
#invisible(readline(prompt="Press [enter] to continue"))
}
 
meanAUC = AUCROC/iter
meanAUC
rocplot(fitted[], labels[])
legend("bottomright",legend=c("WSS, OSI, LSA AUC: 0.40","WSS, OSI AUC: 0.50", "WSS, LSA AUC: 0.42"), col=c("purple","black","green"),
       lty=c(1,1,1))
CIROC/iter
#testcv <- ci.cvAUC(fitted, labels, confidence=0.95)
#Average prediction rate
mean((Error[,1]+Error[,4])/rowSums(Error))
ctable <- as.table(matrix(c(sum(Error[,1]),sum(Error[,2]),sum(Error[,3]),sum(Error[,4])), nrow=2, byrow=TRUE, dimnames=list(c("True R","True U"), c("Ruptured","Unruptured"))))
fourfoldplot(ctable, conf.level = 0.95, margin=1, main = "Confusion Matrix: WSS, LSA")
RupAcc = (ctable[1,1]/(ctable[1,1]+ctable[2,1]))*100
UnrupAcc = (ctable[2,2]/(ctable[2,2]+ctable[1,2]))*100

RupAcc
UnrupAcc 


tpr = TP/(TP+FN)
fpr = FP/(FP+TN)
X <- c(0, tpr, 1)
Y <- c(0, fpr, 1)
caTools::trapz(Y,X)



