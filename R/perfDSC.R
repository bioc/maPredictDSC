                                       

perfDSC=function(pred,gs){

require(ROCR)
require(ROC)
gs=gs[rownames(pred),]

########## CCEM
trueclass=max.col(gs)
predclass=max.col(pred)

z=sum(apply(rbind(pred[trueclass==predclass,]),1,max))-sum(apply(rbind(pred[trueclass!=predclass,]),1,max))
CCEM= (z/dim(gs)[1] + 1)/2



#BCM
delta=diag(max(trueclass))
b=matrix(0,max(trueclass),max(trueclass))
for (k in 1:max(trueclass)){
for (j in 1:max(trueclass)){
b[k,j]=abs(mean(pred[trueclass==k,j])-delta[k,j])
}
}

bcm=sum(diag(b))
BCM=1-bcm/max(trueclass)
###################################


#prec recall
#prec recall
AUPR=NULL
AUROC=NULL
for (myc in 1: dim(gs)[2]){
x=pred[,myc]
tc=gs[,myc][order(x)]
x=x[order(x)]

pred1 <- prediction(x, tc)
perf <- performance(pred1, "prec", "rec")
dat= cbind(recall= unlist(perf@x.values),precision= unlist(perf@y.values))
#dat[dat[,"recall"]==0,"precision"]<-0
dat= na.omit(dat)
recall= dat[,"recall"]
precision= dat[,"precision"]
AUPR=c(AUPR,trapezint(recall,precision,0,1))

roc1 <- rocdemo.sca(truth =tc , data = x, rule = dxrule.sca)
AUROC=c(AUROC,AUC(roc1))
}
AUROC_avg=mean(AUROC,na.rm=TRUE)
AUPR_avg=mean(AUPR,na.rm=TRUE)
#addedl later
if(BCM==1 & CCEM==1 & is.na(AUPR_avg)){
  AUPR_avg=1
}
if(is.nan(AUPR_avg)&length(unique(predclass)==1)){
  AUPR_avg=0.5
}



z=c(BCM,CCEM,AUPR_avg,AUROC_avg)
names(z)=c("BCM","CCEM","AUPR","AUC")
 return(z)
}

