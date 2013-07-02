maPredictDSC=function(ano,celfile.path,annotation,preproc.m="rma",
filter.m="mttest",FCT=1.0,classifier.m="LDA", otherCovariates=NULL,CVP=4,NF=20,by=ifelse(NF>10,2,1), NR=10){

require(affy)
require(limma)
require(gcrma)
require(ROC)
require(class)
require(e1071)
require(caret)
require(annotation,character.only=TRUE)

if("multicore"%in%installed.packages()){
mc=(require(multicore))
}else{mc=FALSE}


anpack=paste(unlist(strsplit(annotation,split=".db")),"ENTREZID",sep="")



#set1
if(!file.exists(paste(preproc.m,"E07tmp.RData",sep=""))){
abatch<-ReadAffy(filenames=ano$files,celfile.path=celfile.path)
if(preproc.m=="rma"){
esetall<-rma(abatch)
}
if(preproc.m=="mas5"){
esetall<-mas5(abatch)
}
if(preproc.m=="gcrma"){
esetall<-gcrma(abatch)
}


anoo=ano[ano$group!="Test",]
anop=ano[ano$group=="Test",]


esetall<-exprs(esetall)
esetp<-exprs(mas5calls(abatch))
presentmin=dim(anoo)[1]/4
oks<-apply(esetp[,anoo$files]=="P",1,sum)>=presentmin
esetall<-esetall[oks,]

nms=unlist(as.list(get(anpack)[rownames(esetall)]))
esetall=esetall[!is.na(nms),]

rownames(esetall)<-paste("F",rownames(esetall),sep="")
rm(esetp)
rm(abatch)


eset=esetall[,anoo$files]
peset=esetall[,anop$files]
rm(esetall)
save(anoo,anop,peset,eset,file=paste(preproc.m,"E07tmp.RData",sep=""))
}else{load(file=paste(preproc.m,"E07tmp.RData",sep=""))}

peset=peset[substr(rownames(peset),1,3)!="FAF",]
eset=eset[substr(rownames(eset),1,3)!="FAF",]

an=ano
grp<-as.character(an$group)
group<-grp

grp=factor(grp);
gi=sort(unique(ano$group))
ginotest=setdiff(gi,"Test")
contr<-paste(ginotest[2],"-",ginotest[1],sep="")


an=anoo

features=rownames(eset)
mydat=t(eset[features,])
mydat=as.data.frame(mydat)
mydat=cbind(mydat,anoo[,otherCovariates])
mydat$CLS=ifelse(an$group==ginotest[1],0,1)
TS=factor(an$group)


set.seed(100)



anoo$group=factor(anoo$group)

xxx=as.list(get(anpack))

sumall=NULL
#pb <- txtProgressBar(min = 2, max = NF, style = 3)
for(NN in seq(2,NF,by=by)){

poolres=NULL

ff=function(x){paste(x,collapse="+")}

for(iter in 1:NR){


resall<-NULL
varall<-NULL


foldlist=createFolds(factor(mydat$CLS), k = CVP, list = TRUE, returnTrain = TRUE)

mff=function(i){
 train <- foldlist[[i]]
 test=setdiff(1:(length(mydat$CLS)),train)


if(filter.m=="mttest"){
design <- model.matrix(formula(paste("~0",paste(c("group",otherCovariates),collapse="+"),sep="+")),data=anoo[train,])
colnames(design)[1:length(levels(TS))]<-c(levels(TS))

#limma
fit <- lmFit(eset[,train], design)
cont.matrix <- makeContrasts(contrasts=contr,levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
j=1
aT1<-topTable(fit2,coef=j, number=60000)
aT1$ID=rownames(aT1)
}


if(filter.m=="ttest"){
design <- model.matrix(formula(paste(paste(c("~group",otherCovariates),collapse="+"),sep="+")),data=anoo[train,])
#ordinary t
fit <- lmFit(eset[,train], design)
tstat.ord <- fit$coef / fit$stdev.unscaled / fit$sigma
aT1<-data.frame(ID=rownames(tstat.ord),logFC=fit$coef[,2],P.Value=pt(tstat.ord[,2],df=fit$df.residual),stringsAsFactors=FALSE)
aT1=aT1[order(aT1$P.Value),]
aT1$ID=rownames(aT1)
}


if(filter.m=="wilcox"){
#ordinary t
design <- model.matrix(formula(paste(paste(c("~group",otherCovariates),collapse="+"),sep="+")),data=anoo[train,])
fit <- lmFit(eset[,train], design)
tstat.ord <- fit$coef / fit$stdev.unscaled / fit$sigma
wf=function(x){wilcox.test(x[TS[train]==levels(TS[train])[1]],x[TS[train]==levels(TS[train])[2]])$p.value}
zx=apply(eset,1,wf)
aT1<-data.frame(ID=rownames(tstat.ord),logFC=fit$coef[,2],P.Value=zx,stringsAsFactors=FALSE)
aT1=aT1[order(aT1$P.Value),]
}
#
aT1$ENTREZ=unlist((xxx[substr(aT1$ID,2,100)]))
aT1=aT1[!duplicated(aT1$ENTREZ),]

#if(sortBy=="FC"){
#tg1<-aT1[abs(aT1$logFC)>log2(cutoff),]
#}else{tg1=aT1}

if(FCT>1){
tg1<-aT1[abs(aT1$logFC)>log2(FCT),]
}else{
tg1<-aT1
}


tg1<-tg1[1:min(NF,dim(tg1)[1]),]
#cat(dim(tg1)[1]);cat("\n");
features=tg1$ID

vr=ff(c(features[1:NN],otherCovariates))

if(classifier.m=="LDA"){
#lda
mod1 <- lda(formula=formula(paste("CLS~1",vr,sep="+")), data = mydat[train,],prior = c(1,1)/2)
mypre=predict(mod1,mydat[test,])$posterior[,2]
predl=as.numeric(as.vector(predict(mod1,mydat[test,])$class))
}

#3-nn
if(classifier.m=="kNN"){
mod1<-knn(mydat[train,c(features[1:NN],otherCovariates)], mydat[test,c(features[1:NN],otherCovariates)], mydat[train,"CLS"], k = 3, prob=TRUE)
attr(mod1,"prob")<-NULL
mypre=as.numeric(as.character(mod1))
predl=as.numeric(as.character(mod1))
}


#svm
if(classifier.m=="svm"){
mod1 <- svm(mydat[train,c(features[1:NN],otherCovariates)], factor(mydat[train,"CLS"]))
mypre <- as.numeric(as.character(predict(mod1, mydat[test,c(features[1:NN],otherCovariates)])))
predl<-mypre
}

if(length(unique(mypre))==1 | length(unique(mydat[test,"CLS"]))==1){ve1=mean(mydat[test,"CLS"]==round(mypre))}else{
roc1 <- rocdemo.sca(truth = mydat[test,"CLS"], data = mypre, rule = dxrule.sca)
ve1= AUC(roc1)
}

ve1

}#end function

ml=list()
for(mk in 1:CVP){
ml[[mk]]<-mk
}

if(mc){
resall=unlist(mclapply(ml,mff))
}else{
resall=unlist(lapply(ml,mff))
}

poolres=c(poolres,resall)
}



sumall=rbind(sumall,c(NN,mean(poolres),sd(poolres)))
#setTxtProgressBar(pb, NN)
}


colnames(sumall)<-c("NN","meanAUC","sdAUC")
#NN=sumall[order(sumall[,"meanAUC"]/(sumall[,"sdAUC"]+0.01),decreasing=TRUE)[1],"NN"]
#best_tAUC=max(sumall[,"meanAUC"]/(sumall[,"sdAUC"]+0.01))
NN=sumall[order(sumall[,"meanAUC"],decreasing=TRUE)[1],"NN"]
best_AUC=max(sumall[,"meanAUC"])

#cat(NN)
#cat("\n");

if(filter.m=="mttest"){
#limma
design <- model.matrix(formula(paste("~0",paste(c("group",otherCovariates),collapse="+"),sep="+")),data=anoo)
colnames(design)[1:length(levels(TS))]<-c(levels(TS))
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(contrasts=contr,levels=design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
j=1
aT1<-topTable(fit2,coef=j, number=60000)
aT1$ID=rownames(aT1)
}

if(filter.m=="ttest"){
#ordinary t
design <- model.matrix(formula(paste(paste(c("~group",otherCovariates),collapse="+"),sep="+")),data=anoo)
fit <- lmFit(eset, design)
tstat.ord <- fit$coef / fit$stdev.unscaled / fit$sigma
aT1<-data.frame(ID=rownames(tstat.ord),logFC=fit$coef[,2],P.Value=pt(tstat.ord[,2],df=fit$df.residual),stringsAsFactors=FALSE)
aT1=aT1[order(aT1$P.Value),]
aT1$ID=rownames(aT1)
}


if(filter.m=="wilcox"){
#ordinary t
design <- model.matrix(formula(paste(paste(c("~group",otherCovariates),collapse="+"),sep="+")),data=anoo)
fit <- lmFit(eset, design)
tstat.ord <- fit$coef / fit$stdev.unscaled / fit$sigma
wf=function(x){wilcox.test(x[TS==levels(TS)[1]],x[TS==levels(TS)[2]])$p.value}
zx=apply(eset,1,wf)
aT1<-data.frame(ID=rownames(tstat.ord),logFC=fit$coef[,2],P.Value=zx,stringsAsFactors=FALSE)
aT1=aT1[order(aT1$P.Value),]
}


aT1$ENTREZ=unlist((xxx[substr(aT1$ID,2,100)]))
aT1=aT1[!duplicated(aT1$ENTREZ),]

if(FCT>1){
tg1<-aT1[abs(aT1$logFC)>log2(FCT),]
}else{
tg1<-aT1
}




tg1<-tg1[1:min(NF,dim(tg1)[1]),]
#cat(dim(tg1)[1]);cat("\n");
features=tg1$ID[1:NN]

vr=ff(c(features,otherCovariates))
mydatp=data.frame(t(peset[features,]),anop[,otherCovariates])

if(classifier.m=="LDA"){
mod1 <- lda(formula=formula(paste("CLS~1",vr,sep="+")), data = mydat,prior = c(1,1)/2)
mypre=predict(mod1,mydatp[,])$posterior[,2]
}

if(classifier.m=="kNN"){
mod1<-knn(mydat[,c(features,otherCovariates)], mydatp[,c(features,otherCovariates)], mydat[,"CLS"], k = 3, prob=TRUE)
attr(mod1,"prob")<-NULL
mypre=as.numeric(as.character(mod1))
}


#svm
if(classifier.m=="svm"){
mod1 <- svm(mydat[,c(features,otherCovariates)], factor(mydat[,"CLS"]))
mypre <- as.numeric(as.character(predict(mod1, mydatp[,c(features,otherCovariates)])))
}




out=data.frame(round(mypre,4),1-round(mypre,4))
#rownames(out)<-unlist(strsplit(rownames(mydatp),split=".cel"))
rownames(out)<-rownames(mydatp)
colnames(out)<-ginotest[2:1]

 
return(list(predictions=out,features=features,model=mod1,performanceTr=sumall,best_AUC=best_AUC))
} 