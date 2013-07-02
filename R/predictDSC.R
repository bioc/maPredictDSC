predictDSC=function(ano,celfile.path,annotation,preprocs=c("rma","gcrma","mas5"),
filters=c("mttest","ttest","wilcox"),classifiers=c("LDA","kNN","svm"),FCT=1.0,
CVP=4,NF=10,by=ifelse(NF>10,2,1), NR=5)
{

if(!all(preprocs%in%c("rma","gcrma","mas5"))){stop("preprocs must be rma,gcrma and/or mas5 !")}
if(!all(filters%in%c("mttest","ttest","wilcox"))){stop("filters must be mttest,ttest and/or wilcox !")}
if(!all(classifiers%in%c("LDA","kNN","svm"))){stop("classifiers must be LDA,kNN and/or svm !")}
if(!(FCT>=1.0 & FCT<1000)){stop("Fold change threshold must be between 1 and 1000")}

tab=table(ano$group)
tab=tab[names(tab)!="Test"]
if(!(length(tab)==2)){stop("The column group in the ano dataframe must have at most two classes (binary phenotype) and should include some Test samples!")}

smallgrpsize=min(tab)

if(!(CVP>1 & CVP<=smallgrpsize)){stop(paste("CVP must be integer >1 and at most ",smallgrpsize))}
if(!(NF<sum(tab))){warning("The maximum number of features is larger than the number of samples!")}
if(!(NR>=1)){stop("The number of repetitions NR needs to be a positive integer!")}



#clean up saved normalized data if already there
for(mm in preprocs){
 if(file.exists(paste(mm,"E07tmp.RData",sep=""))){
 system(paste("rm ",mm,"E07tmp.RData",sep=""))
}
}
modelsDSC=list()
results=list()
for(pre in preprocs){
for(classif in classifiers){
for(filt in filters){
out=maPredictDSC(ano=ano,celfile.path=celfile.path,annotation=annotation,
preproc.m=pre,filter.m=filt,classifier.m=classif,FCT=FCT,CVP=CVP,NF=NF,by=by,NR=NR,
otherCovariates=NULL)

comb=paste(pre,filt,classif,sep="_")
modelsDSC[[comb]]<-out
cat(comb);cat("\t");
cat("\n");

}
}
}
modelsDSC

}