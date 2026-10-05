\name{predictDSC}
\alias{predictDSC}
\title{Phenotype prediction using microarray data: approach of the best overall team 
in the IMPROVER Diagnostic Signature Challenge}
\description{
This function implements the classification pipeline of the best overall team (Team221) 
in the IMPROVER Diagnostic Signature Challenge. The function ofers also eploring other 
combinations of data preprocessing, feature selection and classifier types.  
}
\usage{
predictDSC(ano,celfile.path,annotation,preprocs=c("rma","gcrma","mas5"),
filters=c("mttest","ttest","wilcox"),classifiers=c("LDA","kNN","svm"),FCT=1.0, 
CVP=4,NF=10,by=ifelse(NF>10,2,1), NR=5)
}
\arguments{
  \item{ano}{A data frame with two columns: \code{files} and \code{group} giving the names of the Affymetrix
  .cel files (no full path) and their corresponding groups. Only two groups are allowed as well
  as a third group called "Test". The samples corresponding to these will not be used in 
  training but will be used to normalize the training data with.}

\item{celfile.path}{ The location of the directory where the .cel files are located.}

\item{annotation}{ The names of a package that can be used to map the probesets to the ENTREZ
gene IDS in order to deal with duplicate probesets pre gene. E.g.\code{hgu133plues2.db} }
\item{preprocs}{ A character vector giving the names of the normalization methods to try.
Supported options are "rma","gcrma","mas5"}
\item{filters}{ A character vector giving the names of the methods to use to rank features. 
Supported options are "mttest" for moderated t-test using limma package,"ttest" for regular t-test,
and "wilcox" for wilcoxon test.}
\item{classifiers}{ A character vector giving the names of the classifier types to use for
learning the relation between expression levels and phenotype. Supported options are 
"LDA","kNN","svm".}

\item{FCT}{ A numeric value giving the fold change threshokd to be used to filter out non-relevant
features. Note, setting it to a too large value can produce an error as there need to be at least 
\code{NF} probestes with a fold change larger than FCT in each fold of the cross-validation.}
\item{CVP}{ The number of cross-validation partitions to create (minimum is 2). Do use a CVP 
value which ensures that at least two samples from the smalest group are kept for testing 
at each fold. E.g. If you have 10 samples in the smalest of the 2 groups a CVP of 4 would be maximum.}
\item{NF}{ The maximum number of features that would make sense to consider using as predcitors
in the models. NF should be less than the number of training samples.}
\item{by}{ The size of the step when searching for the number of features to include. By default th esearch starts with the top 2 features, and 
a number of "by" features are added up to NF.}
\item{NR}{ An integer number between 1 and Inf giving the number of times the cross-validation 
should be repeated to ensure a robust solution to the question: how many features to use as predictors
in the model?.}

}

\details{
See cited documents for more details.
}
\value{
 A list object containing one item for each possible combination between the elements of \code{preprocs},
 \code{filters}, and \code{classifiers}. Each item of the list contains the following information:
  \code{predictions} - a data frame with the predicted class membership belief value 
  (posterior probability) for each sample (row) and each class (column).
  \code{features} - Names of the Affy probesets used as predictors by the model. A letter "F" is 
  added as suffix to the probeset names.
  \code{model} - A fitted model object as produced by the lda, svm and kNN functions.
  \code{performanceTr} - A matrix giving the number of features tested (NN) mean AUC over 
  all folds and repetitions (meanAUC), and the standard deviation of AUC values accross folds and 
  repeats of the cross-validation. 
  \code{bestAUC} - The value of mean AUC corresponding to the optimal number of features chosen. 
}
  
  


\references{
Tarca AL, Lauria M, Unger M, Bilal E, Boue S, Dey KK, Hoeng J, Koeppl H,
Martin F, Meyer P, Nandy P, Norel R, Peitsch M, Rice JJ, Romero R,
Stolovitzky G, Talikka M, Xiang Y, Zechner C, and IMPROVER DSC Collaborators.
Strengths and limitations of microarray-based phenotype prediction: lessons
learned from the IMPROVER Diagnostic Signature Challenge.
\emph{Bioinformatics}. 2013;29(22):2892-2899.
doi:10.1093/bioinformatics/btt492. \cr

Tarca AL, Than NG, Romero R.
Methodological approach from the Best Overall Team in the sbv IMPROVER
Diagnostic Signature Challenge.\emph{Systems Biomedicine}. 2013;1(4):217-227.
doi:10.4161/sysb.25980.
}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu>}

\seealso{\code{\link{aggregateDSC}}}

\examples{
library(maPredictDSC)
data(LungCancerACvsSCCGEO)
anoLC
gsLC
table(anoLC$group)

\dontrun{
celfile.path <- getLungCancerData(anoLC)

#run a series of methods combinations
modlist=predictDSC(ano=anoLC,
celfile.path=celfile.path,
annotation="hgu133plus2.db", preprocs=c("rma"),
filters=c("mttest","ttest"),classifiers=c("LDA","kNN"),
CVP=2,NF=4, NR=1,FCT=1.0)

#rank combinations by the performance on training data (AUC)
trainingAUC=sort(unlist(lapply(modlist,"[[","best_AUC")),decreasing=TRUE)
trainingAUC

#optional step; since we know the class of the test samples, let's see how the
#methods combinations perform on the test data

perfTest=function(out){
perfDSC(pred=out$predictions,gs=gsLC)
}
testPerf=t(data.frame(lapply(modlist,perfTest)))
testPerf=testPerf[order(testPerf[,"AUC"],decreasing=TRUE),]
testPerf

#aggregate predictions from top 3 combinations of methods
best3=names(trainingAUC)[1:3]
aggpred=aggregateDSC(modlist[best3])
#test the aggregated model on the test data
perfDSC(aggpred,gsLC)
}
}
\keyword{parametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line