\name{perfDSC}
\alias{perfDSC}
\title{Area Under the Precision-Recall Curve (AUPR), Belief Confusion Metric (BCM) and Correct Class Enrichment Metric (CCEM).}
\description{
This function implements the three metrics used in the IMPROVER Diagnostic Signature Challenge.  
}
\usage{
perfDSC(pred,gs)
}
\arguments{
  \item{pred}{A belief matrix, with rows coresponding to samples and columns to classes. The values
  are between 0 and 1 and sum on each row is 1. It needs to have row names. The belief values
  are the result of a prediction made by a model.}
  \item{gs}{A matrix, with rows coresponding to samples and columns to classes that give the 
  true (gold standard) class membership of samples.}
 



}

\details{
See cited documents for more details.
}
\value{
 A named vector that includes the BCM, CCEM, AUPR_avg and Accuracy.}
  
  


\references{
Adi L. Tarca, Mario Lauria, Michael Unger, Erhan Bilal, Stephanie Boue, Kushal Kumar Dey,
 Julia Hoeng, Heinz Koeppl, Florian Martin, Pablo Meyer, Preetam Nandy, Raquel Norel,
  Manuel Peitsch, Jeremy J Rice, Roberto Romero, Gustavo Stolovitzky, Marja Talikka,
   Yang Xiang, Christoph Zechner, and IMPROVER DSC Collaborators, 
   Strengths and limitations of microarray-based phenotype prediction: Lessons learned 
   from the IMPROVER Diagnostic Signature Challenge. Bioinformatics, submitted 2013. \cr

}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu>}

\seealso{\code{\link{predictDSC}}}

\examples{
#asume a 3 class classification problem; gs is the gold standard and pred are predictions
gs=cbind(A=c(1,1,1,1,0,0,0,0,0,0,0,0),B=c(0,0,0,0,1,1,1,1,0,0,0,0),C=c(0,0,0,0,0,0,0,0,1,1,1,1))
rownames(gs)<-paste("sample",1:12,sep="")
pred=cbind(A=c(0.6,0.9,1,0.3,0,0,0,0,0,0,0,0),B=c(0.4,0.1,0,0.7,1,1,0.7,1,0,0,0,0),C=c(0,0,0,0,0,0,0.3,0,1,1,1,1))
rownames(pred)<-paste("sample",1:12,sep="")
#male sure the sum per row is 1 is both gs and pred
apply(gs,1,sum)
apply(pred,1,sum)
#compute perfromance
perfDSC(pred,gs)



}

\keyword{parametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
