\name{maPredictDSC-internal}
\alias{maPredictDSC}
\title{Internal maPredictDSC functions}
\description{
  Internal maPredictDSC functions. \code{maPredictDSC} is the main function of the package
  called via predictDSC}
\usage{
maPredictDSC(ano,celfile.path,annotation,preproc.m="rma",
filter.m="mttest",FCT=1.0,classifier.m="LDA", otherCovariates=NULL,CVP=4,NF=20,by=ifelse(NF>10,2,1), NR=10)
}
\details{
  These are not to be called directly by the user.
}
\keyword{internal}
