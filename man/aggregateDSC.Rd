\name{aggregateDSC}
\alias{aggregateDSC}
\title{Combine predictions from several fitted models fitted with predictDSC}
\description{
This function simply adds the posterior probabilities for a given class and sample from
several models, and scales the resulting sums so that the sum over the classes is 1.0. 
}
\usage{
aggregateDSC(modlist)
}
\arguments{
  \item{modlist}{An object returned by predictDSC.}
}

\details{
See cited documents for more details.
}
\value{
A a data frame with the predicted class membership belief value 
  (posterior probability) for each sample (row) and each class (column).
  
}
  
  


\references{
Adi L. Tarca, Mario Lauria, Michael Unger, Erhan Bilal, Stephanie Boue, Kushal Kumar Dey,
 Julia Hoeng, Heinz Koeppl, Florian Martin, Pablo Meyer, Preetam Nandy, Raquel Norel,
  Manuel Peitsch, Jeremy J Rice, Roberto Romero, Gustavo Stolovitzky, Marja Talikka,
   Yang Xiang, Christoph Zechner, and IMPROVER DSC Collaborators, 
   Strengths and limitations of microarray-based phenotype prediction: Lessons learned 
   from the IMPROVER Diagnostic Signature Challenge. Bioinformatics, submitted 2013. \cr

Tarca AL, Than NG, Romero R, Methodological Approach from the Best Overall 
 Team in the IMPROVER Diagnostic Signature Challenge, Systems Biomedicine, submitted, 2013. \cr  

}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu>}

\seealso{\code{\link{predictDSC}}}

\examples{
#see function predictDSC for example


}

\keyword{parametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
