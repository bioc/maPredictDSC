\name{getLungCancerData}
\alias{getLungCancerData}
\title{Download the Lung Cancer Example CEL Files for the former LungCancerACvsSCCGEO data package}
\description{
This function downloads 30 .CEL file needed to illustrate predictive modeling starting with raw gene expression .CEL files.  
}
\usage{
getLungCancerData(ano)
}
\arguments{
  \item{ano}{A data frame with at least one columns: \code{files} that should have valid .CEL.gz file names from the
  zenodo record 'https://zenodo.org/api/records/23108090'}
}


\details{
See cited documents for more details.
}
\value{
 Returns the path to the location of the .CEL.gz files that were reterieved from zenodo.org. 
 The files are cached locally and can be read directly by
\code{\link[affy]{ReadAffy}}.
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

\seealso{\code{\link{maPredictDSC}}}

\examples{
library(maPredictDSC)
data(LungCancerACvsSCCGEO)
anoLC
gsLC
table(anoLC$group)

\dontrun{
celfile.path <- getLungCancerData(anoLC)
}
}
\keyword{parametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line