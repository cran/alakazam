% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sequence.R
\name{getAAMatrix}
\alias{getAAMatrix}
\title{Build an AA distance matrix}
\usage{
getAAMatrix(gap = 0)
}
\arguments{
\item{gap}{value to assign to characters in the set \code{c("-", ".")}.}
}
\value{
A \code{matrix} of amino acid character distances with row and column names 
          indicating the character pair.
}
\description{
\code{getAAMatrix} returns a Hamming distance matrix for IUPAC ambiguous
amino acid characters.
}
\examples{
getAAMatrix()

}
\seealso{
Creates an amino acid distance matrix for \link{seqDist}.
          See \link{getDNAMatrix} for nucleotide distances.
}
