% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{seqEqual}
\alias{seqEqual}
\title{Test DNA sequences for equality.}
\usage{
seqEqual(seq1, seq2, ignore = as.character(c("N", "-", ".", "?")))
}
\arguments{
\item{seq1}{character string containing a DNA sequence.}

\item{seq2}{character string containing a DNA sequence.}

\item{ignore}{vector of characters to ignore when testing for equality.
Default is to ignore c("N",".","-","?")}
}
\value{
Returns \code{TRUE} if sequences are equal and \code{FALSE} if they are not.
          Sequences of unequal length will always return \code{FALSE} regardless of
          their character values.
}
\description{
\code{seqEqual} checks if two DNA sequences are identical.
}
\examples{
# Ignore gaps
seqEqual("ATG-C", "AT--C")
seqEqual("ATGGC", "ATGGN")
seqEqual("AT--T", "ATGGC")

# Ignore only Ns
seqEqual("ATG-C", "AT--C", ignore="N")
seqEqual("ATGGC", "ATGGN", ignore="N")
seqEqual("AT--T", "ATGGC", ignore="N")

}
\seealso{
Used by \link{pairwiseEqual} within \link{collapseDuplicates}.
          See \link{seqDist} for calculation Hamming distances between sequences.
}
