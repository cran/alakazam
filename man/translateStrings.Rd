% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Core.R
\name{translateStrings}
\alias{translateStrings}
\title{Translate a vector of strings}
\usage{
translateStrings(strings, translation)
}
\arguments{
\item{strings}{vector of character strings to modify.}

\item{translation}{named character vector or a list of character vectors specifying 
the strings to replace (values) and their replacements (names).}
}
\value{
A modified \code{strings} vector.
}
\description{
\code{translateStrings} modifies a character vector by substituting one or more 
strings with a replacement string.
}
\details{
Does not perform partial replacements. Each translation value must match a complete 
\code{strings} value or it will not be replaced.  Values that do not have a replacement
named in the \code{translation} parameter will not be modified.

Replacement is accomplished using \link{gsub}.
}
\examples{
# Using a vector translation
strings <- LETTERS[1:5]
translation <- c("POSITION1"="A", "POSITION5"="E")
translateStrings(strings, translation)

# Using a list translation
strings <- LETTERS[1:5]
translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
translateStrings(strings, translation)

}
\seealso{
See \link{gsub} for single value replacement in the base package.
}
