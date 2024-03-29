% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Diversity.R
\name{plotAbundanceCurve}
\alias{plotAbundanceCurve}
\title{Plot a clonal abundance distribution}
\usage{
plotAbundanceCurve(
  data,
  colors = NULL,
  main_title = "Rank Abundance",
  legend_title = NULL,
  xlim = NULL,
  ylim = NULL,
  annotate = c("none", "depth"),
  silent = FALSE,
  ...
)
}
\arguments{
\item{data}{\link{AbundanceCurve} object returned by \link{estimateAbundance}.}

\item{colors}{named character vector whose names are values in the 
\code{group} column of \code{data} and whose values are 
colors to assign to those group names.}

\item{main_title}{string specifying the plot title.}

\item{legend_title}{string specifying the legend title.}

\item{xlim}{numeric vector of two values specifying the 
\code{c(lower, upper)} x-axis limits. The lower x-axis 
value must be >=1.}

\item{ylim}{numeric vector of two values specifying the 
\code{c(lower, upper)} y-axis limits. The limits on the 
abundance values are expressed as fractions of 1: use
c(0,1) to set the lower and upper limits to 0\% and 100\%.}

\item{annotate}{string defining whether to added values to the group labels 
of the legend. When \code{"none"} (default) is specified no
annotations are added. Specifying (\code{"depth"}) adds 
sequence counts to the labels.}

\item{silent}{if \code{TRUE} do not draw the plot and just return the ggplot2 
object; if \code{FALSE} draw the plot.}

\item{...}{additional arguments to pass to ggplot2::theme.}
}
\value{
A \code{ggplot} object defining the plot.
}
\description{
\code{plotAbundanceCurve} plots the results from estimating the complete clonal 
relative abundance distribution. The distribution is plotted as a log rank abundance 
distribution.
}
\examples{
# Estimate abundance by sample and plot
abund <- estimateAbundance(ExampleDb, group="sample_id", nboot=100)
plotAbundanceCurve(abund, legend_title="Sample")

}
\seealso{
See \link{AbundanceCurve} for the input object and \link{estimateAbundance} for
generating the input abundance distribution. Plotting is performed with \link{ggplot}.
}
