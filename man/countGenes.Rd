% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gene.R
\name{countGenes}
\alias{countGenes}
\title{Tabulates V(D)J allele, gene or family usage.}
\usage{
countGenes(
  data,
  gene,
  groups = NULL,
  copy = NULL,
  clone = NULL,
  fill = FALSE,
  mode = c("gene", "allele", "family", "asis"),
  remove_na = TRUE
)
}
\arguments{
\item{data}{data.frame with AIRR-format or Change-O style columns.}

\item{gene}{column containing allele assignments. Only the first allele in the
column will be considered when \code{mode} is "gene", "family" or 
"allele". The value will be used as it is with \code{mode="asis"}.}

\item{groups}{columns containing grouping variables. If \code{NULL} do not group.}

\item{copy}{name of the \code{data} column containing copy numbers for each 
sequence. If this value is specified, then total copy abundance
is determined by the sum of copy numbers within each gene.
This argument is ignored if \code{clone} is specified.}

\item{clone}{name of the \code{data} column containing clone identifiers for each 
sequence. If this value is specified, then one gene will be considered 
for each clone. Note, this is accomplished by using the most 
common gene within each \code{clone} identifier. As such,
ambiguous alleles within a clone will not be accurately represented.}

\item{fill}{logical of \code{c(TRUE, FALSE)} specifying when if groups (when specified)
lacking a particular gene should be counted as 0 if TRUE or not (omitted)}

\item{mode}{one of \code{c("gene", "family", "allele", "asis")} defining
the degree of specificity regarding allele calls. Determines whether 
to return counts for genes (calling \code{getGene}), 
families (calling \code{getFamily}), alleles (calling 
\code{getAllele}) or using the value as it is in the column
\code{gene}, without any processing.}

\item{remove_na}{removes rows with \code{NA} values in the gene column if \code{TRUE} and issues a warning. 
Otherwise, keeps those rows and considers \code{NA} as a gene in the final counts 
and relative abundances.}
}
\value{
A data.frame summarizing family, gene or allele counts and frequencies 
          with columns:
          \itemize{
            \item \code{gene}:         name of the family, gene or allele.
            \item \code{seq_count}:    total number of sequences for the gene.
            \item \code{seq_freq}:     frequency of the gene as a fraction of the total
                                       number of sequences within each grouping.
            \item \code{copy_count}:   sum of the copy counts in the \code{copy} column.
                                       for each gene. Only present if the \code{copy} 
                                       argument is specified.
            \item \code{copy_freq}:    frequency of the gene as a fraction of the total
                                       copy number within each group. Only present if 
                                       the \code{copy} argument is specified.
            \item \code{clone_count}:  total number of clones for the gene. Only present if 
                                       the \code{clone} argument is specified.
            \item \code{clone_freq}:   frequency of the gene as a fraction of the total
                                       number of clones within each grouping. Only present if 
                                       the \code{clone} argument is specified.
          }
          Additional columns defined by the \code{groups} argument will also be present.
}
\description{
Determines the count and relative abundance of V(D)J alleles, genes or families within
groups.
}
\examples{
# Without copy numbers
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="family")
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="gene")
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="allele")

# With copy numbers and multiple groups
genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
                    copy="duplicate_count", mode="family")

# Count by clone
genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
                    clone="clone_id", mode="family")

# Count absent genes 
genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", 
                    mode="allele", fill=TRUE)

}
