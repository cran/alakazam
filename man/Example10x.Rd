% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{Example10x}
\alias{Example10x}
\title{Small example 10x Genomics Ig V(D)J sequences from CD19+ B cells isolated from PBMCs of a healthy 
human donor. Down-sampled from data provided by 10x Genomics under a Creative Commons Attribute license,
and processed with their Cell Ranger pipeline.}
\format{
A data.frame with the following AIRR style columns:
  \itemize{
    \item  \code{sequence_id}:                Sequence identifier
    \item  \code{sequence_alignment}:         IMGT-gapped observed sequence.
    \item  \code{germline_alignment}:         IMGT-gapped germline sequence.
    \item  \code{v_call}:                     V region allele assignments.
    \item  \code{d_call}:                     D region allele assignments.
    \item  \code{j_call}:                     J region allele assignments.
    \item  \code{c_call}:                     Isotype (C region) assignment.
    \item  \code{junction}:                   Junction region sequence.
    \item  \code{junction_length}:            Length of the junction region in nucleotides.
    \item  \code{np1_length}:                 Combined length of the N and P regions proximal
                                              to the V region.
    \item  \code{np2_length}:                 Combined length of the N and P regions proximal
                                              to the J region.
    \item  \code{umi_count}:                  Number of unique molecular identifies atttributed to sequence.
    \item  \code{cell_id}:                    Cell identifier.
    \item  \code{locus}:                      Genomic locus of sequence.
}
}
\usage{
Example10x
}
\description{
Small example 10x Genomics Ig V(D)J sequences from CD19+ B cells isolated from PBMCs of a healthy 
human donor. Down-sampled from data provided by 10x Genomics under a Creative Commons Attribute license,
and processed with their Cell Ranger pipeline.
}
\references{
\enumerate{
  \item  Data source: https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_cd19_b
  \item  License: https://creativecommons.org/licenses/by/4.0/
}
}
\keyword{datasets}
