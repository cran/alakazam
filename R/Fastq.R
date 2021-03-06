#' Load sequencing quality scores from a FASTQ file
#' 
#' \code{readFastqDb} adds the sequencing quality scores to a data.frame
#' from a FASTQ file. Matching is done by `sequence_id`.
#'
#' @param    data            \code{data.frame} containing sequence data.
#' @param    fastq_file      path to the fastq file
#' @param    quality_offset  offset value to be used by ape::read.fastq. It is 
#'                           the value to be added to the quality scores 
#'                           (the default -33 applies to the Sanger format and 
#'                           should work for most recent FASTQ files).
#' @param    header          FASTQ file header format; one of \code{"presto"} or 
#'                           \code{"asis"}. Use \code{"presto"} to specify 
#'                           that the fastq file headers are using the pRESTO
#'                           format and can be parsed to extract 
#'                           the \code{sequence_id}. Use \code{"asis"} to skip 
#'                           any processing and use the sequence names as they are.
#' @param    sequence_id     column in \code{data} that contains sequence 
#'                           identifiers to be matched to sequence identifiers in 
#'                           \code{fastq_file}. 
#' @param    sequence        column in \code{data} that contains sequence data.
#' @param    sequence_alignment   column in \code{data} that contains IMGT aligned sequence data.      
#' @param    v_cigar         column in \code{data} that contains CIGAR 
#'                           strings for the V gene alignments.     
#' @param    d_cigar         column in \code{data} that contains CIGAR 
#'                           strings for the D gene alignments.   
#' @param    j_cigar         column in \code{data} that contains CIGAR 
#'                           strings for the J gene alignments.     
#' @param    np1_length      column in \code{data} that contains the number
#'                           of nucleotides between the V gene and first D gene 
#'                           alignments or between the V gene and J gene alignments.
#' @param    np2_length      column in \code{data} that contains the number
#'                           of nucleotides between either the first D gene and J 
#'                           gene alignments or the first D gene and second D gene
#'                           alignments.
#' @param    v_sequence_end  column in \code{data} that contains the 
#'                           end position of the V gene in \code{sequence}.
#' @param    d_sequence_end  column in \code{data} that contains the 
#'                           end position of the D gene in \code{sequence}.                      
#' @param    style           how the sequencing quality should be returned;
#'                           one of \code{"num"}, \code{"phred"}, or \code{"both"}.
#'                           Specify \code{"num"} to store the quality scores as strings of 
#'                           comma separated numeric values. Use \code{"phred"} to have
#'                           the function return the scores as Phred (ASCII) scores. 
#'                           Use \code{"both"} to retrieve both.
#' @param    quality_sequence     specify \code{TRUE} to keep the quality scores for 
#'                                \code{sequence}. If false, only the quality score
#'                                for \code{sequence_alignment} will be added to \code{data}.
#'                                
#' @return   Modified \code{data} with additional fields:
#'           \enumerate{
#'                 \item \code{quality_alignment}:     A character vector with ASCII Phred 
#'                                                     scores for \code{sequence_alignment}.
#'                 \item \code{quality_alignment_num}: A character vector, with comma separated 
#'                                                     numerical quality values for each 
#'                                                     position in \code{sequence_alignment}.
#'                 \item \code{quality}:      A character vector with ASCII Phred 
#'                                                     scores for \code{sequence}.
#'                 \item \code{quality_num}:  A character vector, with comma separated 
#'                                                     numerical quality values for each 
#'                                                     position in \code{sequence}.
#'           }
#' @seealso \link{maskPositionsByQuality} and \link{getPositionQuality}
#' 
#' @examples
#' db <- airr::read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
#' fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
#' db <- readFastqDb(db, fastq_file, quality_offset=-33)
#'
#' @export
readFastqDb <- function(data, fastq_file, quality_offset=-33, 
                        header=c("presto", "asis"), 
                        sequence_id="sequence_id",
                        sequence="sequence",
                        sequence_alignment="sequence_alignment",
                        v_cigar="v_cigar",
                        d_cigar="d_cigar",
                        j_cigar="j_cigar",
                        np1_length="np1_length",
                        np2_length="np2_length",
                        v_sequence_end="v_sequence_end",
                        d_sequence_end="d_sequence_end",
                        style=c("num", "ascii", "both"),
                        quality_sequence=FALSE) {
   
   check_cols <- c(sequence_id, 
                   sequence, 
                   sequence_alignment, 
                   v_cigar, d_cigar, j_cigar, 
                   np1_length, np2_length, 
                   v_sequence_end, d_sequence_end)
   
   alakazam::checkColumns(data, check_cols)
   
   style <- match.arg(style)
   
   # Process the fastq file
   header <- match.arg(header)
   fastq <- ape::read.fastq(fastq_file, offset=quality_offset) #default: -33 (pRESTO)
   fastq_db <- data.frame(
      "quality_num"=as.vector(sapply(attr(fastq, "QUAL"), paste0, collapse=",")),
      stringsAsFactors = F)
   
   fastq_db$quality <- sapply(fastq_db[["quality_num"]], 
                                     function(qual, quality_offset) {
                                        paste0(sapply(strsplit(qual, ",")[[1]], function(x,quality_offset) {
                                           y <- as.numeric(x) - quality_offset
                                           rawToChar(as.raw(y))
                                        },quality_offset), sep="",collapse="")
                                     }, quality_offset)
   
   fastq_db[[sequence_id]] <- attr(fastq, "names")
   if (header=="presto") {
      fastq_db[[sequence_id]] <- gsub("\\|.+","",attr(fastq, "names"))
   }

   # Merge
   by <- sequence_id
   names(by) <- sequence_id
   data <- data %>%
      left_join(fastq_db, by=by)
   
   data <- sequenceAlignmentQuality(data,
                                  sequence=sequence,
                                  sequence_id=sequence_id,
                                  sequence_alignment=sequence_alignment,
                                  quality="quality",
                                  quality_num="quality_num",
                                  v_cigar=v_cigar,
                                  d_cigar=d_cigar,
                                  j_cigar=j_cigar,
                                  np1_length=np1_length,
                                  np2_length=np2_length,
                                  v_sequence_end=v_sequence_end,
                                  d_sequence_end=d_sequence_end,
                                  raw=FALSE)
   
   if (!quality_sequence) {
      data[['quality']] <- NULL
      data[['quality_num']] <- NULL
   }
   
   if (style != "both") {
      if (style == "phred") {
         data[['quality_alignment_num']] <- NULL
      } else {
         data[['quality_alignment']] <- NULL
      }
   }
   
   data 
}

# Thanks!:
# https://drive5.com/usearch/manual/cigar.html &
# https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
# M	 	Match (alignment column containing two letters). This could contain two
#        letters (mismatch) or two identical letters. USEARCH generates CIGAR strings
#        containing Ms rather than X's and ='s (see below).
# N	   Alignment gap 	Next x positions on ref don’t match (Deletion in query?)
# D	 	Deletion (gap in the target sequence).
# I	 	Insertion (gap in the query sequence).
# S	 	Segment of the query sequence that does not appear in the alignment.
#        This is used with soft clipping, where the full-length query sequence
#        is given (field 10 in the SAM record). In this case, S operations specify
#        segments at the start and/or end of the query that do not appear in a
#        local alignment.
# H	 	Segment of the query sequence that does not appear in the alignment.
#        This is used with hard clipping, where only the aligned segment of the
#        query sequences is given (field 10 in the SAM record). In this case, H
#        operations specify segments at the start and/or end of the query that
#        do not appear in the SAM record.
# =	 	Alignment column containing two identical letters. USEARCH can read
#        CIGAR strings using this operation, but does not generate them.
# X	 	Alignment column containing a mismatch, i.e. two different letters.
#        USEARCH can read CIGAR strings using this operation, but does not generate them.

# Tested with IgBlast output, not with IMGT
calcSequenceAlignmentQuality <- function(sequence_db,
                                        sequence="sequence",
                                        sequence_id="sequence_id",
                                        sequence_alignment="sequence_alignment",
                                        quality="quality",
                                        quality_num="quality_num",
                                        v_cigar="v_cigar",
                                        d_cigar="d_cigar",
                                        j_cigar="j_cigar",
                                        np1_length="np1_length",
                                        np2_length="np2_length",
                                        v_sequence_end="v_sequence_end",
                                        d_sequence_end="d_sequence_end",
                                        raw=FALSE) {
   # query sequence
   sequence <- sequence_db[[sequence]]
   quality_phred <- strsplit(sequence_db[[quality]],"")[[1]]
   quality_num_values <- strsplit(sequence_db[[quality_num]],",")[[1]]
   
   v_cigar <- sequence_db[[v_cigar]]
   vd_pseudo_cigar <- NA
   if (!is.na(sequence_db[[np1_length]])) {
      if (sequence_db[[np1_length]]>0) {
         vd_pseudo_cigar <- paste0(sequence_db[[v_sequence_end]],"S",sequence_db[[np1_length]],"X")
      }
   }
   d_cigar <- sequence_db[[d_cigar]]
   dj_pseudo_cigar <- NA
   if (!is.na(sequence_db[[np2_length]])) {
      if (sequence_db[[np2_length]]>0){
         dj_pseudo_cigar <- paste0(sequence_db[[d_sequence_end]],"S",sequence_db[[np2_length]],"X")
      }
   }
   j_cigar <- sequence_db[[j_cigar]]
   cigars <-  c(v_cigar, vd_pseudo_cigar, d_cigar, dj_pseudo_cigar, j_cigar)
   cigars <- cigars[!is.na(cigars)]

   ranges <- bind_rows(lapply(cigars, function(cigar) {
      ops <- GenomicAlignments::explodeCigarOps(cigar)[[1]]
      lengths <- GenomicAlignments::explodeCigarOpLengths(cigar)[[1]]
      keep <- ops %in% c("N", "I") == F
      ops <- ops[keep]
      lengths <- lengths[keep]

      ranges <- data.frame(
         "start"=rep(NA, length(ops)),
         "end"=NA,
         "width"=lengths,
         "operator"=ops,
         stringsAsFactors = F)
      ranges[['start']][1] <- 1

      for (i in 1:nrow(ranges)) {
         if (ranges[['operator']][i] %in% c("S","=","X","D")) {
            ranges[['end']][i] <- ranges[['start']][i]+ranges[['width']][i]-1
         }
         if (i+1<=nrow(ranges)) {
            ranges[['start']][i+1] <- ranges[['end']][i]+1
         }
      }

      ranges <- ranges %>%
         filter(!!rlang::sym("operator")  %in% c("S","D") == FALSE)
      ranges
   }))

   iranges <- IRanges(start=ranges[['start']], end=ranges[['end']], width=ranges[['width']])
   reconstruced_sequence_alignment <- extractAt(BString(sequence), iranges)
   reconstruced_sequence_alignment <- paste0(sapply(reconstruced_sequence_alignment, toString), collapse="")
   positions <- unlist(sapply(1:nrow(ranges), function(i) { ranges$start[i]:ranges$end[i]}))

   quality_df <- data.frame(
      "reconstructed_sequence_alignment"=paste0(sapply(reconstruced_sequence_alignment, toString),collapse=""),
      "sequence_position"=positions,
      "sequence_alignment_position"=NA,
      stringsAsFactors = F
   ) %>%
   mutate (
      !!rlang::sym(quality) := quality_phred[positions],
      !!rlang::sym(quality_num) := as.numeric(quality_num_values)[positions],
      !!rlang::sym(sequence_id) := sequence_db[[sequence_id]])

   # Sanity check. Just to be sure reconstruction is working correctly,
   # and be sure I will later transfer the quality scores correctly
   sequence_alignment <- sequence_db[[sequence_alignment]]
   expected <- gsub("[\\.-]","",sequence_alignment)
   if (reconstruced_sequence_alignment != expected) {
      stop("Reconstructed sequence_alignment from cigar doesn't match db sequence_alignment.")
   }

   # map position numbering: sequence input positions <--> aligned positions
   nt_aln <- strsplit(sequence_alignment,"")[[1]]
   for ( aln_position in 1:length(nt_aln)) {
      if (nt_aln[aln_position] %in% c(".","-") == FALSE ) {
         pos <- sum(nt_aln[1:aln_position] %in% c(".","-") == F)
         quality_df[['sequence_alignment_position']][pos] <- aln_position
         quality_df[['sequence_alignment_nt']][pos] <- nt_aln[aln_position]
      } 
   }
   
   if (raw) {
      quality_df %>%
         select(-!!rlang::sym("reconstructed_sequence_alignment"))
   } else {
      qual_num <- rep(NA, length(nt_aln))
      qual_num[quality_df[['sequence_alignment_position']]] <- quality_df[[quality_num]]
      qual_num <- paste0(qual_num, sep="", collapse=",")
      
      qual_phred <- rep(" ", length(nt_aln))
      qual_phred[quality_df[['sequence_alignment_position']]] <- quality_df[[quality]]
      qual_phred <- paste0(qual_phred, sep="", collapse="")
      
      ret <- data.frame(
                 "quality_alignment_num"=qual_num,
                 "quality_alignment"=qual_phred,
                 stringsAsFactors = F
                 ) %>%
      mutate(!!rlang::sym(sequence_id) := sequence_db[[sequence_id]])
      ret
   }

}


# Retrieve sequencing quality scores from tabular data
# 
# \code{sequenceAlignmentQuality} is used internally by \code{readFastqDb} to 
# process the sequencing quality scores loaded from a \code{fastq} file. 
# 
# Once a repertoire \code{data.frame} has been processed with \link{readFastqDb} and 
# contains the fields \code{quality} and \code{quality_num},
# \code{sequenceAlignmentQuality} can be used to retrieve the quality scores 
# from the already present field \code{quality_num}, without requiring 
# again the \code{fastq} file, and report them as a \code{data.frame} with sequencing 
# qualities per position, not as a string. This is done setting \code{raw=TRUE}. 
# This \code{data.frame} with qualities per position can be used to generate figures, 
# for example.
# 
# @param    data            \code{data.frame} containing sequence data.
# @param    sequence_id     column in \code{data} that contains sequence 
#                           identifiers to be matched to sequence identifiers in 
#                           \code{fastq_file}. 
# @param    sequence        column in \code{data} that contains sequence data.
# @param    sequence_alignment     column in \code{data} that contains 
#                                  IMGT aligned sequence data.      
# @param    quality       column in \code{data} that contains 
#                                  sequencing quality as Phred scores.     
# @param    quality_num   column in \code{data} that contains 
#                                  sequencing quality as a comma separated string.                
# @param    v_cigar         column in \code{data} that contains CIGAR 
#                           strings for the V gene alignments.     
# @param    d_cigar         column in \code{data} that contains CIGAR 
#                           strings for the D gene alignments.   
# @param    j_cigar         column in \code{data} that contains CIGAR 
#                           strings for the J gene alignments.     
# @param    np1_length      column in \code{data} that contains the number
#                           of nucleotides between the V gene and first D gene 
#                           alignments or between the V gene and J gene alignments.
# @param    np2_length      column in \code{data} that contains the number
#                           of nucleotides between either the first D gene and J 
#                           gene alignments or the first D gene and second D gene
#                           alignments.
# @param    v_sequence_end  column in \code{data} that contains the 
#                           end position of the V gene in \code{sequence}.
# @param    d_sequence_end  column in \code{data} that contains the 
#                           end position of the D gene in \code{sequence}.                      
# @param    raw             specify how the sequencing quality should be returned. 
#                           If \code{TRUE}, return a \code{data.frame} with 
#                           quality information per position, where each row is
#                           a position. This \code{data.frame} has columns
#                           "sequence_position", "sequence_alignment_position", 
#                           "quality", "quality_num", "sequence_id"
#                           and "sequence_alignment_nt" If \code{FALSE}, for each sequence,
#                           concatenate the position qualities in a string, and the
#                           quality information to \code{data}
sequenceAlignmentQuality <- function(data, 
                                     sequence_id="sequence_id",
                                     sequence="sequence",
                                     sequence_alignment="sequence_alignment",
                                     quality="quality",
                                     quality_num="quality_num",
                                     v_cigar="v_cigar",
                                     d_cigar="d_cigar",
                                     j_cigar="j_cigar",
                                     np1_length="np1_length",
                                     np2_length="np2_length",
                                     v_sequence_end="v_sequence_end",
                                     d_sequence_end="d_sequence_end",
                                     raw=FALSE) {
   pb <- progressBar(nrow(data))
   qual <- bind_rows(lapply(1:nrow(data),function(i) {
      pb$tick()
      calcSequenceAlignmentQuality(data[i,],
      sequence=sequence,
      sequence_id=sequence_id,
      sequence_alignment=sequence_alignment,
      quality=quality,
      quality_num=quality_num,
      v_cigar=v_cigar,
      d_cigar=d_cigar,
      j_cigar=j_cigar,
      np1_length=np1_length,
      np2_length=np2_length,
      v_sequence_end=v_sequence_end,
      d_sequence_end=d_sequence_end,
      raw=raw)
   }))
   
   if (raw) {
      qual   
   } else {
      data %>%
         dplyr::left_join(qual, by=sequence_id)
   }
}


#' Mask sequence positions with low quality
#' 
#' \code{maskPositionsByQuality} will replace positions that 
#' have a sequencing quality score lower that \code{min_quality} with an
#' \code{"N"} character.
#' 
#' 
#' @param    data          \code{data.frame} containing sequence data.
#' @param    min_quality   minimum quality score. Positions with sequencing quality 
#'                         less than \code{min_qual} will be masked.
#' @param    sequence      column in \code{data} with sequence data to be masked.
#' @param    quality_num   column in \code{data} with quality scores (a
#'                         string of numeric values, comma separated) that can
#'                         be used to mask \code{sequence}. 
#'                         
#' @return   Modified \code{data} data.frame with an additional field containing 
#'           quality masked sequences. The  name of this field is created 
#'           concatenating the \code{sequence} name and \code{"_masked"}.
#'           
#' @seealso \link{readFastqDb} and \link{getPositionQuality}
#' 
#' @examples
#' db <- airr::read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
#' fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
#' db <- readFastqDb(db, fastq_file, quality_offset=-33)
#' maskPositionsByQuality(db, min_quality=90, quality_num="quality_alignment_num")
#' 
#' @export
maskPositionsByQuality <- function(data, min_quality=70,
                                   sequence="sequence_alignment",
                                   quality_num="quality_alignment_num") {
   
   required_cols <- c(sequence,quality_num)
   checkColumns(data, required_cols)
   sequence_masked <- paste0(sequence,"_masked")
   num_masked_seqs <- 0
   data <- bind_rows(lapply(1:nrow(data), function(i) {
      db_row <- data[i,]
      seq_qual <- strsplit(db_row[[quality_num]],",")[[1]]
      low_seq_qual <- which(sapply(seq_qual, function(x) {
         if (x != "NA") {
            as.numeric(x) < min_quality
         } else {
            NA
         }
      }, USE.NAMES = FALSE))
      if (length(low_seq_qual)>0) {
         num_masked_seqs <<- num_masked_seqs + 1
         seq <- strsplit(db_row[[sequence]],"")[[1]]
         seq[low_seq_qual] <- "N"
         seq <- paste0(seq, collapse="")
         db_row[[sequence_masked]] <- seq
      }
      db_row
   }))
   message("Number of masked sequences: ", num_masked_seqs)
   data
} 


#' Get a data.frame with sequencing qualities per position
#' 
#' \code{getPositionQuality} takes a data.frame with sequence quality scores 
#' in the form of a strings of comma separated numeric values, split the quality 
#' scores values by \code{","},  and returns a data.frame with the values
#' for each position.
#' 
#' 
#' @param    data          \code{data.frame} containing sequence data.
#' @param    sequence_id   column in \code{data} with sequence identifiers.
#' @param    sequence      column in \code{data} with sequence data. 
#' @param    quality_num   column in \code{data} with quality scores (as
#'                         strings of numeric values, comma separated) for \code{sequence}.
#'                         
#' @return   \code{data} with one additional field with masked sequences. The 
#'           name of this field is created concatenating \code{sequence} 
#'           and '_masked'.
#'
#' @seealso \link{readFastqDb} and \link{maskPositionsByQuality}
#'           
#' @examples
#' db <- airr::read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
#' fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")
#' db <- readFastqDb(db, fastq_file, quality_offset=-33)
#' head(getPositionQuality(db))
#
#' @export          
getPositionQuality <- function(data, 
                               sequence_id="sequence_id",
                               sequence="sequence_alignment", 
                               quality_num="quality_alignment_num") {
   
   checkColumns(data, c(sequence, quality_num))
   
   bind_rows(lapply(1:nrow(data), function(i) {
      seq_id <- data[[sequence_id]][i]
      seq_len <- nchar(data[[sequence]][i])
      qual_values <- as.numeric(strsplit(data[[quality_num]][i],",")[[1]])
      nt <- strsplit(data[[sequence]],"")[[1]]
      if (seq_len != length(qual_values)) {
         stop("Different length, for sequence: ", seq_id,". seq: ", sequence, ". qual: ", quality_num)
      }
      data.frame(
         "position"=1:seq_len, stringsAsFactors = F) %>%
         mutate(!!rlang::sym(quality_num) := qual_values,
                !!rlang::sym(sequence_id) := seq_id,
                nt = nt)
   }))
   
}
