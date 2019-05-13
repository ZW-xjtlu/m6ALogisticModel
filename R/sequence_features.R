#' @title Generate the sequence feature matrix from GRanges.
#'
#' @description \code{sequence_features} is used to extract the sequence features given GenomicRanges object.
#'
#' @param query_gr a GRanges object to generate the underlying sequence features, all the widths of the ranges must be equal.
#' @param bsgnm a BSgenome object which encode the genome information.
#' @return a data.frame contains the sequence features in its collumns.
#' @import BSgenome
#' @export
#'
#'
sequence_features <- function(query_gr, bsgnm) {

  stopifnot(all(width(query_gr) == width(query_gr)[1]))

  N = width(query_gr)[1]

  bsgnm_view <- Views(bsgnm,query_gr)

  sequences <- as.character( DNAStringSet(bsgnm_view) )

  rm(bsgnm_view)

  sequences_M <- matrix( unlist( strsplit(sequences,"") ) , ncol =  N, byrow = T)

  rm(sequences)
  #vectorized solution
  purine_M <- sequences_M == "A" | sequences_M == "G"
  colnames(purine_M) <-  paste0("purine_",seq_len(N))
  amino_M <- sequences_M == "A" | sequences_M == "C"
  colnames(amino_M) <-  paste0("amino_",seq_len(N))
  weakHyb_M <- sequences_M == "A" | sequences_M == "T"
  colnames(weakHyb_M) <- paste0("weakHyb_",seq_len(N))

  cumFreq_A <- rowCumsums(matrix(as.numeric( sequences_M == "A" ), ncol = N, byrow = F))
  cumFreq_T <- rowCumsums(matrix(as.numeric( sequences_M == "T" ), ncol = N, byrow = F))
  cumFreq_C <- rowCumsums(matrix(as.numeric( sequences_M == "C" ), ncol = N, byrow = F))
  cumFreq_G <- rowCumsums(matrix(as.numeric( sequences_M == "G" ), ncol = N, byrow = F))

  cumFreq_combined <- matrix(0,ncol = N, nrow = length(query_gr))
  cumFreq_combined[sequences_M == "A"] <- cumFreq_A[sequences_M == "A"]
  cumFreq_combined[sequences_M == "T"] <- cumFreq_T[sequences_M == "T"]
  cumFreq_combined[sequences_M == "C"] <- cumFreq_C[sequences_M == "C"]
  cumFreq_combined[sequences_M == "G"] <- cumFreq_G[sequences_M == "G"]

  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined ) <-  paste0("cumFreq_",seq_len(N))

  return(as.data.frame(cbind(purine_M,amino_M,weakHyb_M,cumFreq_combined)))
}
