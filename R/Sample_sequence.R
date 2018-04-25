#' @title Sample nucleotides or motifs on a subset of genome defined by \code{GRanges}.
#'
#' @description \code{sample_sequence} is used to extract the regions of a query sequence mapped to user defined sub-regions of Genome.
#'
#' @param Sequence a character string indicating the query sequence or motifs; Rules in \code{\link{IUPAC_CODE_MAP}} is supported when \code{Fixed} = FALSE.
#'
#' @param Subset_GR the \code{\link{GRanges}} object used to define the subset of the genome.
#'
#' @param BSgnm the \code{\link{BSgenome}} object containing the sequence of the genome.
#'
#' @param Fixed FALSE to support the vague mapping rule of the chagracter string, default is FALSE.
#'
#' @param N number of ranges sampled, by default, it returns all the ranges without sampling.
#'
#' @param Replace whether sample with replacement or not, default is FALSE.
#'
#' @return a \code{GRanges} object contains the (sampled) mapped regions of the query sequence on the defined subset of the genome.
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' exons_hg19 <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' Motif_exons_hg19 <- Sample_sequence_on_TXDB("RRACH",exons_hg19,Hsapiens,N = 100000)
#' #The command above sample 100000 regions of RRACH on exonic regions of hg19.
#'
#' @import BSgenome
#' @import GenomicFeatures
#'
#' @export
sample_sequence <- function(Sequence, Subset_GR, BSgnm,Fixed = F, N = NULL, Replace = F){
  regions_reduced <-  reduce(Subset_GR)
  regions_DNASS <- DNAStringSet( Views(BSgnm,regions_reduced) )
  regions_GRL <- split(regions_reduced, paste0("EX_",1:length(regions_reduced)))
  regions_GRL <- regions_GRL[paste0("EX_",1:length(regions_reduced))]
  MIndx <- vmatchPattern(Sequence,regions_DNASS,fixed = Fixed)
  MIndx_gr <- GRanges(seqnames = rep(names(regions_GRL),elementNROWS(MIndx)),ranges = unlist(MIndx))
  sequence_on_regions <- mapFromTranscripts(MIndx_gr,regions_GRL)
  mcols(sequence_on_regions) = NULL
  if(is.null(N)){}else{
  if(Replace == F) N = min(N,length(sequence_on_regions))
  Indx <- sample.int(length(sequence_on_regions),N,replace = Replace)
  sequence_on_regions <- sequence_on_regions[Indx]
  }
  return(sequence_on_regions)
}
