#' @title Find the longest transcript given an GRangesList object obtained by exonsByTranscripts.
#'
#' @description \code{find_longest_transcript} is used to find the longest transcript per gene.
#'
#' @param exbtx a GRangesList object obtaied from the exonsByTranscripts function defined in package GenomicFeatures.
#' @param txdb the corresponding txdb object.
#' @return a data.frame indicating a given TXID is longest for the gene or not.
#' @importFrom AnnotationDbi select
#' @import GenomicFeatures
#' @export
#'
find_longest_transcript <- function(exbtx,txdb){

  tx_length <- sum(width(exbtx))

  Map_tx_gene <- suppressMessages( select(txdb,names(exbtx),columns = "GENEID",keytype = "TXID") )

  tx_length <- tx_length[ Map_tx_gene$TXID ]

  nogene_indx <- is.na(Map_tx_gene$GENEID)

  tx_length <- tx_length[!nogene_indx]

  Map_tx_gene <- Map_tx_gene[!nogene_indx,]

  longest_indx_lst <- tapply(tx_length, Map_tx_gene$GENEID, function(x) x == max(x))

  #test:# identical( gsub("\\..*","", names(unlist( longest_indx_lst ))), Map_tx_gene$GENEID)

  Map_tx_gene <- Map_tx_gene[order(Map_tx_gene$GENEID),]

  Map_tx_gene$longest <- unlist( longest_indx_lst )

  return(Map_tx_gene)

}
