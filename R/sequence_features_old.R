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
sequence_features_old <- function(query_gr, bsgnm) {

stopifnot(all(width(query_gr) == width(query_gr)[1]))

bsgnm_view <- Views(bsgnm,query_gr)

sequences <- as.character( DNAStringSet(bsgnm_view) )

sequences_lst <- strsplit(sequences,"")

chemical_properties <- function(NT = c("A","T","C","G","N")){
  nt <- match.arg(NT)
  feature_chem = c( nt %in% c("A","G"),
                    nt %in% c("A","C"),
                    nt %in% c("A","T"))
  return(feature_chem)
}

nt_frequency <- function(NT_vec) {
  D <- vector("numeric",length(NT_vec))
  for(i in seq_along(NT_vec)) {
    D[i] <- mean(NT_vec[1:i] == NT_vec[i])
  }
  return(D)
}

features_lst <- lapply(sequences_lst, function(x){
  M_features <- rbind( sapply(x,chemical_properties), nt_frequency(x) )
  return(as.vector(M_features))
})

nt_num <- width(query_gr)[1]

feature_M <- matrix(unlist(features_lst),ncol = nt_num*4,byrow = T)

feature_M <- as.data.frame(feature_M)

colnames( feature_M ) <- paste0(rep(c("purine","amino","weakHyb","cumFreq")),"_",rep( 1:nt_num , each = 4))

return(as.data.frame(feature_M))
}
