#' @title Dealing with ambiguities of the GRanges in query that may mapped to multiple features with different properties.
#'
#' @description \code{properties_map} mitigates the ambiguities by averaging the properties of multiple mapped query.
#' @param query_gr The GRanges to query.
#' @param feature_grl The GRangesList to map on, e.x. the Genes derived by exonsby(txdb ,by = "genes").
#' @param feature_properties The vector of properties of the subject (e.x. the length of the Gene).
#' @param no_map_val The value give to the not mapped query, by default is NA.
#' @param normalize wheather standardize the output mapped properties, i.e. minus the means and dividing the SDs.
#' @return a vector of mapped properties with the length equal to \code{query_gr}.
#'

properties_map <- function(query_gr, feature_gr, feature_properties, no_map_val = NA, normalize = F) {

  fol <- findOverlaps(query_gr,feature_gr)

  return_vec <- rep(no_map_val,length(query_gr))

  features_mapped <- feature_properties[subjectHits(fol)]

  if(is.logical(feature_properties)) {
    Weighted_vec  <- tapply(features_mapped, queryHits(fol), any)
  }else{
    Weighted_vec <- tapply(features_mapped, queryHits(fol), mean)
  }

  return_vec[ as.numeric( names(Weighted_vec) )] <- Weighted_vec

  if(normalize) {
    known_indx <- !is.na(return_vec)
    return_vec[known_indx] = (return_vec[known_indx] - mean(return_vec[known_indx]))/sd(return_vec[known_indx])
  }

  return(return_vec)
}
