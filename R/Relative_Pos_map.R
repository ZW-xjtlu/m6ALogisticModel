#' @title Calculate the relative position of a GRanges object on transcript.
#'
#' @description \code{relative_pos_map} is used to calculate the relative position of a GRanges object on transcript regions defined by a GRangesList object.
#' The query GRanges is assumed to have a single based resolution (or width = 1), otherwise only the start location is considered.
#' @param query_gr The GRanges to query.
#' @param feature_grl The GRangesList to map on.
#' @param no_map_val The value give to the not mapped query, by default is NA.
#' @return a numeric vector with value between [0,1].
#'
relative_pos_map <- function(query_gr, feature_grl, no_map_val = NA, standardize = TRUE) {

  return_vec <- rep(no_map_val,length(query_gr))

  map2tx <- mapToTranscripts(query_gr,feature_grl)

  relat_pos <- start(map2tx)/ sum(width(feature_grl))[map2tx$transcriptsHits]

  relpos_idx <- tapply( relat_pos, map2tx$xHits, mean)

  return_vec[as.numeric(names(relpos_idx))] = relpos_idx

  if(standardize) {
    return((return_vec - mean(return_vec,na.rm = T))/sd(return_vec,na.rm = T))
  } else {
    return(return_vec)
  }
}
