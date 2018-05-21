#' @title calculate the distance to the nearest 5' or 3' subject GRanges
#'
#' @description \code{distance_map} calculate the distance to the nearest 5' or 3' end.
#' @param query_gr The query GRanges.
#' @param subject_gr The subject GRanges.
#' @param end Can be "five" or "three", indicating the distance to the five prime end or three prime end.
#' @param maximum_distance The maximum distance returned, the un-defined distances are set to this value also, default is 2000.
#' @param standardize wheather standardize the output distance, i.e. minus the means and dividing the standard deviations.
#' @return a vector of mapped distances with the length equal to \code{query_gr}.
#'
distance_map <- function(query_gr,
                         subject_gr,
                         end = c("five","three"),
                         maximum_distance = 2000,
                         standardize = T){
  end <- match.arg(end)

  if(end == "five"){
  match_indx <- precede(query_gr,subject_gr,ignore.strand = F)
  }

  if(end == "three"){
    match_indx <- follow(query_gr,subject_gr,ignore.strand = F)
  }

  na_indx <- is.na(match_indx)

  dist_all <- vector("integer",length(query_gr))

  dist_all[!na_indx] <- distance(query_gr[!na_indx],
                            subject_gr[match_indx[!na_indx]])

  dist_all[is.na(dist_all)] <- maximum_distance

  dist_all <- pmin(dist_all,maximum_distance)

  if(standardize) {
    return((dist_all-mean(dist_all))/sd(dist_all))
  } else {
    return(dist_all)
  }
}
