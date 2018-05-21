#' @title Infer looped RNA 2ndary structures given provided hybridized 2ndary structures.
#'
#' @description \code{infer_loop} is used to extract the roughly accurate looped regions beween closed hybridized pairs given regions of hybridization.
#'
#' @param struc_sterm a \code{GRanges} or \code{GRangesList} object that indicates the hybridized structure region on each transcript.
#'
#' @param length_cut_off the allowed maximum pairing distance between 2 hybridized regions, the default setting is 30.
#' Due to the crude algorithm defined by this function, The higher this value is, the more false positive loops are presented in the result.
#'
#' @return a \code{GRanges} object with regions corresponding to the inferred looped region given the provided hybridized region.
#'
infer_loop <- function(struc_sterm, length_cut_off = 30){
  all_range <- reduce(unlist(range(struc_sterm)))
  struc_sterm_r <- reduce(unlist( struc_sterm ))
  dsj_strucs <- disjoin(c(all_range, struc_sterm_r))
  fol <- findOverlaps(dsj_strucs, struc_sterm_r, type = "equal")
  strucs_loop <- dsj_strucs[-1*queryHits(fol)]
  strucs_loop <- strucs_loop[width(strucs_loop) <= 30]
  return(strucs_loop)
}
