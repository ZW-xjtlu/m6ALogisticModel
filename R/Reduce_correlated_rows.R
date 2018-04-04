#' @title Subset the rows/genomic features that are close on genomic coordinate to reduce the row dependencies.
#'
#' @description \code{Reduce_correlated_rows} is a function to define a set of metrics to subset the rows that are highly correlated with its nearby features on genomic scale.
#' @param SE A \code{SummarizedExperiment} with rowRanges being the GRanges of row feaures, and an assay matrix with at least one collumn.
#' @param cor_method The method to define correlations between rows of the assay matrix, can be one in "spearman", "pearson".
#' @param cor_cut_off The correlation cut_off to reduce 2 nearby rows, default is 0.8.
#' @param bin_width The bin width to define the closed/neighbooring row features, default is 101.
#' @param reduction_method The decision criteria on the highly correlated rows:
#'
#' "maxSum" : keep the closed and correlated row features with highest row sums.
#'
#' "maxVar" : keep the closed and correlated row features with highest row variance.
#'
#' "maxInfo": keep the closed and correlated row features with highest row total information defined by \code{information_matrix}.
#'
#' "Random": keep one of the closed and correlated row features randomly.
#'
#' @param information_matrix Wheather to return the clustering index, default is FASLE.
#'
#' @details The correlation values between closed rows are caculated, if the rows are mutually correlated neighboors, only one of them will be kept by the metric defined obove.
#'
#' @return The output is a \code{SummarizedExperiment} with subsetted rows compared with the input.
#'
#' @examples
#' Reduce_correlated_rows(SE_CQN,"spearman",".7",101,"maxInfo",assays(SE_CQN)$IP + assays(SE_CQN)$input)
#'
#' @import SummarizedExperiment

Reduce_correlated_rows <- function(SE,
                                   cor_method = "spearman",
                                   cor_cut_off = 0.8,
                                   bin_width = 101,
                                   reduction_method = "maxSum",
                                   information_matrix = NULL) {

stopifnot(cor_method %in% c("spearman","pearson"))
stopifnot(reduction_method %in% c("maxSum","maxVar","maxInfo","Random"))

bin <- resize( rowRanges(SE) , bin_width )
mcols(bin) = NULL
reducedbin <- reduce(bin)
overlapsum <- findOverlaps(bin, reducedbin)
values(bin) <- DataFrame(group = subjectHits(overlapsum), M = assay(SE), mcols(bin))
bingroup <- as.data.frame(mcols(bin))

bingroup$idx = 1:dim(bingroup)[1]

bin_lst <- split(bingroup[,-1],bingroup[,1])

func_cor <- function(df){
  if(dim(df)[1] == 1) return(0) else {
    idx <- which(!colnames(df) %in% c("modstart","modName","geneName","geneType","gene_id","idx"))
    spearman_M <- cor(t(df[,idx]),method = cor_method)
    idx_M <- spearman_M > cor_cut_off
    idx_c <- combn(1:dim(idx_M)[1],2)
    idx_v <- idx_M[t(idx_c)]

    fac_v <- 1:dim(idx_M)[1]
    x = NULL
    y = NULL
    for(i in 1:length(idx_v)) {
      x <- idx_c[1,i]
      y <- idx_c[2,i]
      if(is.na(idx_v[i])){} else {
        if(idx_v[i]) {
          if (fac_v[x] != fac_v[y]) {
            replacement = fac_v[x]
            replaced = fac_v[y]
            fac_v[which(fac_v == replaced)] = replacement
          }
        }
      }
    }
    return(fac_v)
  }
}

func_decision <- function(df){
  idx <- func_cor(df)
  df_lst <- split(df,idx)

 Indx_SE <- sapply(df_lst,function(x) {

  if(reduction_method  == "maxSum"){
    return(x[ which.max( rowSums( x[,-1*ncol(x)]  ,na.rm = T)),"idx"])
  }
  if(reduction_method == "maxVar"){
    return(x[ which.max( rowVars( x[,-1*ncol(x)]  ,na.rm = T)),"idx"])
  }
  if(reduction_method == "maxInfo"){
    return(x[ which.max( rowSums( information_matrix[x[,ncol(x)],]  ,na.rm = T)),"idx"])
  }
  if(reduction_method == "Random"){
    return(sample(x[,ncol(x)],1))
  }

}
)

 names(Indx_SE) = NULL
 return(Indx_SE)
}


Keep_row_indx <- unlist(sapply(bin_lst,func_decision))

return(SE[Keep_row_indx,])
}
