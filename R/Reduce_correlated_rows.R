#' @title Subset the rows/genomic features that are close on genomic coordinate to reduce the dependencies between rows.
#'
#' @description \code{reduce_correlated_rows} is a function to define a set of metrics to subset the rows that are highly correlated with its nearby features on genomic scale.
#' @param SE A \code{SummarizedExperiment} with rowRanges being the GRanges of row feaures, and an assay matrix with at least one collumn.
#' @param cor_method The method to define correlations between rows of the assay matrix, can be one in "spearman" and "pearson".
#' @param cor_cut_off The correlation cut off threshold used to group 2 nearby rows, default is 0.8.
#' @param bin_width The bin width to define the closed/neighbooring row features, default is 101.
#' @param reduction_method The decision criteria for the grouped correlated rows:
#'
#' "maxSum" : keep the closed and correlated row features with the highest row sums.
#'
#' "maxMad" : keep the closed and correlated row features with the highest row Median Absolute Deviation.
#'
#' "maxInfo": keep the closed and correlated row features with highest row total information defined by \code{information_matrix}.
#'
#' "random": keep one of the correlated row features randomly.
#'
#' @param information_matrix The information matrix used when argument \code{reduction_method} = "maxInfo".
#'
#' @details The correlation between closed row features are caculated, the rows have mutually correlated neighboors are grouped together, only one of the member in the group will be kept using one of the metric defined obove.
#'
#' @return The output is a \code{SummarizedExperiment} object with subsetted rows compared with the input.
#'
#' @examples
#' reduce_correlated_rows(SE_CQN,"spearman",".7",101,"maxInfo",assays(SE_CQN)$IP + assays(SE_CQN)$input)
#'
#' @import SummarizedExperiment
#' @export

reduce_correlated_rows <- function(SE,
                                   cor_method = "spearman",
                                   cor_cut_off = 0.8,
                                   bin_width = 101,
                                   reduction_method = "maxSum",
                                   information_matrix = NULL) {

stopifnot(cor_method %in% c("spearman","pearson"))
stopifnot(reduction_method %in% c("maxSum","maxMad","maxInfo","random"))
stopifnot(is.numeric(cor_cut_off))

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
    spearman_M <- cor(na.omit(t(df[,idx])),method = cor_method)
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
  if(reduction_method == "maxMad"){
    return(x[ which.max( rowMads( as.matrix(x[,-1*ncol(x)])  ,na.rm = T)),"idx"])
  }
  if(reduction_method == "maxInfo"){
    return(x[ which.max( rowSums( cbind( information_matrix[x[,ncol(x)],])  ,na.rm = T)),"idx"])
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

message(paste0( nrow(SE) - length(Keep_row_indx) , " rows are dropped.") )

return(SE[Keep_row_indx,])
}
