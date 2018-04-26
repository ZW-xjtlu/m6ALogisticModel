#' @title Plot the joint distribution between rows of a genomic assay.
#'
#' @description \code{plot_column_joint} is a function used to evaluate the clustering quality between m6A sites.
#' @param SE A \code{SummarizedExperiment} with features annotated by \code{\link{predictors.annot}}, the colnames of the SummarizedExperiment should be sample names.
#' @param HDER The subtitle and the file name of the plot.
#' @param K Number of centers used in K medoids clustering.
#' @param ROW_STAND Whether to standardize rows before clustering, default is TRUE.
#' @param RETURN_INDX Whether to return the clustering index, default is TRUE.
#'
#' @details By default, a K medoids clustering will be applied between rescaled row entires with metric of euclidean; then, a simplified heat map will be plotted.
#' Finally, a report by multinomial GLM is conducted using the clustering label and features.
#'
#' The suggested quick check row number for SE should be less than 10000, otherwise this function would be time consuming.
#'
#' About silhouette plot (from wikipedia).
#'
#' The average silhouette of the data is another useful criterion for assessing the natural number of clusters.
#' The silhouette of a data instance is a measure of how closely it is matched to data within its cluster and how loosely it is matched to data of the neighbouring cluster,
#' i.e. the cluster whose average distance from the datum is lowest.
#' A silhouette close to 1 implies the datum is in an appropriate cluster, while a silhouette close to −1 implies the datum is in the wrong cluster.
#' Optimization techniques such as genetic algorithms are useful in determining the number of clusters that gives rise to the largest silhouette.
#' It is also possible to re-scale the data in such a way that the silhouette is more likely to be maximised at the correct number of clusters
#'
#' @return A clustering index, also a plot will be saved under a file named by \code{HDER}.
#'
#' @examples
#'
#'eval_row_joint(SE_features_added, "Row_joint_CQN")
#'
#' @seealso \code{\link{predictors.annot}},  \code{\link{glm_multinomial}}, \code{\link{go_multinomial}}
#'
#'
#' @import ggplot2
#' @import cluster
#' @export
plot_row_joint <- function(SE, HDER = "Row_joint", K = 3, ROW_STAND = T, RETURN_INDX = F, PROVIDE_INDX = NULL) {
stopifnot(class(SE)=="RangedSummarizedExperiment")
stopifnot(!is.null(mcols(SE)))
stopifnot(is.null(PROVIDE_INDX) | (length(PROVIDE_INDX) == nrow(SE)))
#Cluster: K-means

if(ROW_STAND){
assay_M <- t(scale(t(assay(SE))))
}else{
assay_M <- assay(SE)
}

row_cluster <- pam(x = dist(assay_M), K)


Mean_row_clusters <- apply(assay_M,2, function(x)tapply(x, row_cluster$clustering, mean, na.rm = T))

mat.melted <- melt(Mean_row_clusters)

colnames(mat.melted) = c("Clusters","Samples","mean_value")
#Plot simplified heat map with ggplot2
#Extract the table of the cluster statistics
Df_x = data.frame(Clusters =  1:K ,
                   Obs_number = c(table(row_cluster$clustering)),
                    MeanSh_widths = round( c( row_cluster$silinfo$clus.avg.widths),3),
                  X = min( 5, ncol(SE)/2))

Df_x$Label = paste0("Mean sil width: ",Df_x$MeanSh_widths,"; Obs #:",Df_x$Obs_number)

p1 <- ggplot(mat.melted) + geom_tile(aes(y = Clusters, x = Samples, fill = mean_value)) +
  scale_fill_gradient(low = "Blue", high = "Yellow") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 310, vjust =1, hjust = 0),
        plot.margin = margin(1,1,1,1,"cm")) +
  labs(title = paste0("K Medoids Clustering With K = ",K,"; ",HDER),
       subtitle = paste0("total mean silinfo width: ",round(row_cluster$silinfo$avg.width,3))) + geom_label(data = Df_x,mapping = aes(x = X,y = Clusters - .25,label = Label),size = 3,label.size = .05)

fig_height_p1 = 3.4 + .3*K + .05 * max(nchar(as.character(colnames(SE))))
fig_width_p1 = 4 + .2 * ncol(SE)
ggsave(paste0(HDER,"_kmedoids.pdf"),p1,width = fig_width_p1,height = fig_height_p1)

if(RETURN_INDX) return(row_cluster$clustering)
}

