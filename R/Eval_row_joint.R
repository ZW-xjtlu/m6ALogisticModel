#' @title Evaluate the joint distribution between rows of a genomic assay.
#'
#' @description \code{Plot_column_joint} is a function used to evaluate the clustering quality between m6A sites.
#' @param SE A \code{SummarizedExperiment} with features annotated by \code{\link{predictors.annot}}, the colnames of the SummarizedExperiment should be sample names.
#' @param HDER The subtitle and the file name of the plot.
#' @param K The number of centers used in K medoids clustering.
#' @param ROW_STAND Wheather standardize rows before clustering, default is TRUE.
#' @param RETURN_INDX Wheather to return the clustering index, default is FASLE.
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
#' About deviance in GLM (from stackexchange).
#'
#'Let LL = loglikelihood
#'Here is a quick summary of what you see from the summary(glm.fit) output,
#'Null Deviance = 2(LL(Saturated Model) - LL(Null Model)) on df = df_Sat - df_Null
#'Residual Deviance = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Proposed
#'The Saturated Model is a model that assumes each data point has its own parameters (which means you have n parameters to estimate.)
#'The Null Model assumes the exact "opposite", in that is assumes one parameter for all of the data points, which means you only estimate 1 parameter.
#'The Proposed Model assumes you can explain your data points with p parameters + an intercept term, so you have p+1 parameters.
#'If your Null Deviance is really small, it means that the Null Model explains the data pretty well. Likewise with your Residual Deviance.
#'What does really small mean? If your model is "good" then your Deviance is approx Chi^2 with (df_sat - df_model) degrees of freedom.
#'If you want to compare you Null model with your Proposed model, then you can look at
#'(Null Deviance - Residual Deviance) approx Chi^2 with df Proposed - df Null = (n-(p+1))-(n-1)=p
#'Are the results you gave directly from R? They seem a little bit odd, because generally you should see that the degrees of freedom reported on the Null are always higher than the degrees of freedom reported on the Residual. That is because again, Null Deviance df = Saturated df - Null df = n-1 Residual Deviance df = Saturated df - Proposed df = n-(p+1)
#'
#' @return Plots saved under a directory named by \code{HDER}.
#'
#' @examples
#'library(SummarizedExperiment)
#'library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'library(BSgenome.Hsapiens.UCSC.hg19)
#'library(fitCons.UCSC.hg19)
#'library(phastCons100way.UCSC.hg19)
#'
#'Feature_List_expanded_hg19 = list(
#'  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
#'  YTHDC1_TREW = YTHDC1_TREW_gr,
#'  YTHDF1_TREW = YTHDF1_TREW_gr,
#'  YTHDF2_TREW = YTHDF2_TREW_gr,
#'  miR_targeted_genes = miR_targeted_genes_grl,
#'  TargetScan = TargetScan_hg19_gr,
#'  Verified_miRtargets = verified_targets_gr,
#'  METTL3_TREW = METTL3_TREW,
#'  METTL14_TREW = METTL14_TREW,
#'  WTAP_TREW = WTAP_TREW,
#'  METTL16_CLIP = METTL16_CLIP,
#'  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
#'  FTO_CLIP = FTO_CLIP,
#'  FTO_eCLIP = FTO_eCLIP
#')
#'
#'SE_features_added <- predictors.annot(se = SE_CQN,
#'                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                      bsgnm = Hsapiens,
#'                                      fc = fitCons.UCSC.hg19,
#'                                      pc = phastCons100way.UCSC.hg19,
#'                                      struct_hybridize = Struc_hg19,
#'                                      feature_lst = Feature_List_expanded_hg19,
#'                                      HK_genes_list = HK_hg19_eids)
#'
#'Eval_row_joint(SE_features_added, "Row_joint_CQN")
#'
#' @seealso \code{\link{predictors.annot}}
#'
#'
#' @import ggplot2
#' @import cluster
#' @import nnet
#' @import SummarizedExperiment
#' @importFrom reshape2 melt
#' @export
Eval_row_joint <- function(SE, HDER = "Row_joint", K = 3, ROW_STAND = T, RETURN_INDX = F) {
stopifnot(class(SE)=="RangedSummarizedExperiment")
stopifnot(!is.null(mcols(SE)))
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

#Fit multinomial GLM to analysis the clustering label
GLM_df = mcols(SE)
GLM_df$Y = factor( paste0("Cluster ",row_cluster$clustering) )
GLM_df$Y = relevel(GLM_df$Y, "Cluster 1")
Null_model <-  multinom(Y ~ 1, data = GLM_df)
Proposed_model <-  multinom(Y ~ ., data = GLM_df)
stat_test <- summary(Proposed_model)
z <- stat_test$coefficients/stat_test$standard.errors

#Prepare a data frame for the ggplot.
Plot_df_Z <- melt(z)
Plot_df_Z$Stat <- "Wald_Z"
Plot_df_Estimate <- melt(stat_test$coefficients)
Plot_df_Estimate$Stat <- "log(OddsRatio)"

if(ncol(Plot_df_Z) == 2) {
  Plot_df_Z = cbind("Cluster 2", rownames(Plot_df_Z), Plot_df_Z )
  colnames(Plot_df_Z) = c("Clusters","Covariates","values","Statistics")
  Plot_df_Estimate = cbind("Cluster 2", rownames(Plot_df_Estimate), Plot_df_Estimate)
  colnames(Plot_df_Estimate) = c("Clusters","Covariates","values","Statistics")
}

Plot_df <- rbind(Plot_df_Z,Plot_df_Estimate)

colnames(Plot_df) <- c("Clusters","Covariates","values","Statistics")


Plot_df$Covariates = gsub("TRUE","",Plot_df$Covariates)

Indx <- names( sort( tapply(Plot_df$values[Plot_df$Statistics == "Wald_Z"] ,
                            Plot_df$Covariates[Plot_df$Statistics == "Wald_Z"],
                            function(x) sum(abs(x))) , decreasing = F))

Plot_df$Covariates = factor(Plot_df$Covariates, levels = Indx)

Plot_df$values_abs = abs(Plot_df$values)
Plot_df$Sign = "positive"
Plot_df$Sign[Plot_df$values < 0] = "negative"
Plot_df$Group = paste0(Plot_df$Statistics,":",Plot_df$Clusters)
Indx_Group = unique(Plot_df$Group)

n = length(Indx_Group)
indx <- 1:n
for(i in 1:n){
  if(i%%2 != 0)  {
    indx[i] = (i%%2 + i)/2
  }else{
    indx[i] =  n/2  + i/2
  }
}

Plot_df$Group = factor( Plot_df$Group,
                                levels = Indx_Group[indx] )

Z = Plot_df$values[Plot_df$Statistics == "Wald_Z"]

p <- (1 - pnorm(abs(z), 0, 1))*2

Sig_indx <- p.adjust(p,method = "fdr") < .05

Plot_df$FDR_sig = "Insig"
Plot_df$FDR_sig[Plot_df$Statistics == "Wald_Z"][Sig_indx] = "< .05"

p2 <- ggplot(Plot_df) +
  geom_bar(stat = "identity",aes(x = Covariates, y = values_abs, fill = Sign, colour = FDR_sig), width = .75,size = .35, linetype = 1) +
  facet_wrap(~Group,nrow = 1,scales = "free_x") +
  theme_classic() + coord_flip() + theme(axis.text.y = element_text(size = 7),
                                         plot.margin = margin(1,1.5,1,1,"cm")) +
  scale_fill_brewer(direction = -1) + scale_colour_manual(values = c("red", 0)) +
  labs(title = paste0("Multinomial logistic model for clustering result: ",HDER),
       subtitle = paste0("Proposed model Chisq statistics: ",
                        round(  Null_model$deviance-Proposed_model$deviance, 3),
                         " on ",length(Indx) * (K-1) - 1 ," df"),
       y = "abs(values)")

fig_width_p2 = 5 + 2.5*(K-1) + .01 * max(nchar(Indx))

fig_height_p2 = 4 + .08 * length(Indx)

ggsave(paste0(HDER,"_GLMestimates.pdf"),p2,width = fig_width_p2,height = fig_height_p2)

#Calculate the total statistical significance of the proposed model.
Stat_df <- data.frame(
NULL_Deviance = Null_model$deviance,
Residual_Deviance = Proposed_model$deviance
)

Stat_df$Reduced_prop = 1 - Stat_df$Residual_Deviance/Stat_df$NULL_Deviance
Stat_df$Cost_df = length(Indx) * (K-1) - 1
Stat_df$Chisq_stat = Stat_df$NULL_Deviance - Stat_df$Residual_Deviance
write.table(t(Stat_df),paste0("Model_report_",HDER,".txt"),col.names = F)

if(RETURN_INDX){
  return(row_cluster$clustering)
}
}

