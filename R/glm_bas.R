#' @title Conduct bayesian model selection, inference, and visualization for logistic regression modeling.
#'
#' @description \code{glm_bas} is used to perform logistic regression analysis based on previously generated features, it reports logistic regression statistics after bayesian model selection,
#' and it can plot diagrams of logit estimates, wald test / likelihood ratio test statistics, and the goodness of fit (deviance) across various samples.
#'
#' @details \code{glm_bas} build logistic regression model on the features and targets defined by the annotated \code{\link{SummarizedExperiment}} object returned by \code{\link{predictors.annot}}. The model selection is conducted by bayesian model selection method defined by the package \code{\link{BAS}}.
#'
#' @param se A \code{\link{SummarizedExperiment}} object containing the matrix of response variables, each row should represent one modification site, and each collumn should represent a sample or condition. The entries of \code{assay} are integer values of 1 or 0, with 1 indicating the positive class and 0 indicating the negative class used in logistic regression, the uncertain values should be set as \code{NA}.
#'
#' The sample names should be defined by the \code{colnames} of the \code{SummarizedExperiment}, alternatively, it can be defined by the first collumn of the \code{colData}, or by the user defined collumns in \code{colData} by the argument \code{Sample_Names_coldata}.
#'
#' The features or design matrix should be included in the \code{mcols} of the \code{SummarizedExperiment} object.
#'
#' @param beta_prior prior distribution for regression coefficients, default setting is "robust"; see \code{\link{bas.glm}}
#' @param model_prior Family of the prior distribution on the models, default setting is beta.binomial(1,1); see \code{\link{bas.glm}}
#' @param MCMC_iterations Number of models to sample, default setting is 100000,the maximum number is 1^e8.
#' @param decision_method Decision method used in the bayesian model selection, default setting is 'BPM'; see \code{\link{predict.bas}}
#' @param top The number of top models used for BMA related decision making, default setting is NULL (BMA over all models); see \code{\link{predict.bas}}
#' @param save_dir The path to save the statistics and the diagrams of the logistic models, the default path is named "LogisticModel".
#' @param sample_names_coldata Provided column names in \code{colData} for the sample labels.
#' @param group_list Optional, a \code{list} indicating the grouping of features in the output diagrams, by default it uses \code{\link{group_list_default}}.
#'
#' @return
#' Folders under the directory specified by \code{save_dir will} be created, the reports and the diagrams will be saved within them.
#'
#' @seealso Use \code{\link{predictors.annot}} to annotate features.
#'
#' @examples
#' library(SummarizedExperiment)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(fitCons.UCSC.hg19)
#' library(phastCons100way.UCSC.hg19)
#'
#' Feature_List_hg19 = list(
#' HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
#' YTHDC1_TREW = YTHDC1_TREW_gr,
#' YTHDF1_TREW = YTHDF1_TREW_gr,
#' YTHDF2_TREW = YTHDF2_TREW_gr,
#' miR_targeted_genes = miR_targeted_genes_grl,
#' #miRanda = miRanda_hg19_gr,
#' TargetScan = TargetScan_hg19_gr,
#' Verified_miRtargets = verified_targets_gr
#' )
#'
#'
#' SE_features_added <- predictors.annot(se = SE_example,
#'                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                      bsgnm = Hsapiens,
#'                                      fc = fitCons.UCSC.hg19,
#'                                      pc = phastCons100way.UCSC.hg19,
#'                                      struct_hybridize = Struc_hg19,
#'                                      feature_lst = Feature_List_hg19,
#'                                      HK_genes_list = HK_hg19_eids)
#'
#'
#' glm_bas(
#' SE_features_added,
#' MCMC_iterations = 10000,
#' decision_method = "BPM",
#' save_dir = "LogisticModel_x",
#' sample_names_coldata = "ID"
#' )
#'
#' #To Do:
#' #1. make a option called "no model selection", to quikly infer the effect sizes.
#' #1.5 make the response variable could be TRUE or FALSE for logistic model.
#' #2. Isolate the meta plot functions from the original ones.
#' #3. Add choice of Poisson GLM and Gaussian OLM.
#'
#' @import GenomicRanges
#' @import ggplot2
#' @import SummarizedExperiment
#'
#' @export
#'
#'
#'

glm_bas <- function(
  se,
  beta_prior = robust(),
  model_prior = beta.binomial(1,1),
  MCMC_iterations = 100000,
  decision_method = "BPM",
  top = NULL,
  save_dir = "LogisticModel",
  sample_names_coldata = colnames(colData(se))[1],
  group_list = group_list_default
){
if(!require(BAS)){
  stop("You need to first install BAS to use this function.")
}
#Create saving directories
if(dir.exists(save_dir)) {} else{
  dir.create(save_dir)
}
save_dir = gsub("/$","",save_dir)

#Check wheather the matrix in assay has correct values.
Target_matrix = assay( se )
if (any(!unique(as.vector(Target_matrix)) %in% c(NA,1,0))) stop("The assay matrix contains values other than 0,1,and NA")

#Build logistic models one by one.
if(!is.null(colnames(Target_matrix))) {
  idx_sample = colnames(Target_matrix)
} else {
  idx_sample = colData(se)[[sample_names_coldata]]
  if(anyNA(idx_sample)) stop("A index for the names of samples/coditions are required by either collumn names or first collumn of colData.")
}

if(any(duplicated(idx_sample))) stop("The sample index contains duplicated names.")

Features = SummarizedExperiment::mcols( se )

if(ncol(Features) == 0) stop("There are 0 features, please run predictors.annot first.")

colnames(Target_matrix) = idx_sample

for (i in idx_sample) {
Design =  as.data.frame( cbind(Target_matrix[,i],Features)[!is.na(Target_matrix[,i]),] )
colnames(Design)[1] = "Y"
Design$Y = Design$Y > 0
cat(paste0("Running model selection for collumn: ",i,"...\n"))

BAS_I =  BAS::bas.glm(Y ~ .,
                   family = binomial(link = "logit"),
                   data = Design,
                   n.models= 2^8,
                   betaprior = beta_prior,
                   method="MCMC",
                   MCMC.iterations = MCMC_iterations,
                   modelprior= model_prior)

MPIP <- data.frame( MPIP = BAS_I$probne0,
                   Covariate = c("Intercept",colnames(Features)) )

save_dir_i <- paste0(save_dir,"/",i)

if(dir.exists(save_dir_i)) {} else{dir.create(save_dir_i)}

write.csv(MPIP, paste0(save_dir_i,"/marginal_posterior_inclusion_prob.csv"),row.names = F)

indx_best_modelI <- predict(BAS_I, estimator = decision_method, top = top)$bestmodel

if (sum(indx_best_modelI) == 0) {

Final_moddel <- glm(Design$Y ~ 1, family = binomial(link = "logit"))

} else {

Final_moddel <- glm(Y ~ ., data = Design[,c(indx_best_modelI+1)], family = binomial(link = "logit"))

}

summary_glm <- summary(Final_moddel)

Deviance_df <- data.frame(Deviances = c(summary_glm$deviance,summary_glm$null.deviance),
                          Dof = c(summary_glm$df.residual,summary_glm$df.null) ,
                          Model = c("Residual","NULL"))

write.csv(Deviance_df, paste0(save_dir_i,"/Deviance_and_Dof.csv"),row.names = F)

Critical_value <- qnorm((.05/2)/length(indx_best_modelI),lower.tail = F)

plot_df <-  data.frame(summary_glm$coefficients)[,c("Estimate","z.value")]
plot_df$X_lab = gsub("_"," ", gsub("TRUE","", rownames(plot_df)))
if(plot_df$X_lab[1] == "(Intercept)") plot_df$X_lab[1] = "Intercept"
plot_df$X_lab = factor(plot_df$X_lab,levels = plot_df$X_lab[order(plot_df$Estimate,decreasing = F)])

plot_df <- reshape2::melt(plot_df,id.vars = "X_lab")
plot_df$Cv <- NA
t_idx <- which(plot_df$variable == "z.value")
plot_df$Cv[t_idx[1]] = Critical_value
plot_df$Cv[t_idx[2]] = -1*Critical_value

p1 <- ggplot(plot_df,aes(x = X_lab, y = value)) + geom_bar(stat = "identity", width = .4, fill = "red", colour = "red", size = 0.1) + geom_hline(aes(yintercept = Cv), alpha = .5, linetype = 2, size = .35) + coord_flip() + facet_grid(~variable,scales = "free") + theme_classic() + labs(title = "Logistic model on genomic features",subtitle = i, x = "predictors")

suppressWarnings( ggsave(paste0(save_dir_i,"/",i,"_logistic_model.pdf"), p1, width = 4.8, height = 1.7 +  (nrow(plot_df)/2)*.1 ) )

save_df <-  data.frame(summary_glm$coefficients)[,c("Estimate","z.value","Pr...z..")]
colnames( save_df ) = c("Estimate","z.statistics","pvalue")
save_df$Bf_adj_p = p.adjust(save_df$pvalue,method = "bonferroni")
write.csv(save_df,paste0(save_dir_i,"/","Statistics.csv"))
write.csv(plot_df,paste0(save_dir_i,"/","Plot.csv"))

}

#Integrate the data generated to visualize them.


#Construct Plot dataframe for meta plots

Plot_lst <- lapply(idx_sample,function(x) read.csv(paste0(save_dir,"/",x,"/Plot.csv")))

names( Plot_lst ) = idx_sample

Predictors <- c("Intercept" , colnames(Features) )

PLOT_DF = data.frame(
  Predictors = rep(Predictors,length(Plot_lst)),
  Condations = factor( rep(names( Plot_lst ),
                           each = length(Predictors)) ,
                       levels = idx_sample)
)

for(i in names(group_list)) {
  PLOT_DF$Group[PLOT_DF$Predictors %in% group_list[[i]]] = i
}

PLOT_DF$Group = factor(PLOT_DF$Group, levels = names(group_list))

PLOT_DF$Predictors = gsub("_"," ",PLOT_DF$Predictors)

PLOT_DF$Predictors = factor(PLOT_DF$Predictors, levels =  gsub("_"," ",unlist( group_list )))

PLOT_DF$z.statistics = 0

PLOT_DF$logit = 0

PLOT_DF$Cv = NA

PLOT_DF$MPIP = 0

for (i in names( Plot_lst )) {
  Plot_df_i <- Plot_lst[[i]]
  z_indx <- Plot_df_i$variable == "z.value"
  Plot_df_i2 <- Plot_df_i[z_indx,c("X_lab","value","Cv")]
  idx <- match(as.character( Plot_df_i2$X_lab ),  PLOT_DF$Predictors[PLOT_DF$Condations == i])
  idx_plot_df <- idx + min(which(PLOT_DF$Condations == i)) - 1
  PLOT_DF$z.statistics[idx_plot_df] <- Plot_df_i2$value
  PLOT_DF$logit[idx_plot_df] <- Plot_df_i$value[!z_indx]
  PLOT_DF$Cv[which(PLOT_DF$Condations == i)] <- Plot_df_i2$Cv[1]
  MPIP_i <- read.csv(paste0(save_dir,"/",i,"/","marginal_posterior_inclusion_prob.csv"))
  idx2 <- match(gsub("_"," ",MPIP_i$Covariate),PLOT_DF$Predictors[PLOT_DF$Condations == i])
  PLOT_DF$MPIP[idx2 + min(which(PLOT_DF$Condations == i)) - 1] <- MPIP_i$MPIP[idx2]
}

PLOT_DF$Direction = "Not selected"
PLOT_DF$Direction[PLOT_DF$z.statistics > 0] = "Positive"
PLOT_DF$Direction[PLOT_DF$z.statistics < 0] = "Negative"

PLOT_DF$Direction  = factor(PLOT_DF$Direction,levels = c("Positive","Negative","Not selected"))

PLOT_DF$ABS_Z = abs(PLOT_DF$z.statistics)

Z_STAT <- ggplot(PLOT_DF,aes(x = Predictors, y = ABS_Z)) + geom_bar(stat = "identity", width = .5, colour = 0, aes(fill = Direction)) + facet_grid(Group ~ Condations,scales = "free",space = "free_y") + coord_flip() + theme_linedraw() + geom_hline(aes(yintercept = Cv), alpha = .5, linetype = 2, size = .35) +
  theme(panel.grid.minor.x = element_line(linetype = 0),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        plot.margin = margin(1.5,1.5,1,1,"cm")) + labs(y = "Absolute values of Wald test z satistics", title = "Test statistics of logistic regression analysis") + scale_fill_manual(values = c("#1b9e77","#7570b3",NA))


ggsave(paste0(save_dir,"/","statistics-bar.pdf"),Z_STAT,width = 2 + 1.6*length(idx_sample),height = 1.5 + .2*ncol(Features))

##Plot logits

LOGIT <- ggplot(PLOT_DF,aes(x = Predictors, y = abs(logit))) + geom_bar(stat = "identity", width = .5, colour = 0, aes(fill = Direction)) + facet_grid(Group ~ Condations,scales = "free",space = "free_y") + coord_flip() + theme_linedraw() +
  theme(panel.grid.minor.x = element_line(linetype = 0),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        plot.margin = margin(1.5,1.5,1,1,"cm")) + labs(y = "Absolute values of logit estimates", title = "Logit estimates of logistic regression analysis") + scale_fill_manual(values = c("#d95f02","#7570b3",NA))

ggsave(paste0(save_dir,"/","effectsize-bar.pdf"), LOGIT, width = 2 + 1.6*length(idx_sample),height = 1.5 + .2*ncol(Features))

PLOT_DF$Predictors = factor(PLOT_DF$Predictors,levels = rev(levels( PLOT_DF$Predictors )))

##Plot marginal posterior inclusion probabilities (For model selection)
MPIP <- ggplot(PLOT_DF,aes(y = Predictors, x = Condations)) + geom_tile(aes(fill = MPIP)) + theme_linedraw() +
  theme(panel.grid.minor.x = element_line(linetype = 0),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.justification = "center",
        plot.margin = margin(1,2,1,.5,"cm")) + labs(x = "Conditions") + scale_fill_continuous(low = "#fff7bc", high = "#d95f0e")

ggsave(paste0(save_dir,"/","model-selection.pdf"), MPIP, width = 2.2 + .25*length(idx_sample),height = 1 + .2*ncol(Features))

##Plot model goodness of fit
Plot_df_Dev <- data.frame(Conditions = rep( idx_sample, 2))
Plot_df_Dev$Deviance = 0
Plot_df_Dev$Class = rep(c("Residual","Explained"),each = length(idx_sample))
Plot_df_Dev$Dof
Plot_df_Dev$lab_y_pos = 0

for (i in 1:length(idx_sample)) {
  GOF <- read.csv(paste0(save_dir,"/",idx_sample[i],"/","Deviance_and_Dof.csv"))
  Plot_df_Dev$Deviance[i] <- GOF$Deviances[1]
  Explained_deviance_i <- max(0, GOF$Deviances[2] - GOF$Deviances[1])
  Plot_df_Dev$Deviance[i+length(idx_sample)] <- Explained_deviance_i
  Plot_df_Dev$Dof[i] <- paste0(GOF$Dof[2]-GOF$Dof[1])
  Plot_df_Dev$Dof[i+length(idx_sample)] <- GOF$Dof[1]
  Plot_df_Dev$lab_y_pos[i] =  GOF$Deviances[1] + Explained_deviance_i/2
  Plot_df_Dev$lab_y_pos[i+length(idx_sample)] = GOF$Deviances[1]/2
}


Gof <- ggplot(Plot_df_Dev,aes(x = Conditions, label = Dof)) +
  geom_bar(stat = "identity",aes(fill = Class,y = Deviance), colour = "black", width = 1) +
  theme(axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1)) +
  geom_label(aes(y = lab_y_pos),size = 3) + scale_fill_brewer(palette = "Spectral") + labs(title = "Goodness of fits of the logistic models")

ggsave(paste0(save_dir,"/","Goodness-of-fit.pdf"), Gof, width = 1.8 + .5*length(idx_sample),height = 3.3)

cat("Summary plots are generated.\n")
}
