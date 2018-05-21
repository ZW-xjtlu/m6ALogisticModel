#' @title Build and visualize a logistic model to analyze the multinomial response varirable.
#'
#' @description \code{glm_multinomial} is a function used to build a multinomial logistic model on a given categorical Y.
#' @param Y A \code{vector} that defines the observed multinomial response variable.
#' @param PREDICTORS A \code{data.frame} that defines the model design, which will include the predictors / features in the collumns.
#' @param HDER The subtitle and the file name of the plot.
#' @param MODE Define the 0 outcomes of each logistic model, should be one of "against background","against pivot",
#' default is "against background".
#' @param CUT_OFF The cut off of the occurence of the less abundence class in binary features, if the less frequent class is less than this threshold, the feature will be dropped, default is 5. This is important when we want to have a reliable asymptotics result in Wald test.
#'
#' @details By default, K independent logistic models will be fit against the outcomes in Y that is not in the Kth outcome. In the mode of "against pivot", the reference level in the provided Y will be treated as the pivot / reference class.
#'
#' About deviance in GLM (from stackexchange):
#'
#' Let LL = loglikelihood
#' Here is a quick summary of what you see from the summary(glm.fit) output,
#' Null Deviance = 2(LL(Saturated Model) - LL(Null Model)) on df = df_Sat - df_Null
#' Residual Deviance = 2(LL(Saturated Model) - LL(Proposed Model)) df = df_Sat - df_Proposed
#' The Saturated Model is a model that assumes each data point has its own parameters (which means you have n parameters to estimate.)
#' The Null Model assumes the exact "opposite", in that is assumes one parameter for all of the data points, which means you only estimate 1 parameter.
#' The Proposed Model assumes you can explain your data points with p parameters + an intercept term, so you have p+1 parameters.
#' If your Null Deviance is really small, it means that the Null Model explains the data pretty well. Likewise with your Residual Deviance.
#' What does really small mean? If your model is "good" then your Deviance is approx Chi^2 with (df_sat - df_model) degrees of freedom.
#' If you want to compare you Null model with your Proposed model, then you can look at
#' (Null Deviance - Residual Deviance) approx Chi^2 with df Proposed - df Null = (n-(p+1))-(n-1)=p
#' Are the results you gave directly from R? They seem a little bit odd, because generally you should see that the degrees of freedom reported on the Null are always higher than the degrees of freedom reported on the Residual. That is because again, Null Deviance df = Saturated df - Null df = n-1 Residual Deviance df = Saturated df - Proposed df = n-(p+1)
#'
#' About multinomial logistic regression:
#'
#' Please check the \href{https://en.wikipedia.org/wiki/Multinomial_logistic_regression#As_a_set_of_independent_binary_regressions}{wikipedia page}
#'
#' Notice that I cannot sum the chisq statistics in background reference scenario, because the K models are dependent (not independent as in the pivot case).
#'
#' @import ggplot2
#' @import nnet
#' @importFrom reshape2 melt
#' @export

glm_multinomial <- function(Y,
                            PREDICTORS,
                            HDER = "multinomial",
                            MODE = c("against background","against pivot"),
                            CUT_OFF = 5) {

 stopifnot(length(Y) == nrow(PREDICTORS))
 MODE <- match.arg(MODE)

  #Fit multinomial GLM to analysis the clustering label

  indx_no_info <- sapply( PREDICTORS, function(x) { if (is.logical(x)){
    return( (sum(x) <= CUT_OFF | sum(!x) <= CUT_OFF) )
  } else {
    return(F)
  }
  } )

  if(any(indx_no_info)){
    warning(paste0("dropped dummy variable feature(s): ", paste0( gsub("TRUE","",names(indx_no_info[which(indx_no_info)])) , collapse=", "), "; the threshold is defined by the CUT_OFF argument (sum of TRUE / FALSE entries <= ",CUT_OFF,")."),call. = F, immediate. = T)
  }

  Model = PREDICTORS[,!indx_no_info]
  Model$Y = Y

  if( MODE == "against pivot") {

  Null_model <-  multinom(Y ~ 1, data = Model)
  Proposed_model <-  multinom(Y ~ ., data = Model)
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
  levels(Plot_df$FDR_sig) = c("< .05", "Insig")

  Proposed_model_chisq <- Null_model$deviance-Proposed_model$deviance
  Proposed_model_chisq = round(Proposed_model_chisq, 3)
  Model_dof <- length(Indx) * (K-1) - 1
  }

  if (MODE == "against background") {

    K = length(unique(Model$Y))
    Response_levels = as.character( unique(Model$Y))

    Plot_df_lst = list()

    Proposed_model_chisq = vector("numeric",length = K)
    names(Proposed_model_chisq ) = Response_levels
    Model_dof = vector("numeric",length = K)
    names(Model_dof ) = Response_levels

    for(i in Response_levels ) {
    Model_i = Model
    Model_i$Y = Model$Y == i
    logit_model_i <-  glm(Y ~ ., data = Model_i, family = binomial(link = "logit"))

    Proposed_model_chisq[i] <- logit_model_i$null.deviance - logit_model_i$deviance
    Model_dof[i] <- logit_model_i$df.null - logit_model_i$df.residual

    stat_test_i <- summary(logit_model_i)
    Plot_df_i <- data.frame(Clusters = i,
                            Covariates = rep(rownames(stat_test_i$coefficients),2),
                            values = c(stat_test_i$coefficients[,"Estimate"],stat_test_i$coefficients[,"z value"]),
                            Statistics = rep(c("log(OddsRatio)","Wald_Z"),each = nrow(stat_test_i$coefficients)),
                            pvalues = rep(stat_test_i$coefficients[,"Pr(>|z|)"],2)
                            )
    Plot_df_lst[[i]] <- Plot_df_i
    }

    Plot_df <- Reduce(rbind,Plot_df_lst)

    FDR_sig_indx <-  p.adjust( Plot_df$pvalues[Plot_df$Statistics == "Wald_Z"] , method = "fdr") < .05

    FDR_sig_indx <- unlist( lapply( split(FDR_sig_indx,rep(1:K,each = nrow(stat_test_i$coefficients))), function(x)rep(x,2)) )

    Plot_df$FDR_sig = "insig"

    Plot_df$FDR_sig[FDR_sig_indx] <- "fdr < .05"

    levels(Plot_df$FDR_sig) = c("fdr < .05", "insig")

    Plot_df$Covariates = gsub("TRUE","",Plot_df$Covariates)

    Indx <- names( sort( tapply(Plot_df$values[Plot_df$Statistics == "Wald_Z"] ,
                                Plot_df$Covariates[Plot_df$Statistics == "Wald_Z"],
                                function(x) sum(abs(x))) , decreasing = F) )

    Plot_df$Covariates = factor(Plot_df$Covariates, levels = Indx)

    Plot_df$values_abs = abs(Plot_df$values)
    Plot_df$Sign = "positive"
    Plot_df$Sign[Plot_df$values < 0] = "negative"
    Plot_df$Group = paste0(Plot_df$Statistics,":",Plot_df$Clusters)
    Indx_Group = unique(Plot_df$Group)

    n = length(Indx_Group)

    Plot_df$Group = factor( Plot_df$Group,
                            levels = Indx_Group )

    Proposed_model_chisq = round(Proposed_model_chisq, 3)
    Proposed_model_chisq = paste0( Proposed_model_chisq, collapse = ", ")
    Model_dof = paste0( Model_dof, collapse = ", ")
  }

  K = length(unique(Y))

  p <- ggplot(Plot_df) +
    geom_bar(stat = "identity",aes(x = Covariates, y = values_abs, fill = Sign, colour = FDR_sig), width = .75,size = .35, linetype = 1) +
    facet_wrap(~Group,nrow = 1,scales = "free_x") +
    theme_classic() + coord_flip() + theme(axis.text.y = element_text(size = 7),
                                           plot.margin = margin(1,1.5,1,1,"cm")) +
    scale_fill_brewer(direction = -1) + scale_colour_manual(values = c(ifelse(all(Plot_df$FDR_sig == "Insig"),0,"red"), 0)) +
    labs(title = paste0("Multinomial logistic model ",MODE,": ",HDER),
         subtitle = paste0("Proposed model deviance reduction: ",
                           Proposed_model_chisq,
                           " on ", Model_dof ," df"),
         y = "abs(values)")

  fig_width_p = 5 + 2.5*(K-1) + .01 * max(nchar(Indx))

  fig_height_p = 4 + .08 * length(Indx)

  ggsave(paste0(HDER,"_GLMestimates.pdf"),p,width = fig_width_p,height = fig_height_p)

  #Calculate the total statistical significance of the proposed model.
  #Stat_df <- data.frame(
  #  NULL_Deviance = Null_model$deviance,
  #  Residual_Deviance = Proposed_model$deviance
  #)

  #Stat_df$Reduced_prop = 1 - Stat_df$Residual_Deviance/Stat_df$NULL_Deviance
  #Stat_df$Cost_df = length(Indx) * (K-1) - 1
  #Stat_df$Chisq_stat = Stat_df$NULL_Deviance - Stat_df$Residual_Deviance
  #write.table(t(Stat_df),paste0("Model_report_",HDER,".txt"),col.names = F)


}
