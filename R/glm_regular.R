#' @title Build and visualize a generalized linear model.
#'
#' @description \code{glm_regular} is a function used to build a regular generalized linear model based on genomic features.
#' @param Y A \code{vector} that defines the response variable. It can be a 2-column integer matrix given the binomial family：the first column gives the number of successes and the second the number of failures.
#' @param PREDICTORS A \code{data.frame} that defines the model design, which will include the predictors / features in the collumns.
#' @param HDER The subtitle and the file name of the plot.
#' @param family Define the family of the glm, should be one of "gaussian", "binomial", and "poisson",
#' @param CUT_OFF The cut off of the occurence of the less abundence class in binary features, if the less frequent class is less than this threshold, the feature will be dropped, default is 5. This is important when we want to have a reliable asymptotics result in Wald test.
#' @param Exclude_intercept Whether to omit the intercept term when plot the estimates and statistics, this should be applied when the intercept estimates is too big relative to other predictors, default is FALSE.
#'
#' @details The function will fit a linear model based on the provided predictors and the response variable. The coefficient estimates and the wald test statistics will be save in a graph.
#' @import ggplot2
#' @export

glm_regular <- function(Y,
                        PREDICTORS,
                        HDER = "glm",
                        family = c("gaussian","binomial","poisson"),
                        CUT_OFF = 5,
                        Critical_value = 0.05,
                        Exclude_intercept = F) {

  stopifnot( is.null(ncol(Y)) | family == "binomial" )

  stopifnot( is.null(ncol(Y)) | ncol(Y) == 2 )

  stopifnot(length(Y) == nrow(PREDICTORS) | nrow(Y) == nrow(PREDICTORS))

  family <- match.arg(family)

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

  if( family == "gaussian") {

    Fitted_model <-  glm(Y~.,family = gaussian(link = "identity"),data = Model)

  }

  if( family == "binomial") {

    Fitted_model <- glm(Y~., family = binomial(link = "logit"),data = Model)

    }

  if( family == "poisson") {

    Fitted_model <- glm(Y~., family = poisson(link = "log"),data = Model)
  }

  summary_glm <- summary(Fitted_model)

  Critical_value <- qnorm((.05/2)/(ncol(Model)),lower.tail = F)

  plot_df <-  data.frame(summary_glm$coefficients)[,c("Estimate","z.value")]
  plot_df$X_lab = gsub("_"," ", gsub("TRUE","", rownames(plot_df)))
  if(plot_df$X_lab[1] == "(Intercept)") plot_df$X_lab[1] = "Intercept"

  if(Exclude_intercept) {
    plot_df =  plot_df[plot_df$X_lab != "Intercept",]
    plot_df$X_lab <- as.character(plot_df$X_lab)
  }

  plot_df$X_lab = factor(plot_df$X_lab,levels = plot_df$X_lab[order(plot_df$Estimate,decreasing = F)])

  plot_df <- reshape2::melt(plot_df,id.vars = "X_lab")
  plot_df$Cv <- NA
  t_idx <- which(plot_df$variable == "z.value")
  plot_df$Cv[t_idx[1]] = Critical_value
  plot_df$Cv[t_idx[2]] = -1*Critical_value

  library(ggplot2)

  p1 <- ggplot(plot_df,aes(x = X_lab, y = value)) +
    geom_bar(stat = "identity", width = .4, fill = "blue2", colour = "red", size = 0.1) +
    geom_hline(aes(yintercept = Cv), alpha = .5, linetype = 2, size = .35) +
    coord_flip() + facet_grid(~variable,scales = "free") +
    theme_classic() +
    labs(title = paste0(family, " linear model on genomic features" ),
         subtitle = HDER,
         x = "predictors")

  suppressWarnings( ggsave( paste0(HDER,"_",family,"_glm.pdf"), p1, width = 4.8, height = 1.7 +  (nrow(plot_df)/2)*.1 ))

}
