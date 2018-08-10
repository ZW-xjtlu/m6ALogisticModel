#' @title Run GO analysis over methylation sets with multi-class labeling.
#'
#' @description \code{go_multinomial} is a function used to run GO analysis on multi-class labeled ranged data. The back ground used is the shared genes in all of the ranges.
#' @param Y A \code{vector} that defines the observed multinomial response variable.
#' @param row_ranges A \code{GRanges} object that defines the genomic location of each sites on genome.
#' @param HDER A \code{character} defining the subtitle and the file name of the plot.
#' @param txdb A \code{txdb} object downloaded from bioconductor.
#' @param orgDb An \code{OrgDb} object defined by AnnotationDbi.
#' @param category A character specifying the gene ontology category, can be one of "BP", "CC", and "MF", default "BP".
#' @param gene_key A character specifying the type of the gene ID, the available types of the keys can be find using \code{keytypes(org.Hs.eg.db)}, default "ENTREZID".
#' @param min_bg_count Term minimum number of occurence in background genes; default 10.
#' @param max_bg_count Term maximum number of occurence in background genes; default 500.
#' @param min_gs_count Term minimum number of occurence in gene set genes; default 10.
#' @param max_gs_count Term maximum number of occurence in gene set genes; default 500.
#' @param pvalue_correction Method used for multiple hypothesis adjustment, can be one in "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", and "none".
#' @param GO_Slim Whether to conduct enrichment analysis only on GO slim terms (a certain subset to GO terms), default FALSE.
#' @param EASE_Score Whether or not use EASE score method. (a more conservative hypergeomatric test calculation used in DAVID) for more details please refer to \href{https://david.ncifcrf.gov/helps/functional_annotation.html#fisher}{here}, default TRUE.
#' @param top_terms Top number of terms (sorted by p values) for each multinomial class returned in the analysis, default is 10.
#' @param save_chart Whether to save the chart of GO enrichment for each calss, default FALSE.
#' @param background_ranges A \code{GRanges} object that defines the GO enrichment background, default is not used.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @import GenomicFeatures
#' @export

go_multinomial <- function(Y,
                           row_ranges,
                           txdb,
                           orgDb,
                           HDER = "GO_multinomial",
                           category = "BP",
                           gene_key = "ENTREZID",
                           min_bg_count = 10,
                           min_gs_count = 10,
                           max_bg_count = 500,
                           max_gs_count = 500,
                           pvalue_correction = "BH",
                           GO_Slim = FALSE,
                           EASE_Score = TRUE,
                           top_terms = 10,
                           save_chart = FALSE,
                           background_ranges = NULL
                           ){

stopifnot(length(Y) == length(row_ranges))

stopifnot(is.vector(Y))

if(!require(golite)){
  stop("You need to install package 'golite' first. You can consider to use the command: devtools::install_github('/ZhenWei10/golite').")
}

gene_by_tx <- transcriptsBy(txdb,by = "gene")

genes_lst <- tapply(row_ranges, Y, function(x) names(subsetByOverlaps(gene_by_tx,x)) )

if(!is.null(background_ranges)) {

genes_bg <- names( subsetByOverlaps(gene_by_tx, background_ranges) )

stopifnot(all( unique(unlist(genes_lst)) %in% genes_bg ))

} else {

genes_bg <-  unique( unlist(genes_lst) )

}


goea_lst <- goea(gene_set = genes_lst,
                 back_ground = genes_bg,
                 orgDb = orgDb,
                 category = category,
                 gene_key = gene_key,
                 min_bg_count = min_bg_count,
                 max_bg_count = max_bg_count,
                 min_gs_count = min_gs_count,
                 max_gs_count = max_gs_count,
                 EASE_Score = EASE_Score,
                 pvalue_correction = pvalue_correction,
                 interpret_term = TRUE,
                 show_gene_name = F,
                 GO_Slim = GO_Slim,
                 Slim_ss = NULL,
                 Exclude_self = T)

if(save_chart) {
for(i in names(goea_lst)) {
write.csv(goea_lst[[i]],paste0(HDER,"_chart_",i,".csv"))
}
}

goea_lst_sub <- lapply(goea_lst, function(x) x[1: min(top_terms,nrow(x)),])

Plot_df <- Reduce( rbind, goea_lst_sub )

unique_terms <- rev(unique(as.character(Plot_df$definition)))

Plot_df$definition = factor(Plot_df$definition , levels = unique_terms)

Plot_df$Clusters <- rep(names(goea_lst_sub),each = top_terms)

#Add the unselected rows

goea_lst_sub_comp <- Map(function(x,y) {
 return( y[(y$definition %in% unique_terms) & !(y$definition %in% x),] )
},  lapply(goea_lst_sub, function(x) x$definition), goea_lst)

Plot_df_sup <- Reduce( rbind, goea_lst_sub_comp )
Plot_df_sup$Clusters <-  rep(names(goea_lst_sub_comp),elementNROWS(goea_lst_sub_comp))

if(is.factor(Y)) {Plot_df_sup$Clusters = factor(Plot_df_sup$Clusters,
                                                levels = levels(Y))}

Plot_df_sup$definition = factor(Plot_df_sup$definition , levels = unique_terms)
Plot_df <- rbind(Plot_df,Plot_df_sup)

Plot_df$nlog2padj <- -1 * log2 (Plot_df[[paste0("adj_",pvalue_correction)]])

Plot_df$stat_sig <- "Insig"

Plot_df$stat_sig[Plot_df[[paste0("adj_",pvalue_correction)]] < .05] <- "padj < .05"

colnames(Plot_df)[colnames(Plot_df) == "nlog2padj"] = "-log2 padj"

mean_padj <- round( tapply(Plot_df$`-log2 padj`,Plot_df$Clusters,mean), 2 )


p <- ggplot(Plot_df) + geom_point(aes(x = Clusters,
                                 y = definition,
                                 size = `-log2 padj`,
                                 colour = Clusters,
                                 shape = stat_sig)) + theme_bw() +
               labs(x = "Clusters", y = paste0("GO ",category," terms"),
                    title = paste0("GO ",ifelse(GO_Slim,"slim ",""),"enrichment analysis"),
                    subtitle = paste0("mean(top",top_terms ," -log2 padj): ",paste0(mean_padj, collapse = ", "))) +
                theme(axis.text.y = element_text(size = 7),
                      plot.margin = margin(1,1.5,1,1,"cm")) +
               scale_shape_manual(values = c(16,18))

if (length(unique(Y)) <= 8) { p <- p + scale_colour_brewer(palette = "Dark2") }

K = length(goea_lst)

fig_width_p = 2 + 1*K + .033 * max(nchar(unique_terms))

fig_height_p = 4 + .13 * length(unique_terms)

ggsave(paste0(HDER,"_goea.pdf"),p,width = fig_width_p,height = fig_height_p)

return(p)
}
