#' @title RNA secondary structures filtering
#' @description \code{rBS_2ndStructure} is used to filter thermal stable RNA secondary structures from RNA bisulfite sequencing data. 
#'
#' @details The RNA secondary structures are predicted by RNAfold in ViennaRNA package. 
#' The structures are MEA secondary structures predicted under gammar = 0.1 and temperature = 70 degree.
#' 
#' The maximum pairing distance is set to be 150 nt. At the current version of the function, 
#' only structures on exons and mitochondria chromosome are included. Therefore, the sites not overlapped with exons or chrM will also be filtered.
#' 
#' @seealso To directly use the \code{GRangesList} of RNA secondary structures: \code{Struc_mm10},  \code{Struc_hg19}.
#' 
#' To generate RNA secondary structures with self specified genome sequence and transcript annotations:
#'  \code{\link{tx_seq_extraction}},  \code{\link{RNAfold}}, and  \code{\link{rfold_assembly_tx}}.
#' 
#' @param rBS_gr A \code{GRanges} object containing the methylation site information of RNA bisulfite sequencing.
#' 
#'  If the GRanges has a logical vector in its first column of the meta data column indicating the methylation states, 
#'  a report will be produced based on the association between the methylation and the RNA secondary structure.
#' @param Genome A character string which indicates the genome assembly used in rBS_gr, it can be either "mm10" or "hg19".
#' @param RNA_struc_dir Optinal; the directory of the self generated rds file containing the RNA secondary structures.
#' @param TXDB Optinal; the txdb object provided by the users.
#' @return A filtered GRanges object that has no overlapping with the predicted RNA structures.
#' The output ranges are also filtered by the overlapping with exons or the mitochondria chromosome. 
#' @examples
#' rBS_gr = GRanges(seqnames = rBS_df$`#SeqID`, strand = rBS_df$refStrand, ranges = IRanges(start = rBS_df$refPos, width = 1))
#' rBS_gr$mcols = p.adjust(rBS_df$`p-value_mState`,method = "BH") < .05
#' rBS_gr_filtered <- rBS_2ndStructure_Filter(rBS_gr,"hg19")
#' 
#' @export
rBS_2ndStructure_Filter <-
function(rBS_gr,Genome,RNA_struc_dir = NULL,TXDB = NULL) {
  
  if(!is.null(TXDB)) {
    exbytx <- exonsBy( TXDB )
  }
  
  if(!is.null(RNA_struc_dir)) { 
    Tx_struc <- readRDS(RNA_struc_dir)
  } else {
    Tx_struc <- eval(parse(text = paste0("Struc_",Genome)))
    exbytx <- eval(parse(text = paste0("exbytx_",Genome)))
    }
 
  idx_ex <- suppressWarnings(rBS_gr %over% exbytx) | (as.character(seqnames(rBS_gr)) == "chrM")

  rBS_gr = rBS_gr[idx_ex,]
  
  MEA_struc_70_ex <- suppressWarnings(rBS_gr %over% Tx_struc)
  
  if(ncol(mcols(rBS_gr)) >0) {
    if(is.logical(mcols(rBS_gr)[,1])) {
    
  fdr_mState = mcols(rBS_gr)[,1]
  
  chrM_idx <- which(as.logical(seqnames(rBS_gr) == 'chrM'))
  
  fisher_rs <- table(fdr_mState[-chrM_idx],MEA_struc_70_ex[-chrM_idx] ) %>% fisher.test
  
  fisher_M_rs <- table(fdr_mState[chrM_idx],MEA_struc_70_ex[chrM_idx]) %>% fisher.test
  
   cat( paste0("On nucleus chromosome:\n",
                "Fisher exact test p value: ",round(fisher_rs$p.value,4),"\n",
                "Odds ratio of enrichment: ",round(fisher_rs$estimate,4),"\n",
                "95% confidence interval of the enrichment: [",round(fisher_rs$conf.int[1],4),
                       ", ",round(fisher_rs$conf.int[2],4),"]","\n","\n",
                "On mitochondria chromosome:\n",
                "Fisher exact test p value: ",round(fisher_M_rs$p.value,4),"\n",
                "Odds ratio of enrichment: ",round(fisher_M_rs$estimate,4),"\n",
                "95% confidence interval of the enrichment: [",round(fisher_M_rs$conf.int[1],4),
                       ", ",round(fisher_M_rs$conf.int[2],4),"]")
  )
    }
  }
   
  return(rBS_gr[MEA_struc_70_ex])
}
