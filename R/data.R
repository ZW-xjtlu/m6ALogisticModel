#' MEA RNA 2ndary structures of mm10 exons predicted by RNAfold
#'
#' The full length transcripts are extracted from mm10;
#' the transcript annotation file was downloaded from the refSeq archive-2015-07-17-14-33-26.
#' The temperature of RNAfold prediction was set at 70 degree, and
#' the gammar for MEA structure was set at 0.1. The maximum paring distance was set at 150bp.
#'
#' For the transcripts longer than 8000bp, the predictions were conducted in windows of 2000bp with steps of 1000bp.
#'
#'
#'
#' @format A GRangesList object of length 34597:
#' \describe{
#'   Each element of the GRangesList represents a transcript, exept that the first element refers to the mitochondria chromosome.
#'   The GRanges within each of the GRangesList element represent the genomic locations of the MEA RNA 2ndary structures predicted by the method above.
#' }
#' @usage
#' rBS_gr = GRanges(seqnames = rBS_df$`#SeqID`, strand = rBS_df$refStrand, ranges = IRanges(start = rBS_df$refPos, width = 1))
#' rBS_gr_filtered <- rBS_gr[!rBS_gr \%over\% Struc_mm10]
"Struc_mm10"

#' MEA RNA 2ndary structures of hg19 exons predicted by RNAfold
#'
#' The full length transcripts are extracted from hg19;
#' the transcript annotation file was downloaded from the refSeq archive-2015-07-17-14-32-32.
#' The temperature of RNAfold prediction was set at 70 degree, and
#' the gammar for MEA structure was set at 0.1. The maximum paring distance was set at 150bp.
#'
#' For the transcripts longer than 8000bp, the predictions were conducted in windows of 2000bp with steps of 1000bp.
#'
#' @format A GRangesList object of length 53350:
#' \describe{
#'   Each element of the GRangesList represents a transcript, exept that the first element refers to the mitochondria chromosome.
#'   The GRanges within each of the GRangesList element represent the genomic locations of the MEA RNA 2ndary structures predicted by the method above.
#' }
#' @usage
#' rBS_gr = GRanges(seqnames = rBS_df$`#SeqID`, strand = rBS_df$refStrand, ranges = IRanges(start = rBS_df$refPos, width = 1))
#' rBS_gr_filtered <- rBS_gr[!rBS_gr \%over\% Struc_hg19]
"Struc_hg19"

#' An alternative group list
#'
#' Defined by the following Code:
#'
#'
#' group_list_expanded = list(
#' UTR5 = c("UTR5", "Pos_UTR5", "length_UTR5"),
#' CDS = c("CDS", "Pos_CDS", "length_CDS"),
#' UTR3 = c("UTR3", "Pos_UTR3", "length_UTR3"),
#' Exon = c("exons", "Pos_exons", "long_exon","Last_exons_50bp"),
#' Gene = c("Pos_Tx","length_gene_ex","length_gene_full","Isoform_num","sncRNA","lncRNA","HK_genes"),
#' LandMarks = c("m6Am","Start_codons","Stop_codons"),
#' Motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"),
#' Structure = c("struct_hybridize","struct_loop"),
#' Evolution = c("PC_1nt","PC_101nt","FC_1nt","FC_101nt"),
#' miRNA_RBP = c("HNRNPC_eCLIP", "YTHDC1_TREW", "YTHDF1_TREW", "YTHDF2_TREW", "miR_targeted_genes","TargetScan","Verified_miRtargets"),
#' Writers_Erasers = c("METTL3_TREW","METTL14_TREW","WTAP_TREW","METTL16_CLIP","ALKBH5_PARCLIP","FTO_CLIP","FTO_eCLIP"),
#' Batch = c("GC_cont_genes","GC_cont_101bp","Intercept")
#' )
#'
#' The list include an extra module called Writers_Erasers compared with \Code{group_list_default}, these data are published CLIP datasets of m6A writers and erasers.
#'
#' @format A List object.
#'
#' @usage
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
#'
#'SE_features_added <- predictors.annot(se = SE_example,
#'                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                      bsgnm = Hsapiens,
#'                                      fc = fitCons.UCSC.hg19,
#'                                      pc = phastCons100way.UCSC.hg19,
#'                                      struct_hybridize = Struc_hg19,
#'                                     feature_lst = Feature_List_expanded_hg19,
#'                                      HK_genes_list = HK_hg19_eids)
#'
#'
#'logistic.modeling(
#'  SE_features_added,
#'  MCMC_iterations = 100000,
#'  decision_method = "BPM",
#'  save_dir = "LogisticModel_X",
#'  sample_names_coldata = "ID",
#' group_list = group_list_expanded
#')
#'
#' @note For other provided list for grouping, please check \code{\link{group_list_default}}
"group_list_expanded"
