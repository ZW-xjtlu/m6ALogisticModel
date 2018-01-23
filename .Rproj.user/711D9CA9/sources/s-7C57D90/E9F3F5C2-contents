#' @title Generate predictors/features for a range based RNA modification data.
#'
#' @description \code{predictors.annot} is used to generate features given a \code{SummarizedExperiment} object of RNA modification / target.
#'
#' @details This function retreave transcript related features that are previous known to be related with m6A modifications based on
#' provided \code{rowRanges} of the \code{SummarizedExperiment},
#' and it return features in forms of meta data collums of the \code{SummarizedExperiment}.
#'
#' The features that must be included:
#'
#'
#'
#'   ###1. Transcript regions ### ---- The entries are logical / dummy variables.
#'
#' - UTR5: 5'UTR.
#'
#' - UTR3: 3'UTR.
#'
#' - CDS: Coding Sequence.
#'
#' - Stop_codons: Stop codon (301 bp center).
#'
#' - Start_codons: Start codon (201 bp center).
#'
#' - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
#'
#' - Exons: Exonic regions.
#'
#' - Last_exons_50bp: Start 50bp of the last exon of a transcript.
#'
#'
#'
#'   ###2. Relative positions ### ---- The entries fall into the scale of [0,1].  If the site is not mapped to any range on the right, the value is set to 0. (can be viewed as an interactive term on top of the region model.)
#'
#' - Pos_UTR5: Relative positioning on 5'UTR.
#'
#' - Pos_UTR3: Relative positioning on 3'UTR.
#'
#' - Pos_CDS: Relative positioning on Coding Sequence.
#'
#' - Pos_Tx: Relative positioning on Transcript.
#'
#' - Pos_exons: Relative positioning on exons.
#'
#'
#'    ###3. Region length ###
#'
#' - long_UTR3: Long 3'UTR (length > 400bp).
#'
#' - long_exon: Long exon (length > 400bp).
#'
#' - Gene_length_ex: standardized gene length of exonic regions (z score).
#'
#' - Gene_length_all: standardized gene length of all regions (z score).
#'
#'
#'
#'   #####=============== The following features that are optional ===============#####
#'
#'    ###4. Motif ###
#'
#' by default it includes the following motifs search c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"): i.e. instances of RRACH.
#'
#'    ###5. Evolutionary fitness ###
#'
#' - PC 1nt: standardized PC score 1 nt.
#'
#' - PC 201nt: standardized PC score 101 nt.
#'
#' - FC 1nt: standardized Fitness consequences scores 1nt.
#'
#' - FC 5nt: standardized Fitness consequences scores 101nt.
#'
#'    ###6. User specified features by argument \code{feature_lst} ###
#'
#' The entries are logical / dummy variables, specifying wheather overlapping with each GRanges or GRanges list.
#'
#'    ###7.Gene attribute ###
#'
#' - sncRNA: small noncoding RNA (<= 200bp)
#'
#' - lncRNA: long noncoding RNA (> 200bp)
#'
#' - Isoform_num: Transcript isoform numbers standardized by z score.
#'
#' - HK_genes: mapped to house keeping genes, such as defined by paper below.
#'
#' Eisenberg E, Levanon EY (October 2013). "Human housekeeping genes, revisited". Trends in Genetics. 29
#'
#'   ###7.Batch effect ###
#'
#' - GC_cont_genes: GC content of each gene.
#'
#' - GC_cont_101bp: GC content of 101bp local region of the sites.
#'
#'
#' @param se A \code{\link{SummarizedExperiment}} object containing the \code{rowRanges} for modifications. \code{colData} and \code{assay} are not neccessarily specified for this function.
#'
#' @param txdb \code{TxDb} object for annotating the corresponding \code{rowRanges}, this is either obtained from bioconductor or converted from the annotation files by \code{GenomicFeatures::makeTxDbFromGFF}.
#' @param bsgnm \code{\link{BSgenome}} object for genomic sequence annotation, this should be downloaded from bioconductor.
#' @param fc,pc Optional; \code{GScores} objects for annotations of standardized Fitness consequences scores and UCSC phastCons conservation scores.
#'
#' Gulko B, Melissa J. Hubisz, Gronau I and Siepel A (2015). “Probabilities of fitness consequences for point mutations across the human genome.” Nature Genetics, 47, pp. 276-283.
#'
#' Siepel A and al. e (2005). “Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes.” Genome Research, 15, pp. 1034-1050.
#'
#' @param struct_hybridize Optional; A \code{\link{GRanges}} or \code{\link{GRangesList}} object indicating the hybridized region on the transcribed or exonic regions.
#'
#' The precomputed MEA 2ndary structures could be find at the data attached in this package: \code{\link{Struc_hg19}} and \code{\link{Struc_mm10}}.
#'
#' @param feature_lst Optional; A list of \code{\link{GRanges}} for user defined features, the names of the list will correspond to the names of features.
#'
#' @param motif A character vector indicating the motifs centered by the modification nucleotite, the motif will not be attached if the \code{rowRanges} of \code{se} is not single nucleotide resolution (with all width = 1).
#'
#' By default, the motif selected is RRACH: c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC").
#'
#' @param HK_genes_list Optional; A character string of the Gene IDs of the House Keeping genes. The Gene IDs should correspond to the Gene IDs used by the provided \code{TxDb} object.
#'
#' The entrez gene IDs of the house keeping genes of mm10 and hg19 are included in this package: \code{HK_hg19} and \code{HK_mm10}.
#'
#' @return This function will return a \code{\link{SummarizedExperiment}} object with a \code{mcols} of a feature or design matrix.
#'
#' @examples
#' ### ==== For hg19 ==== ###
#'
#' library(SummarizedExperiment)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(fitCons.UCSC.hg19)
#' library(phastCons100way.UCSC.hg19)
#'
#' gr <- GRanges(seqnames =  'chr3',
#' IRanges(100:200, width=1),
#' strand='-')
#'
#' SE <- SummarizedExperiment()
#' rowRanges(SE) = gr
#'
#' Feature_lst_hg19 = list(
#' HNRNPC_eCLIP = readRDS("eCLIP_HNRNPC.rds"),
#' YTHDC1_TREW = readRDS("YTHDC1_TREW_gr.rds"),
#' YTHDF1_TREW = readRDS("YTHDF1_TREW_gr.rds"),
#' YTHDF2_TREW = readRDS("YTHDF2_TREW_gr.rds"),
#' miR_targeted_genes = readRDS("miR_targeted_genes_grl.rds"),
#' miRanda = readRDS("miRanda_hg19_gr.rds"),
#' TargetScan = readRDS("TargetScan_hg19_gr.rds"),
#' Verified_miRtargets = readRDS("verified_targets.rds")
#' )
#'
#' SE_features_added <- predictors.annot(se = SE,
#'                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                          bsgnm = Hsapiens,
#'                            fc = fitCons.UCSC.hg19,
#'                            pc = phastCons100way.UCSC.hg19,
#'                          struct_hybridize = Struc_hg19,
#'                        feature_lst = Feature_lst_hg19,
#'                      HK_genes_list = HK_hg19)
#'
#' mcols(SE_features_added) ###Check the generated feature matrix.
#'
#' @seealso \code{\link{logistic.modeling}} to perform model selection, statistics calculation, and visualization across multiple samples.
#'
#' @import BSgenome
#' @import dplyr
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @export


predictors.annot <- function(se,
                             txdb,
                             bsgnm,
                             fc = NULL,
                             pc = NULL,
                             struct_hybridize = NULL,
                             feature_lst = NULL,
                             motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"),
                             HK_genes_list = NULL
                             ) {
row_gr = rowRanges(se)
genes_txdb <- genes(txdb)
exbytx_txdb <- exonsBy(txdb,by = "tx")
exbg_txdb <- exonsBy(txdb,by = "gene")
exs_txdb <- exons(txdb)
Stop_codons <- resize( cds(txdb) , 1, fix = "end" )
Start_codons <- resize( cds(txdb) , 1, fix = "start" )
TSS <- transcripts(txdb) %>% resize(.,1,fix = "start")
A_idx <- vcountPattern("A", DNAStringSet( Views(Hsapiens, TSS))) > 0
UTR3 <- threeUTRsByTranscript(txdb)
UTR5 <- fiveUTRsByTranscript(txdb)

Feature_matrix = data.frame(UTR5 = row_gr%over%UTR5)

#Report Features
i = 1
Speak <- function(Fname,I) {cat( paste0("Feature " ,i," : ",Fname," is generated.\n"))
                        return(I+1)}
i  = Speak("5'utr",i)

Feature_matrix$UTR3 <- row_gr%over%UTR3

i  = Speak("3'utr",i)

Feature_matrix$CDS <- row_gr%over%cds(txdb)

i  = Speak("cds",i)

 # - Stop_codons: Stop codon (201 bp center).

Feature_matrix$Stop_codons <- row_gr%over%(Stop_codons + 100)

i  = Speak("stop codons 201bp",i)

 # - Start_codons: Start codon (201 bp center).
 #

Feature_matrix$Start_codons <- row_gr%over%(Start_codons + 100)
i  = Speak("start codons 201bp",i)

 # - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
 #
Feature_matrix$m6Am <- row_gr%over%TSS[A_idx]

i  = Speak("m6Am (Transcription strat site 1bp with A)",i)

 # - Exons: Exonic regions.
Feature_matrix$exons <- row_gr%over%exbg_txdb

i  = Speak("exon",i)

 #
 # - Last_exons_50bp: Start 50bp of the last exon of a transcript.
 #

exbytx_unlist <- unlist( exbytx_txdb )
ex_ranks <- exbytx_unlist$exon_rank
names(ex_ranks) = 1:length(ex_ranks)
Idx_last_exon <-  tapply( ex_ranks , names(exbytx_unlist), function(x) names(x)[which.max(x)]) %>% as.numeric
Last_exons_50bp <- exbytx_unlist[Idx_last_exon] %>% resize(.,50,fix = "start")
Feature_matrix$Last_exons_50bp <- row_gr%over%Last_exons_50bp

i  = Speak("Start 50bp of the last exon",i)

 #
 #
 #  ###2. Relative positions### ---- The entries fall into the scale of [0,1], if the entries not mapped to the range, the value is 0.
 #
 # - Pos_UTR5: Relative positioning on 5'UTR.

Feature_matrix$Pos_UTR5 <- Relative_Pos_map(row_gr,UTR5,0)

i  = Speak("relative positioning on 5'utr",i)

 #
 # - Pos_UTR3: Relative positioning on 3'UTR.
 #

Feature_matrix$Pos_UTR3 <- Relative_Pos_map(row_gr,UTR3,0)

i  = Speak("relative positioning on 3'utr",i)

 # - Pos_CDS: Relative positioning on Coding Sequence.

Feature_matrix$Pos_CDS <- Relative_Pos_map(row_gr,cdsBy(txdb,by = "tx"),0)

i  = Speak("relative positioning on cds",i)

 #
 # - Pos_Tx: Relative positioning on Transcript.
 #
 # - Pos_exons: Relative positioning on exons.

exs_grl <- GenomicRanges::split(exs_txdb,1:length(exs_txdb))
Feature_matrix$Pos_exons <- Relative_Pos_map(row_gr,exs_grl,0)
i  = Speak("relative positioning on exon",i)

 #
 #   ###3. Region length###
 #
 # - long_UTR3: Long 3'UTR (length > 400bp).

long_UTR3 <- UTR3[sum(width(UTR3)) > 400]
Feature_matrix$long_UTR3 <- row_gr%over%long_UTR3
i  = Speak("Long 3'UTR length > 400bp ",i)

 #
 # - long_exons: Long exon (length > 400bp).
 #

long_exs_txdb <- exs_txdb[width(exs_txdb) > 400]
Feature_matrix$long_exon <- row_gr%over%long_exs_txdb
i  = Speak("Long exon (length > 400bp)",i)

 # - Gene_length_ex: standardized gene length of exonic regions (z score).
 #

Feature_matrix$Gene_length_ex <- Properties_map(query_gr = row_gr,
                                                 feature_gr = exbg_txdb,
                                                 feature_properties = sum(width(exbg_txdb)),
                                                 no_map_val = 0,
                                                 normalize = T)


i  = Speak("gene length-exons (z-score)",i)

 #
 # - Gene_length_all: standardized gene length of all regions (z score).
 #
 #

Feature_matrix$Gene_length_full <-  Properties_map(query_gr = row_gr,
                                                                feature_gr = genes_txdb,
                                                                feature_properties = width(genes_txdb),
                                                     no_map_val = 0,
                                                   normalize = T)

i  = Speak("gene length full transcript (z-score)",i)

 #  #####=============== The following features that are optional ===============#####
 #
 #   ###4. Motif###
 #
 # by default it includes the following motifs search c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"): i.e. instances of RRACH.
 #

if (any( width(row_gr) != 1 )) {
  warning("At least 1 range with width > 1, the motifs are not attached.")
}else{

is_motif <- function(motif,dss,exact = T) vcountPattern(DNAString(motif), dss, fixed = exact) > 0
DSS <- DNAStringSet( Views(Hsapiens,row_gr + 2) )

for (m_i in motif)  {
Feature_matrix[[m_i]] <- is_motif(m_i,DSS,T)
i  = Speak(paste0("Motif: ",m_i),i)
}

}

 #   ###5. Evolutionary fitness###
 #
 # - PC 1nt: standardized PC score 1 nt.

if(!is.null(pc)) {
  Feature_matrix$PC_1nt <- scores(pc,row_gr)$scores

  Feature_matrix$PC_1nt[is.na(Feature_matrix$PC_1nt)] = mean(na.omit(Feature_matrix$PC_1nt))

  Feature_matrix$PC_1nt = Feature_matrix$PC_1nt >= 0.9

  i  = Speak("Phast cons scores 1nt >= 0.9",i)

  Feature_matrix$PC_101nt <- scores(pc,row_gr+50)$scores

  Feature_matrix$PC_101nt[is.na(Feature_matrix$PC_101nt)] = mean(na.omit(Feature_matrix$PC_101nt))

  Feature_matrix$PC_101nt = Feature_matrix$PC_101nt >= 0.8

  i  = Speak("Phast cons scores 101nt >= 0.8",i)
}

#Feature 21. fitness consequences scores 1nt.

if(!is.null(fc)) {
  suppressWarnings(  Feature_matrix$FC_1nt <- scores(fc,row_gr)$scores )

  Feature_matrix$FC_1nt[is.na(Feature_matrix$FC_1nt)] = mean(na.omit(Feature_matrix$FC_1nt))

  Feature_matrix$FC_1nt = Feature_matrix$FC_1nt >= .7

  i  = Speak("fitness consequences scores 1nt >= 0.7",i)

  suppressWarnings(  Feature_matrix$FC_101nt <- scores(fc,row_gr+50)$scores )

  Feature_matrix$FC_101nt[is.na(Feature_matrix$FC_101nt)] = mean(na.omit(Feature_matrix$FC_101nt))

  Feature_matrix$FC_101nt = Feature_matrix$FC_101nt >= .6

  i  = Speak("fitness consequences scores 101nt >= 0.6",i)
}

 #
 #   ###6. User specified features by argument \code{feature_lst}###
 #
 # The entries are logical / dummy variables, specifying wheather overlapping with each GRanges or GRanges list.

if(!is.null(feature_lst)) {
    for(i_flst in names(feature_lst)) {
        suppressWarnings( Feature_matrix[[i_flst]] <- row_gr %over% feature_lst[[i_flst]] )
      i  = Speak(paste0 ("annotation feature --- ", i_flst ),i)
 }
}

 #
 #   ###7.Gene attribute###
 #
 # - sncRNA: small noncoding RNA (<= 200bp)
 #

coding_tx <- names(cdsBy(txdb, by = "tx"))
exbytx_nc <- exbytx_txdb[!names(exbytx_txdb)%in%coding_tx]
lnc_idx <- sum(width(exbytx_nc)) > 200
Feature_matrix$sncRNA  <- row_gr%over%exbytx_nc[!lnc_idx]
i  = Speak("snc RNA (<= 200bp)",i)

 # - lncRNA: long noncoding RNA (> 200bp)
 #
Feature_matrix$lncRNA  <- row_gr%over%exbytx_nc[lnc_idx]
i  = Speak("snc RNA (> 200bp)",i)

 # - Isoform_num: Transcript isoform numbers standardized by z score.
 #

txbygenes <- transcriptsBy(txdb,by = "gene")
Feature_matrix$Isoform_num <- Properties_map(query_gr = row_gr,
                                              feature_gr = txbygenes,
                                              feature_properties =  pmin( elementNROWS(txbygenes), 20),
                                              no_map_val = 0,
                                             normalize = T)

i  = Speak("Isoform number z score",i)

 # - HK_genes: mapped to house keeping genes, such as defined by paper below.
 #
if(!is.null(HK_genes_list)){
  Feature_matrix$HK_genes <- row_gr%over%exbg_txdb[HK_genes_list]
  i  = Speak("House keeping genes",i)
}
 # Eisenberg E, Levanon EY (October 2013). "Human housekeeping genes, revisited". Trends in Genetics. 29
 #
 #  ###7.Batch effect###
 #
 # - GC_cont_genes: GC content of each gene.
exbg_gr <- unlist(exbg_txdb)
exs_GC_freq <-  letterFrequency( DNAStringSet( Views(Hsapiens, exbg_gr) ) , letters="CG")
Genes_GC_freq <- tapply( exs_GC_freq, names(exbg_gr), sum )
Genes_length <- tapply( width(exbg_gr), names(exbg_gr), sum )
Genes_GC_cont = Genes_GC_freq/Genes_length

Feature_matrix$GC_cont_genes <- Properties_map(query_gr = row_gr,
                                                feature_gr = exbg_txdb,
                                                feature_properties = Genes_GC_cont,
                                                no_map_val = .5,
                                               normalize = T)

i  = Speak("Gene level GC content z score",i)

 # - GC_cont_101bp: GC content of 101bp local region of the sites.

Feature_matrix$GC_cont_101bp <- as.numeric( letterFrequency( DNAStringSet( Views(Hsapiens,row_gr+50) ) , letters="CG", as.prob = T) )

Feature_matrix$GC_cont_101bp <- (Feature_matrix$GC_cont_101bp - mean(Feature_matrix$GC_cont_101bp))/sd(Feature_matrix$GC_cont_101bp)

i  = Speak("101bp GC content z score",i)

mcols(se) = Feature_matrix

return(se)
}
