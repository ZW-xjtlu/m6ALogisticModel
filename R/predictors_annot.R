#' @title Generate predictors/features for a range based RNA modification data.
#'
#' @description \code{predictors_annot} is used to generate features given a \code{SummarizedExperiment} object of RNA modification / target.
#'
#' @details This function retreave transcript related features that are previous known to be related with m6A modifications based on
#' provided \code{rowRanges} of the \code{SummarizedExperiment},
#' and it return features in forms of meta data collums of the \code{SummarizedExperiment}.
#'
#' The features that must be included:
#'
#'
#'   ###1. Transcript regions ### ---- The entries are logical / dummy variables.
#'
#' - UTR5: 5'UTR.
#'
#' - UTR3: 3'UTR.
#'
#' - cds: Coding Sequence.
#'
#' - Stop_codons: Stop codon (301 bp center).
#'
#' - Start_codons: Start codon (201 bp center).
#'
#' - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
#'
#' - Exons: Exonic regions.
#'
#' - last_exons_50bp: Start 50bp of the last exon of a transcript.
#'
#'
#'
#'   ###2. Relative positions ### ---- The entries fall into the scale of [0,1].  If the site is not mapped to any range on the right, the value is set to 0. (can be viewed as an interactive term on top of the region model.)
#'
#' - pos_UTR5: Relative positioning on 5'UTR.
#'
#' - pos_UTR3: Relative positioning on 3'UTR.
#'
#' - pos_cds: Relative positioning on Coding Sequence.
#'
#' - pos_Tx: Relative positioning on Transcript.
#'
#' - pos_exons: Relative positioning on exons.
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
#' - PC 1bp: standardized PC score 1 nt.
#'
#' - PC 201bp: standardized PC score 101 nt.
#'
#' - FC 1bp: standardized Fitness consequences scores 1bp.
#'
#' - FC 5nt: standardized Fitness consequences scores 101bp.
#'
#'    ###6. User specified features by argument \code{feature_lst} ###
#'
#' The entries are logical / dummy variables, specifying whether overlapping with each GRanges or GRanges list.
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
#'
#' @param bsgnm \code{\link{BSgenome}} object for genomic sequence annotation, this should be downloaded from bioconductor.
#'
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
#' @param motif_clustering A character vector indicating the motif used to generate the features for the clustering indexes, Default: "DRACH".
#'
#' @param annot_clustering A \code{\link{GRanges}} object to generate clustering features. Default: NULL.
#'
#' The resulting clustering features will be named \code{clust_f100}, \code{clust_f1000}, \code{dist_nearest_p200}, and \code{dist_nearest_p2000}
#'
#' @param hk_genes_list Optional; A character string of the Gene IDs of the House Keeping genes. The Gene IDs should correspond to the Gene IDs used by the provided \code{TxDb} object.
#'
#' The entrez gene IDs of the house keeping genes of mm10 and hg19 are included in this package: \code{HK_hg19_eids} and \code{HK_mm10_eids}.
#'
#' @param isoform_ambiguity_method Can be "longest_tx" or "average". The former keeps only the longest transcript as the transcript annotation.
#' The later will use the average feature entries for multiple mapping of the transcript isoform.
#'
#' @param genes_ambiguity_method Can be "drop_overlap" or "average". The former will not annotate the modification sites overlapped with > 1 genes (By returning NA).
#' The later will use the average feature entries for mapping of multiple genes.
#'
#' @param standardization A logical indicating whether to standardize the continous features; Default TRUE.
#'
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
#' SE_features_added <- predictors_annot(se = SummarizedExperiment(rowRanges = hg19_miCLIP_gr),
#' txdb = txdb,
#' bsgnm = Hsapiens,
#' fc = fitCons.UCSC.hg19,
#' pc = phastCons100way.UCSC.hg19,
#' struct_hybridize = Struc_hg19,
#' feature_lst = Additional_features_hg19,
#' hk_genes_list = HK_hg19_eids,
#' motif = c("AAACA","AGACA","AAACT","AGACT","AAACC","AGACC",
#'           "GAACA","GGACA","GAACT","GGACT","GAACC","GGACC",
#'           "TAACA","TGACA","TAACT","TGACT","TAACC","TGACC"),
#' motif_clustering = "DRACH",
#' standardization = F,
#' genes_ambiguity_method = "average")
#'
#'
#' mcols(SE_features_added) ###Check the generated feature matrix.
#'
#'
#' @seealso \code{\link{glm_bas}}, \code{\link{glm_multinomial}}, \code{\link{glm_regular}} to perform model selection, statistics calculation, and visualization across multiple samples.
#'
#' @import BSgenome
#' @import Biostrings
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom GenomicScores score
#' @export


predictors_annot <- function(se,
                               txdb,
                                 bsgnm,
                                 fc = NULL,
                                 pc = NULL,
                                 struct_hybridize = NULL,
                                 feature_lst = NULL,
                                 motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"),
                                 motif_clustering = "DRACH",
                                 annot_clustering = NULL,
                                hk_genes_list = NULL,
                             isoform_ambiguity_method = c("longest_tx","average"),
                          genes_ambiguity_method = c("drop_overlap","average"),
                       standardization = TRUE
) {

  #Pre_check
  stopifnot(nrow(se)!=0)
  isoform_ambiguity_method <- match.arg(isoform_ambiguity_method)
  genes_ambiguity_method <- match.arg(genes_ambiguity_method)
  row_gr = rowRanges(se)
  stopifnot(length(row_gr)!=0)

  #Extract genomic features.

  genes_txdb <- GenomicFeatures::genes(txdb)
  exbytx_txdb <- exonsBy(txdb,by = "tx")

  #Handel the isoform ambiguity.

  if(isoform_ambiguity_method == "longest_tx"){
    Longest_tx_table <- find_longest_transcript(exbytx_txdb,txdb)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
  } else {
    Kept_tx_indx <- T
  }

  #Filter the transcripts into the longest ones by their genes
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]

  #Remove the overlapped transcripts that belong to multiple genes.
  if(genes_ambiguity_method == "drop_overlap") {
    exbytx_txdb <- exbytx_txdb[countOverlaps(exbytx_txdb,exbytx_txdb) == 1]
  }

  #PS. exbg_txdb is used independently of the isoform ambiguity (but the gene dis-ambiguity is still used).
  exbg_txdb <- exonsBy(txdb,by = "gene")

  if(genes_ambiguity_method == "drop_overlap") {
    keep_indx <- countOverlaps(exbg_txdb,exbg_txdb) == 1
    removed_exbg_txdb <- exbg_txdb[!keep_indx]
    exbg_txdb <- exbg_txdb[keep_indx]
    rm(keep_indx)
  }

  #Extract all the ambiguity removed exonic regions.
  exs_txdb <- unlist(exbytx_txdb)

  txid <- names(exbytx_txdb)

  cds <- cdsBy(txdb,by = "tx")

  cds <- cds[ names(cds) %in% txid ]

  Stop_codons <- resize( unlist( range(cds) ) , 3, fix = "end" )

  Start_codons <- resize( unlist( range(cds) ) , 3, fix = "start" )

  TSS <- resize(unlist(range(exbytx_txdb)),1,fix = "start")

  A_idx <- vcountPattern("A", DNAStringSet( Views(bsgnm, TSS))) > 0

  TSS <- resize(TSS,100,fix = "start")

  UTR3 <- threeUTRsByTranscript(txdb)

  UTR3 <- UTR3[ names(UTR3) %in% txid ]

  UTR5 <- fiveUTRsByTranscript(txdb)

  UTR5 <- UTR5[ names(UTR5) %in% txid ]

  Feature_matrix = data.frame(UTR5 = row_gr%over%UTR5)

  Intron <- intronsByTranscript(txdb)

  #Define a function for standardizing features
  scale_i <- function(x,stand = T){
    if(stand) return((x - mean(x))/sd(x))
     else return(x)
    }

  #Report Features
  i = 1

  Speak <- function(Fname,I) {cat( paste0("feature " ,i," : ",Fname," is generated.\n"))
    return(I+1)}

  i  = Speak("5'utr",i)

  Feature_matrix$UTR3 <- row_gr%over%UTR3

  i  = Speak("3'utr",i)

  Feature_matrix$cds <- row_gr%over%cds

  i  = Speak("cds",i)

  # - Stop_codons: Stop codon (201 bp center).

  Feature_matrix$Stop_codons <- row_gr%over%(Stop_codons + 100)

  i  = Speak("stop codons 203bp",i)

  # - Start_codons: Start codon (201 bp center).
  #

  Feature_matrix$Start_codons <- row_gr%over%(Start_codons + 100)
  i  = Speak("start codons 203bp",i)

  # - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
  #
  Feature_matrix$TSS <- row_gr%over%TSS

  i  = Speak("downstream 100bp of transcription start site",i)

  # - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
  #

  Feature_matrix$TSS_A <- row_gr%over%TSS[A_idx]

  i  = Speak("downstream 100bp of transcription start site with sequence A",i)

  # - Annotate various exonic regions.

  Feature_matrix$exon_stop <- row_gr%over%subsetByOverlaps( exbg_txdb, Stop_codons )
  i  = Speak("exon with stop codon",i)

  #Alternative exons: exons that could be introns in some transcripts
  exbytx_unlist <- unlist( exbytx_txdb )

  Feature_matrix$alternative_exon <- row_gr%over%subsetByOverlaps(exs_txdb, unlist(Intron)+1, type = "within", maxgap=0L)

  i  = Speak("alternatively spliced exons",i)

  #Constitutive exons: exons that always present in all transcripts

  Feature_matrix$constitutive_exon <- row_gr%over%subsetByOverlaps(exbytx_unlist, unlist(Intron)+1, type = "within", maxgap = 0L, invert = T )

  i  = Speak("constitutively spliced exons",i)


  #Internal exons: exons that are not the first or the last exons
  ex_ranks <- exbytx_unlist$exon_rank
  names(ex_ranks) = 1:length(ex_ranks)

  Idx_last_exon <- tapply( ex_ranks ,names(exbytx_unlist) ,function(x) names(x)[which.max(x)])
  Idx_last_exon <- as.numeric( Idx_last_exon[unique(names(exbytx_unlist))] )
  Indx_last_exon <- vector("logical",length(ex_ranks))
  Indx_last_exon[Idx_last_exon] <- T

  Indx_first_exon <- ex_ranks == 1

  Feature_matrix$internal_exon <- row_gr%over%exbytx_unlist[!(Indx_last_exon|Indx_first_exon)]

  i  = Speak("internal exons",i)

  #
  # - long_exons: Long exon (length > 400bp).
  #

  long_exs_txdb <- exs_txdb[width(exs_txdb) > 400]
  Feature_matrix$long_exon <- row_gr%over%long_exs_txdb
  i  = Speak("long exons (length > 400bp)",i)

  #
  ### below tries to annotate the features proposed by Ke's paper
  #

  last_exons <- exbytx_unlist[Indx_last_exon]

  last_exons_400bp <-  resize( last_exons[width(last_exons) >= 400],400,fix = "start")

  Feature_matrix$last_exon <- row_gr%over%last_exons

  i  = Speak("the last exon",i)

  Feature_matrix$last_exon_400bp <- row_gr%over%last_exons_400bp

  i  = Speak("5' start 400 bp of the last exon",i)

  Feature_matrix$last_exon_sc400 <- row_gr%over% subsetByOverlaps( last_exons_400bp, Stop_codons)

  i  = Speak("5' start 400 bp of the last exon including stop codons",i)

  Feature_matrix$intron <- row_gr%over%Intron[ names(Intron) %in% txid ]
  i  = Speak("intron",i)

  #
  #
  #  ###2. Relative positions### ---- The entries fall into the scale of [0,1], if the entries not mapped to the range, the value is 0.
  #
  # - pos_UTR5: Relative positioning on 5'UTR.

  Feature_matrix$pos_UTR5 <- relative_pos_map(row_gr,UTR5,0,standardization)

  i  = Speak("relative positioning on 5'utr",i)

  #
  # - pos_UTR3: Relative positioning on 3'UTR.
  #

  Feature_matrix$pos_UTR3 <- relative_pos_map(row_gr,UTR3,0,standardization)

  i  = Speak("relative positioning on 3'utr",i)

  # - pos_cds: Relative positioning on Coding Sequence.

  Feature_matrix$pos_cds <- relative_pos_map(row_gr,cds,0,standardization)

  i  = Speak("relative positioning on cds",i)

  #
  # - pos_Tx: Relative positioning on Transcript.


  # - pos_exons: Relative positioning on exons.

  exs_grl <- GenomicRanges::split(exs_txdb,1:length(exs_txdb))
  Feature_matrix$pos_exons <- relative_pos_map(row_gr,exs_grl,0,standardization)
  i  = Speak("relative positioning on exon",i)

  #
  # ### Splicing. ###
  #

  Splicing_junctions <- reduce( c(resize(exbytx_unlist,1,fix = "start"),
                          resize(exbytx_unlist,1,fix = "end")), min.gapwidth=0L)

  Splicing_junctions <- subsetByOverlaps( Splicing_junctions ,c(resize(TSS,1,fix = "start"),
                                          resize(unlist(range(exbytx_txdb)),1,fix = "end")), invert = T)

  Feature_matrix$dist_sj_5_p2000 <- distance_map(row_gr,
                                  Splicing_junctions,
                                  end = "five",
                                  maximum_distance = 2000,
                                  standardize = standardization)

  i  = Speak("distance to the upstream (5' end) splicing junction",i)

  Feature_matrix$dist_sj_3_p2000 <- distance_map(row_gr,
                                  Splicing_junctions,
                                  end = "three",
                                  maximum_distance = 2000,
                                  standardize = standardization)

  i  = Speak("distance to the downstream (3' end) splicing junction",i)

  #
  #   ###3. Region length###
  #
  # - long_UTR3: Long 3'UTR (length > 400bp).

  Feature_matrix$length_UTR3 <- properties_map(query_gr = row_gr,
                                               feature_gr = UTR3,
                                               feature_properties = sum(width(UTR3)),
                                               no_map_val = NA,
                                               normalize = standardization)

  Feature_matrix$length_UTR3[is.na(Feature_matrix$length_UTR3)] = 0

  i  = Speak("3'UTR length ",i)

  Feature_matrix$length_UTR5 <- properties_map(query_gr = row_gr,
                                               feature_gr = UTR5,
                                               feature_properties = sum(width(UTR5)),
                                               no_map_val = NA,
                                               normalize = standardization)

  Feature_matrix$length_UTR5[is.na(Feature_matrix$length_UTR5)] = 0

  i  = Speak("5'UTR length ",i)

  Feature_matrix$length_cds <- properties_map(query_gr = row_gr,
                                              feature_gr = cds,
                                              feature_properties = sum(width(cds)),
                                              no_map_val = NA,
                                              normalize = standardization)

  Feature_matrix$length_cds[is.na(Feature_matrix$length_cds)] = 0

  i  = Speak("cds length ",i)


  # - Gene_length_ex: standardized gene length of exonic regions (z score).
  #

  txbg_txdb <- transcriptsBy(txdb, by = "gene")

  Feature_matrix$length_gene_ex <- properties_map(query_gr = row_gr,
                                                  feature_gr = txbg_txdb,
                                                  feature_properties = sum(width(reduce(exbg_txdb))),
                                                  no_map_val = NA,
                                                  normalize = standardization)

  Feature_matrix$length_gene_ex[is.na(Feature_matrix$length_gene_ex)] = 0

  i  = Speak("gene length-exons" ,i)

  #
  # - Gene_length_all: standardized gene length of all regions (z score).
  #
  #

  Feature_matrix$length_gene_full <-  properties_map(query_gr = row_gr,
                                                     feature_gr = txbg_txdb,
                                                     feature_properties = sum(width(reduce(txbg_txdb))),
                                                     no_map_val = NA,
                                                     normalize = standardization)

  Feature_matrix$length_gene_full[is.na(Feature_matrix$length_gene_full)] = 0

  rm(txbg_txdb)

  i  = Speak("gene length full transcript",i)

  #
  # - length_tx_exon: the length of the transcripts exons
  #
  #

  tx <- transcripts(txdb)

  Feature_matrix$length_tx_exon <-  properties_map(query_gr = row_gr,
                                                     feature_gr = tx,
                                                     feature_properties = sum(width(exbytx_txdb)),
                                                     no_map_val = NA,
                                                     normalize = standardization)

  Feature_matrix$length_tx_exon[is.na(Feature_matrix$length_tx_exon)] = 0

  i  = Speak("transcript exonic length",i)

  #
  # - length_tx_full: the length of the transcripts
  #
  #

  Feature_matrix$length_tx_full <-  properties_map(query_gr = row_gr,
                                                   feature_gr = tx,
                                                   feature_properties = width(tx),
                                                   no_map_val = NA,
                                                   normalize = standardization)

  Feature_matrix$length_tx_full[is.na(Feature_matrix$length_tx_full)] = 0

  i  = Speak("transcript full length",i)

  rm(tx)

  #  #####=============== The following features that are optional ===============#####
  #
  #   ###4. Motif###
  #
  # by default it includes the following motifs search c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"): i.e. instances of RRACH.
  #

  if (any( width(row_gr) != 1 )) {
    warning("At least 1 range with width > 1, the motifs are not attached.")
  }else{

    is_motif <- function(motif,dss,exact = F) vcountPattern(DNAString(motif), dss, fixed = exact) > 0
    DSS <- DNAStringSet( Views(bsgnm,row_gr + 2) )

    for (m_i in motif)  {
      Feature_matrix[[m_i]] <- is_motif(m_i,DSS,F)
      i  = Speak(paste0("motif --- ",m_i),i)
    }

  }

  #
  # ## Clustering effect.
  #
  if(!is.null(annot_clustering)) {

  self_count <- countOverlaps( row_gr , annot_clustering)

  Feature_matrix$clust_f1000 <- as.numeric( scale_i(countOverlaps( row_gr + 1000 , annot_clustering) - self_count , standardization) )

  i  = Speak("clustering indicators --- number of neighboors within 1000bp flanking regions",i)

  Feature_matrix$clust_f100 <- as.numeric( scale_i(countOverlaps( row_gr + 100 , annot_clustering) - self_count , standardization) )

  i  = Speak("clustering indicators ---  number of neighboors within 100bp flanking regions",i)

  rm(self_count)

  dist_self <- rep(2000, length(row_gr))

  match_obj <- distanceToNearest(row_gr, annot_clustering)

  dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance

  rm(match_obj)

  #For row_gr that is overlapping with annot_clustering
  if(any(row_gr %over% annot_clustering)){

  match_obj <- distanceToNearest(annot_clustering)

  fol <- findOverlaps(row_gr, annot_clustering)

  dist_tmp <- rep(2000, queryLength(fol))

  dist_tmp <- mcols(match_obj)$distance[match(subjectHits(fol),queryHits(match_obj))]

  dist_self[queryHits(fol)] <- dist_tmp

  rm(fol, match_obj, dist_tmp)

  }

  Feature_matrix$dist_nearest_p2000 <- as.numeric( scale_i( pmin(dist_self, 2000) , standardization) )

  i  = Speak("clustering indicators --- distance to the nearest neigboors (peaked at 2000bp)", i)

  Feature_matrix$dist_nearest_p200 <- as.numeric( scale_i( pmin(dist_self,200) , standardization) )

  i  = Speak("clustering indicators --- distance to the nearest neigboors (peaked at 200bp)", i)

  rm(dist_self)

  }

  if(!is.null(motif_clustering)) {

    tx_reduced <- reduce(transcripts(txdb))

    row_gr_expanded <- reduce(row_gr+2000)

    motif_regions <- intersect(row_gr_expanded,tx_reduced)

    rm(row_gr_expanded, tx_reduced)

    motif_transcripts <- sample_sequence( motif_clustering,
                                          motif_regions,
                                          bsgnm )

    rm(motif_regions)

    Feature_matrix[[paste0("clust_", motif_clustering, "_f1000")]] <- as.numeric( scale_i(countOverlaps( row_gr+1000 , motif_transcripts) , standardization) )

    i  = Speak(paste0("motif clustering --- number of ", motif_clustering, " motif neighboors within 1000bp flanking regions"),i)

    Feature_matrix[[paste0("clust_", motif_clustering, "_f100")]] <- as.numeric( scale_i(countOverlaps( row_gr+100 , motif_transcripts) , standardization) )

    i  = Speak(paste0("motif clustering ---  number of ", motif_clustering, " motif neighboors within 100bp flanking regions"),i)

    motif_transcripts <- subsetByOverlaps( motif_transcripts, row_gr, invert = T) #Remove all the self motifs

    match_obj <- distanceToNearest(row_gr, motif_transcripts)

    dist_motif <- rep(2000,length(row_gr))

    dist_motif[queryHits(match_obj)] <- mcols(match_obj)$distance

    dist_self <- rep(2000,length(row_gr))

    match_obj <- distanceToNearest(row_gr)

    dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance

    less_index <- dist_self < dist_motif

    dist_motif[less_index] <- dist_self[less_index]

    Feature_matrix[[paste0("dist_",motif_clustering,"_p2000")]] <- as.numeric( scale_i( pmin(dist_motif,2000) , standardization))

    i  = Speak(paste0("motif clustering --- distance to the nearest ",motif_clustering," motif (peaked at 2000bp)"),i)

    Feature_matrix[[paste0("dist_",motif_clustering,"_p200")]] <- as.numeric( scale_i( pmin(dist_motif,200) , standardization))

    i  = Speak(paste0("motif clustering --- distance to the nearest ",motif_clustering," motif (peaked at 200bp)"),i)

    rm(match_obj)
    rm(dist_self)
    rm(dist_motif)
    rm(less_index)
    rm(motif_transcripts)
  }

  #   ###5. Evolutionary fitness###
  #
  # - PC 1bp: standardized PC score 1 nt.

  if(!is.null(pc)) {
    Feature_matrix$PC_1bp <- score(pc,row_gr)

    Feature_matrix$PC_1bp[is.na(Feature_matrix$PC_1bp)] = mean(na.omit(Feature_matrix$PC_1bp))

    Feature_matrix$PC_1bp = as.numeric( scale_i( Feature_matrix$PC_1bp , standardization) )

    i  = Speak("phast cons scores 1bp",i)

    Feature_matrix$PC_101bp <- score(pc,row_gr+50)

    Feature_matrix$PC_101bp[is.na(Feature_matrix$PC_101bp)] = mean(na.omit(Feature_matrix$PC_101bp))

    Feature_matrix$PC_101bp = as.numeric( scale_i( Feature_matrix$PC_101bp , standardization) )

    i  = Speak("phast cons scores 101bp",i)
  }

  #Feature 21. fitness consequences scores 1bp.

  if(!is.null(fc)) {
    suppressWarnings(  Feature_matrix$FC_1bp <- score(fc,row_gr) )

    Feature_matrix$FC_1bp[is.na(Feature_matrix$FC_1bp)] = mean(na.omit(Feature_matrix$FC_1bp))

    Feature_matrix$FC_1bp = as.numeric( scale_i( Feature_matrix$FC_1bp , standardization) )

    i  = Speak("fitness consequences scores 1bp z score",i)

    suppressWarnings(  Feature_matrix$FC_101bp <- score(fc,row_gr+50) )

    Feature_matrix$FC_101bp[is.na(Feature_matrix$FC_101bp)] = mean(na.omit(Feature_matrix$FC_101bp))

    Feature_matrix$FC_101bp =  as.numeric( scale_i( Feature_matrix$FC_101bp , standardization) )

    i  = Speak("fitness consequences scores 101bp z score",i)
  }

  if(is.null(struct_hybridize)) {} else {

    Feature_matrix$struct_hybridize <- row_gr%over%struct_hybridize

    i  = Speak("RNA structure --- predicted hybridized region",i)

    Feature_matrix$struct_loop <- row_gr%over%infer_loop( struct_hybridize )

    i  = Speak("RNA structure --- inferred loop structures between hybridized region",i)

  }

  #
  #   ###6. User specified features by argument \code{feature_lst}###
  #
  # The entries are logical / dummy variables, specifying whether overlapping with each GRanges or GRanges list.

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

  coding_tx <- names(cds)
  exbytx_nc <- exbytx_txdb[!names(exbytx_txdb)%in%coding_tx]
  lnc_idx <- sum(width(exbytx_nc)) > 200
  Feature_matrix$sncRNA  <- row_gr%over%exbytx_nc[!lnc_idx]
  i  = Speak("snc RNA (<= 200bp)",i)

  # - lncRNA: long noncoding RNAs (> 200bp)
  #
  Feature_matrix$lncRNA  <- row_gr%over%exbytx_nc[lnc_idx]
  i  = Speak("lnc RNA (> 200bp)",i)

  # -lincRNA: Long intervening/intergenic noncoding RNAs.

  Feature_matrix$lncRNA  <- row_gr%over%subsetByOverlaps(exbytx_nc[lnc_idx],
                                                         subsetByOverlaps( genes_txdb, cds ),
                                                         invert = TRUE)

  i  = Speak("linc RNA (> 200bp) lncRNA which do not overlap protein-coding genes",i)

  # - Isoform_num: Transcript isoform numbers standardized by z score.
  #

  txbygenes <- transcriptsBy(txdb,by = "gene")
  Feature_matrix$isoform_num <- properties_map(query_gr = row_gr,
                                               feature_gr = txbygenes,
                                               feature_properties =  pmin( elementNROWS(txbygenes), 20),
                                               no_map_val = 0,
                                               normalize = standardization)

  i  = Speak("isoform number z score",i)


  Feature_matrix$exon_num <- properties_map(query_gr = row_gr,
                                               feature_gr = exbytx_txdb,
                                               feature_properties =  elementNROWS(exbytx_txdb),
                                               no_map_val = 0,
                                               normalize = standardization)

  i  = Speak("exon number z score",i)

  # - HK_genes: mapped to house keeping genes, such as defined by paper below.
  #
  if(!is.null(hk_genes_list)){
    Feature_matrix$HK_genes <- row_gr%over%exbg_txdb[names(exbg_txdb) %in% hk_genes_list]
    i  = Speak("house keeping genes",i)
  }

  # Eisenberg E, Levanon EY (October 2013). "Human housekeeping genes, revisited". Trends in Genetics. 29
  #
  #  ###7.Batch effect###
  #
  # - GC_cont_genes: GC content of each gene.
  exbg_gr <- unlist(exbg_txdb)
  exs_GC_freq <-  letterFrequency( DNAStringSet( Views(bsgnm, exbg_gr) ) , letters="CG")
  Genes_GC_freq <- tapply( exs_GC_freq, names(exbg_gr), sum )
  Genes_length <- tapply( width(exbg_gr), names(exbg_gr), sum )
  Genes_GC_cont = Genes_GC_freq/Genes_length

  Feature_matrix$GC_cont_genes <- properties_map(query_gr = row_gr,
                                                 feature_gr = exbg_txdb,
                                                 feature_properties = Genes_GC_cont,
                                                 no_map_val = .5,
                                                 normalize = standardization)

  i  = Speak("gene level GC content z score",i)

  # - GC_cont_101bp: GC content of 101bp local region of the sites.

  Feature_matrix$GC_cont_101bp <- as.numeric( letterFrequency( DNAStringSet( Views(bsgnm,row_gr+50) ) , letters="CG", as.prob = T) )

  if(standardization) Feature_matrix$GC_cont_101bp <- (Feature_matrix$GC_cont_101bp - mean(Feature_matrix$GC_cont_101bp))/sd(Feature_matrix$GC_cont_101bp)

  i  = Speak("101bp GC content z score",i)

  if(standardization){

  # - GC_cont_101bp_abs: Absolute value of the GC content of 101bp local region of the sites.

  Feature_matrix$GC_cont_101bp_abs <- abs(Feature_matrix$GC_cont_101bp)

  i  = Speak("absolute value of the 101bp GC content z score",i)

  }

  #Finally, mask the feature values that mapped to ambiguious genes as NA.

  if(genes_ambiguity_method == "drop_overlap") {
  Feature_matrix[row_gr %over% removed_exbg_txdb,] = NA
  }
  mcols(se) = Feature_matrix

  return(se)
}
