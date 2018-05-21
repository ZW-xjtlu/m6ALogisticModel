#A new package: m6ALogisticModel

#2 main functions

#function 1: Predictors_annotation.
#Input: Summarized experiment. Output: Summarized experiment with annotated mcols.
#- Improvement, Some of the non-dummy data features are rescaled into normalized range. (z - scores)
#- The Pos features should be not re-scaled by z, since assuming the begining and end of region must have reversed logit, and suggesting the mid point as the splitting point is not making biological sense in this domain.
#Length of genes/3'UTR/Isoform number/GC content (?) should be rescaled into Z (mean 0), SD 1.

#function 2: m6ALogisticModel.

#Motivation to use logistic model: the transcript features are usually confounded with each other.
#Example 1. TSSA and Start codon.
#Example 2. Start of the last exon and Stop codon.
#Example 3. Long exon and Stop codon.
#Example 4. Long exon and Long genes.
#Example 5. GC content, structure, and transcript start.

#Motivation to use asymetrical distance: the negative and negative match is less important than target target match, since the negative match may only contributed by the sgared Cell line expression.

#Below are used while developing predictors.annot

library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

se_combinded <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/D_InfMerge_2017_12_17/se_combinded.rds")

rowRanges(se_combinded) = resize(rowRanges(se_combinded),1,fix = "center")

setwd("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4")

Feature_List_hg19 <- list(
    HNRNPC_eCLIP = readRDS("eCLIP_HNRNPC.rds"),
    YTHDC1_TREW = readRDS("YTHDC1_TREW_gr.rds"),
    YTHDF1_TREW = readRDS("YTHDF1_TREW_gr.rds"),
    YTHDF2_TREW = readRDS("YTHDF2_TREW_gr.rds"),
    miR_targeted_genes = readRDS("miR_targeted_genes_grl.rds"),
    miRanda = readRDS("miRanda_hg19_gr.rds"),
    TargetScan = readRDS("TargetScan_hg19_gr.rds"),
    Verified_miRtargets = readRDS("verified_targets.rds")
  )

HK_hg19_eids = names(readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/HK_genes_gr.rds"))

library(m6ALogisticModel)

SE_features_added <- predictors.annot(se = se_combinded,
                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      bsgnm = Hsapiens,
                                      fc = fitCons.UCSC.hg19,
                                      pc = phastCons100way.UCSC.hg19,
                                      struct_hybridize = rBS2ndStructure::Struc_hg19,
                                      feature_lst = Feature_List_hg19,
                                      HK_genes_list = HK_hg19_eids)

se = se_combinded
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
bsgnm = Hsapiens
fc = fitCons.UCSC.hg19
pc = phastCons100way.UCSC.hg19
struct_hybridize =  rBS2ndStructure::Struc_hg19
feature_lst = Feature_List_hg19
HK_genes_list = HK_hg19_eids
motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC")



###Code part


genes_txdb <- genes(txdb)
Feature_matrix = data.frame(Gene_length_all = Properties_mapp(query_gr = row_gr,
                                                              feature_gr = genes_txdb,
                                                              feature_properties = width(genes_txdb),
                                                              no_map_val = 0))

cat("Feature 1: Gene length-All is attached.\n")


#Feature 2. Gene length exons.
exbg_txdb <- exonsBy(txdb,by = "gene")
Feature_matrix$Gene_length_ex <- Properties_mapp(query_gr = row_gr,
                                                 feature_gr = exbg_txdb,
                                                 feature_properties = sum(width(exbg_txdb)),
                                                 no_map_val = 0)
cat("Feature 2: Gene length-exons is attached.\n")


#Feature 3. 5'Cap m6Am (TSS 1nt with A).
TSS <- transcripts(txdb) %>% resize(.,1,fix = "start")

A_idx <- vcountPattern("A", DNAStringSet( Views(Hsapiens, TSS))) > 0

Feature_matrix$m6Am <- row_gr%over%TSS[A_idx]

cat("Feature 3:  5'Cap m6Am is attached.\n")

#Feature 4. Start 50bp of the last exon.
exbytx <- exonsBy(txdb,by = "tx")
exbytx_unlist <- unlist( exbytx )
ex_ranks <- exbytx_unlist$exon_rank
names(ex_ranks) = 1:length(ex_ranks)
Idx_last_exon <-  tapply( ex_ranks , names(exbytx_unlist), function(x) names(x)[which.max(x)]) %>% as.numeric
Last_exons_50bp <- exbytx_unlist[Idx_last_exon] %>% resize(.,50,fix = "start")
Feature_matrix$Last_exons_50bp <- row_gr%over%Last_exons_50bp

cat("Feature 4:  Start 50bp of the last exon is attached.\n")

#Feature 5. CDS
Feature_matrix$CDS <- row_gr%over%cds(txdb)

cat("Feature 5:  CDS is attached.\n")

#Feature 6. 5'UTR.
Feature_matrix$UTR5 <- row_gr%over%fiveUTRsByTranscript(txdb)

cat("Feature 6:  5'UTR is attached.\n")

#Feature 7. 3'UTR.
Feature_matrix$UTR3 <- row_gr%over%threeUTRsByTranscript(txdb)

cat("Feature 7:  3'UTR is attached.\n")

#Feature 8. Long exon.
exs <- exons(txdb)
long_exs <- exs[width(exs) > 400]
Feature_matrix$long_exon <- row_gr%over%long_exs

cat("Feature 8:  Long exons (> 400nt) is attached.\n")

#Feature 9. Long 3'UTR.
UTR3 <- threeUTRsByTranscript(txdb)
long_UTR3 <- UTR3[sum(width(UTR3)) > 400]
Feature_matrix$long_UTR3 <- row_gr%over%long_UTR3

cat("Feature 9:  Long 3'UTR is attached.\n")

#Feature 10. Relative positioning on tx

Feature_matrix$Pos_Tx <- Relative_Pos_map(row_gr,exbytx,0)

cat("Feature 10:  Transcript relative position is attached.\n")

#Feature 11. Relative positioning UTR5.

Feature_matrix$Pos_UTR5 <- Relative_Pos_map(row_gr,fiveUTRsByTranscript(txdb),0)

cat("Feature 11:  5'UTR relative position is attached.\n")

#Feature 12. Relative positioning CDS.

Feature_matrix$Pos_CDS <- Relative_Pos_map(row_gr,cdsBy(txdb,by = "tx"),0)

cat("Feature 12:  CDS relative position is attached.\n")

#Feature 13. Relative positioning UTR3.

Feature_matrix$Pos_UTR3 <- Relative_Pos_map(row_gr,threeUTRsByTranscript(txdb),0)

cat("Feature 13:  3'UTR relative position is attached.\n")

#Feature 14. Relative positioning Exons.
exs_hg19 <- exons(txdb)
exs_grl <- split(exs_hg19,1:length(exs_hg19))
Feature_matrix$Pos_exons <- Relative_Pos_map(row_gr,exs_grl,0)

cat("Feature 14:  Exons relative position is attached.\n")

#Feature 15. Stop codon (301 bp center).
Stop_codons <- resize( cds(txdb) , 1, fix = "end" ) + 150

Feature_matrix$Stop_codons <- row_gr%over%Stop_codons

cat("Feature 14:  Stop codons 150bp flanking is attached.\n")

#Feature 16. Start codon (201 bp center).
Start_codons <- resize( cds(txdb) , 1, fix = "start" ) + 100

Feature_matrix$Start_codons <- row_gr%over%Start_codons

cat("Feature 16:  Start codons 100bp flanking is attached.\n")

#Feature 17. Struc hybridize.
Feature_matrix$Struc_hybridize <- row_gr%over%rBS2ndStructure::Struc_hg19

cat("Feature 17:  Hybridized structure region.\n")

#Feature 18. Struc loop.
Feature_matrix$Struc_loop  <- row_gr%over%Infer_loop(rBS2ndStructure::Struc_hg19)

cat("Feature 18:  Between paired hybridize structure region.\n")

#Feature 19. PC score 1nt.
Feature_matrix$PC_1nt <- scores(phastCons100way.UCSC.hg19,row_gr)$scores
Feature_matrix$PC_1nt[is.na(Feature_matrix$PC_1nt)] = mean(na.omit(Feature_matrix$PC_1nt))

cat("Feature 19:  Phast cons score of a single nucleotide.\n")

#Feature 20. PC score 201nt.
Feature_matrix$PC_201nt <- scores(phastCons100way.UCSC.hg19,row_gr+100)$scores
Feature_matrix$PC_201nt[is.na(Feature_matrix$PC_201nt)] = mean(na.omit(Feature_matrix$PC_201nt))

cat("Feature 20:  Phast cons score of 100 flanked nucleotide.\n")

#Feature 21. fitness consequences scores 1nt.
suppressWarnings( Feature_matrix$FC_1nt <- scores(fitCons.UCSC.hg19,row_gr)$scores )
Feature_matrix$FC_1nt[is.na(Feature_matrix$FC_1nt )] = mean(na.omit(Feature_matrix$FC_1nt))

cat("Feature 21:  fitness consequences scores of 1 nucleotide.\n")

#Feature 22. fitness consequences scores 5nt.
suppressWarnings( Feature_matrix$FC_5nt <- scores(fitCons.UCSC.hg19,row_gr + 2)$scores )
Feature_matrix$FC_5nt[is.na(Feature_matrix$FC_5nt )] = mean(na.omit(Feature_matrix$FC_5nt))

cat("Feature 22:  fitness consequences scores of 5 nucleotides (consensus motif).\n")

#Feature 23. lnc RNA.
coding_tx <- names(cdsBy(txdb, by = "tx"))
exbytx_nc <- exbytx[!names(exbytx)%in%coding_tx]
lnc_idx <- sum(width(exbytx_nc)) > 200
Feature_matrix$lncRNA  <- row_gr%over%exbytx_nc[lnc_idx]

cat("Feature 23:  lnc-RNA (non-coding RNA with length > 200bp).\n")

#Feature 24. snc RNA.
Feature_matrix$sncRNA  <- row_gr%over%exbytx_nc[!lnc_idx]

cat("Feature 24:  snc-RNA (small non-coding RNA with length <= 200bp).\n")

#Feature 25. Iso form numbers.
txbygenes <- transcriptsBy(txdb,by = "gene")
Feature_matrix$Isoform_num <- Properties_mapp(query_gr = row_gr,
                                              feature_gr = txbygenes,
                                              feature_properties = elementNROWS(txbygenes),
                                              no_map_val = 0)
cat("Feature 25:  Transcript Iso-form numbers.\n")

#Feature 26-37. Alternative motifs:
DSS <- DNAStringSet( Views(Hsapiens,row_gr + 2) )

is_motif <- function(motif,dss,exact = T) vcountPattern(DNAString(motif), dss, fixed = exact) > 0

#26. -AAACA
Feature_matrix$AAACA <- is_motif("AAACA",DSS,T)

#27. -GAACA
Feature_matrix$GAACA <- is_motif("GAACA",DSS,T)

#28. -AGACA
Feature_matrix$AGACA <- is_motif("AGACA",DSS,T)

#29. -GGACA
Feature_matrix$GGACA <- is_motif("GGACA",DSS,T)

#30. -AAACT
Feature_matrix$AAACT <- is_motif("AAACT",DSS,T)

#31. -GAACT
Feature_matrix$GAACT <- is_motif("GAACT",DSS,T)

#32. -AGACT
Feature_matrix$AGACT <- is_motif("AGACT",DSS,T)

#33. -GGACT
Feature_matrix$GGACT <- is_motif("GGACT",DSS,T)

#34. -AAACC
Feature_matrix$AAACC <- is_motif("AAACC",DSS,T)

#35. -GAACC
Feature_matrix$GAACC <- is_motif("GAACC",DSS,T)

#36. -AGACC
Feature_matrix$AGACC <- is_motif("AGACC",DSS,T)

#37. -GGACC
Feature_matrix$GGACC <- is_motif("GGACC",DSS,T)

cat("Feature 26-37:  Specific motifs are attached.\n")

#Feature 38-39. Batch features.
Feature_matrix$GC_cont_101bp <- as.numeric( letterFrequency( DNAStringSet( Views(Hsapiens,row_gr+50) ) , letters="CG", as.prob = T) )

cat("Feature 38:  Batch effect --- GC contents of 101 bp flanking region is attached.\n")

exbg_gr <- unlist(exbg_txdb)
exs_GC_freq <-  letterFrequency( DNAStringSet( Views(Hsapiens, exbg_gr) ) , letters="CG")
Genes_GC_freq <- tapply( exs_GC_freq, names(exbg_gr), sum )
Genes_length <- tapply( width(exbg_gr), names(exbg_gr), sum )
Genes_GC_cont = Genes_GC_freq/Genes_length

Feature_matrix$GC_cont_genes <- Properties_mapp(query_gr = row_gr,
                                                feature_gr = exbg_txdb,
                                                feature_properties = Genes_GC_cont,
                                                no_map_val = .5)

cat("Feature 39:  Batch effect --- GC contents of gene's exonic regions is attached.\n")

#Feature 40+ Features provided by user defined GRanges/GRanges List.
for(i in names(Feature_lst_gr)) {
  suppressWarnings( Feature_matrix[[i]] <- row_gr %over% Feature_lst_gr[[i]] )
  cat(paste0 ("Additional annotation feature: ", i , " is attached.\n"))
}

mcols(row_gr) = Feature_matrix
return(row_gr)




###Meta for logistic modeling

SE_Learning <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/J_GLM_Learning_2018_1_15/SE_Learning.rds")

SE <- SE_Learning[,colData(SE_Learning)$Role == "eraser"]

SE_ASSAY = assay(SE)
SE_ASSAY[SE_ASSAY == "0"] = NA
SE_ASSAY[SE_ASSAY == "-1"] = 0
assay(SE) = SE_ASSAY
colnames( SE ) = paste0("Eraser",1:ncol(SE_ASSAY))

SE <- predictors.annot(se = SE,
                       txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                       bsgnm = Hsapiens,
                       fc = fitCons.UCSC.hg19,
                       pc = phastCons100way.UCSC.hg19,
                       struct_hybridize = rBS2ndStructure::Struc_hg19,
                       feature_lst = Feature_List_hg19,
                       HK_genes_list = HK_hg19_eids)


SE_features_added <- predictors.annot(
  se = SE,
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  bsgnm = Hsapiens,
  fc = fitCons.UCSC.hg19,
  pc = phastCons100way.UCSC.hg19,
  struct_hybridize = rBS2ndStructure::Struc_hg19,
  feature_lst = Feature_lst_hg19,
  HK_genes_list = HK_hg19_eids
)


group_list_default = list(
  miRNA_RBP = Feature_lst_hg19
)

#Try to demonstrate predict evolutionary conservation with other features.

# Case study to demonstrate the usage of this tool, classify the differences between evolutionary conserved and evolutionary unconserved (by FC and PC) m6A sites reported by RMBase2.

Gtcoord_hg19 <- readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_hg19.rds")
Guitar::GuitarPlot(list(Start_codons),Gtcoord_hg19)
Guitar::GuitarPlot(list(Stop_codons),Gtcoord_hg19)
Guitar::GuitarPlot(list(TSS),Gtcoord_hg19)
Guitar::GuitarPlot(list(TSS[A_idx]),Gtcoord_hg19)
Guitar::GuitarPlot(list(Last_exons_400bp),Gtcoord_hg19)

Guitar::GuitarPlot(list(subsetByOverlaps( Last_exons_400bp, Stop_codons)),Gtcoord_hg19)
