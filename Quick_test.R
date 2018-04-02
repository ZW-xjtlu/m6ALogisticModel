#Quick test for eval row joint.
library(SummarizedExperiment)
SE_CQN <- readRDS("/Users/zhenwei/Documents/GitHub/mRNA-cor/results_tidy/m6A_DESEQ2_cqn.rds")
assays(SE_CQN)$m6Alog2FC[assays(SE_CQN)$geneExpression < 8.5] = NA
indx_keep <- rowSums( is.na(assay(SE_CQN)) ) <= 10
SE_CQN <- SE_CQN[indx_keep,]
assays(SE_CQN)$m6Alog2FC = scale(assays(SE_CQN)$m6Alog2FC)

#Annotate features
library(m6ALogisticModel)

#Eval_row_joint
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

Feature_List_expanded_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

SE_features_added <- predictors.annot(se = SE_CQN,
                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      bsgnm = Hsapiens,
                                      fc = fitCons.UCSC.hg19,
                                      pc = phastCons100way.UCSC.hg19,
                                      struct_hybridize = Struc_hg19,
                                      feature_lst = Feature_List_expanded_hg19,
                                      HK_genes_list = HK_hg19_eids)


SE = SE_features_added
row_Mads <- rowMads(assay(SE),na.rm = T)
SE = SE[row_Mads > .8,]

Eval_row_joint(SE,"Row_Joint_test",K = 2)
