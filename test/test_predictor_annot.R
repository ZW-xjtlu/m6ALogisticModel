library(testthat)
test_that("predictors annotation", {
  library(m6ALogisticModel)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(fitCons.UCSC.hg19)
  library(phastCons100way.UCSC.hg19)

  Additional_features_hg19 = list(
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

  SE_CQN <- readRDS("SE_CQN_filtered.rds")

  SE_features_added <- predictors_annot(se = SE_CQN,
                                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        bsgnm = Hsapiens,
                                        fc = fitCons.UCSC.hg19,
                                        pc = phastCons100way.UCSC.hg19,
                                        struct_hybridize = Struc_hg19,
                                        feature_lst = Additional_features_hg19,
                                        hk_genes_list = HK_hg19_eids,
                                        standardization = FALSE)

  expect_that(SE_features_added, is_a("SummarizedExperiment"))
})

