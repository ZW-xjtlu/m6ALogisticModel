library(testthat)
test_that("meth k means cluster", {
 library(m6ALogisticModel)
 SE_CQN <- readRDS("SE_CQN_filtered.rds")
 Y <- plot_row_joint(SE_CQN,RETURN_INDX = T)
 expect_that(Y, is_a("integer"))
 expect_that(length(Y) == nrow(SE_CQN), is_true())
})

test_that("multinomial model", {
  glm_multinomial(Y,mcols(SE_CQN),"GLM_bg")
  })

test_that("go multinomial", {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  go_multinomial(Y,
                 row_ranges = rowRanges(SE_CQN),
                 txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                 orgDb = org.Hs.eg.db,
                 HDER = "GO_test",
                 GO_Slim = T
                 )
})
