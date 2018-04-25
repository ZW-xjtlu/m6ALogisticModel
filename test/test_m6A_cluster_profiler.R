library(testthat)
test_that("meth k means cluster", {
 library(m6ALogisticModel)
 SE_CQN <- readRDS("SE_CQN_filtered.rds")[1:500,]
 Y <- plot_row_joint(SE_CQN,RETURN_INDX = T)
 expect_that(Y, is_a("integer"))
 expect_that(length(Y) == nrow(SE_CQN), is_true())
})

test_that("multinomial model", {
  glm_multinomial(Y,mcols(SE_CQN),"test","x")
  })

