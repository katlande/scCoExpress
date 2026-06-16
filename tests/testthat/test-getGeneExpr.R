library(testthat)

test_that(
  "Confirm gene background expression function works as expected",
  {
    # confirm background gene expression df is generated:
    bk <- getGeneExpr(pbmc)
    expect_shape(bk, nrow = 13714)
    # confirm quantile plotting works:
    expect_no_error(CheckMaxQuantile(bk))
  }
  
)

