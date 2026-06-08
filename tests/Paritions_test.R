library(testthat)

test_that(
  "Confirm all variations of partition mode work as expected",
  {
    # confirm no NAs occur when running CoExpress on default settings:
    df <- CoExpress(pbmc, GOIs)
    expect_false(anyNA(df$MOC_Ratio))
    # check if CoExpress runs when skip.extremes == T
    expect_shape(CoExpress(pbmc, GOIs, skip.extremes=T), nrow=15)
    # query input against HLA-DQB2 and expect HLA-DBQ1 to be the most co-expressed (biological truth):
    df2 <- CoExpress(pbmc, GOIs, gene2="HLA-DQB2")
    expect_true(df2$GeneA[order(df2$MOC_Ratio, decreasing=T)][1] == "HLA-DQB1")
    # confirm plotting does not cause errors:
    expect_no_error(ShowResiduals(df))
    expect_no_error(plotCoExpr(df))
  }
  
)


