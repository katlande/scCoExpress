library(testthat)

test_that(
  "Confirm all variations of local mode work as expected",
  {
    # confirm no NAs occur when running CoExpress on default local settings:
    dl <- CoExpress(pbmc, GOIs, nPartitions = NULL)
    expect_false(anyNA(dl$MOC_Ratio))
    expect_shape(dl, nrow=15)
    # query input against HLA-DQB2 and expect HLA-DBQ1 to be the most co-expressed (biological truth):
    dl2 <- CoExpress(pbmc, GOIs, gene2="HLA-DQB2", nPartitions = NULL)
    expect_true(dl2$GeneA[order(dl2$MOC_Ratio, decreasing=T)][1] == "HLA-DQB1")
  }
  
)


