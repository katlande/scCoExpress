# Internal function
#' @noRd
bkgdMOCLocal <- function(df, gene_expr, SeuratObj, bkgd=5, assay=NULL, slot="data"){
  bkgd <- ceiling(bkgd/3)
  
  message("Checking backgrounds MOCs for each gene combination...")
  tryCatch(lapply(1:nrow(df), function(i){
    # get background genes with similar expression to gene A:
    a <- df$GeneA[i]
    aV <- which(gene_expr$gene == a)
    if(aV-bkgd <= 0){
      start_a <- aV+1
      end_a <- aV+(bkgd*2)
      bkgd_a <- gene_expr$gene[start_a:end_a]
    } else if(aV+bkgd > nrow(gene_expr)){
      start_a <- aV-(bkgd*2)
      end_a <- aV-1
      bkgd_a <- gene_expr$gene[start_a:end_a]
    } else{
      bkgd_a <- gene_expr$gene[c((aV-bkgd):(aV-1), (aV+1):(aV+bkgd))]
    }
    # bkgd_a <- sample(bkgd_a)
    
    # get background genes with similar expression to gene B:
    b <- df$GeneB[i]
    bV <- which(gene_expr$gene == b)
    if(bV-bkgd <= 0){
      start_b <- bV+1
      end_b <- bV+(bkgd*2)
      bkgd_b <- gene_expr$gene[start_b:end_b]
    } else if(bV+bkgd > nrow(gene_expr)){
      start_b <- bV-(bkgd*2)
      end_b <- bV-1
      bkgd_b <- gene_expr$gene[start_b:end_b]
    } else{
      bkgd_b <- gene_expr$gene[c((bV-bkgd):(bV-1), (bV+1):(bV+bkgd))]
    }
    # bkgd_b <- sample(bkgd_b)
    
    lapply(1:3, function(q){
      # re-order
      bkgd_b <- sample(bkgd_b)
      bkgd_a <- sample(bkgd_a)
      
      # prevent self-to-self comparisons:
      j <- 0
      while(any(apply(data.frame(a=bkgd_a, b=bkgd_b), 1, function(x){identical(as.character(x)[1], as.character(x)[2])}))){
        bkgd_a <- sample(bkgd_a)
        j <- j+1
        if(j>50){
          warning(paste0("Cannot find enough unique background gene pairs for ", a, " and ", b, "!"))
        }
        if(j > 100){
          stop(paste("Could not find enough valid gene pairs to create a null background! There are likely not enough genes in this dataset to find", permutations, "permutations of null."))
        }
      }
      
      # get the mean MOC for all the background comparisons:
      k <- as.numeric(unlist(lapply(1:length(bkgd_b), function(x){
        calcMOC(SeuratObj, assay, slot, bkgd_a[as.numeric(x)], bkgd_b[as.numeric(x)])
      }) ))
      
      return(k)
    }) -> resample_vec_moc
    return(unlist(resample_vec_moc))
  })) -> mocs
  
  return(mocs)
}
