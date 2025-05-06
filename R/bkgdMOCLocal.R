bkgdMOCLocal <- function(df, gene_expr, SeuratObj, bkgd=5, assay=NULL, slot="data"){
  bkgd <- ceiling(bkgd/3)
  
  message("Checking backgrounds MOCs for each gene combination...")
  mocs <- list()
  for(i in 1:nrow(df)){
    
    # find the most similarly expressed background genes for colocalization targets:
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
    bkgd_a <- sample(bkgd_a)
    
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
    bkgd_b <- sample(bkgd_b)
    
    # Re-sample the orders 3 times:
    resample_vec_moc <- c()
    for(q in 1:3){
      # randomly re-order the vector until you have no identical gene comparisons, if that does happen (rare):
      j <- 0
      while(any(apply(data.frame(a=bkgd_a, b=bkgd_b), 1, function(x){identical(as.character(x)[1], as.character(x)[2])}))){
        bkgd_a <- sample(bkgd_a)
        j <- j+1
        if(j>50){
          warning(paste0("Background for a & b genes is highly similar! It is taking a long time to find unique background gene pairs for ", a, " and ", b, "..."))
          j <- 0
        }
      }
      
      # gene the mean MOC for all the background comparisons:
      k <- as.numeric(unlist(lapply(1:length(bkgd_b), function(x){
        calcMOC(SeuratObj, assay, slot, bkgd_a[as.numeric(x)], bkgd_b[as.numeric(x)])
      }) ))
      
      resample_vec_moc <- c(resample_vec_moc, k)
      bkgd_b <- sample(bkgd_b)
      bkgd_a <- sample(bkgd_a)
    }
    
    mocs <- append(mocs, list(resample_vec_moc))
    
  }
  
  return(mocs)
}
