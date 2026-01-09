bkgdMOCParition <- function(SeuratObj, gene_expr, permutations=5, paritions=4, assay=NULL, slot="data", top, bottom){
  
  # high and low expressiong cut-offs:
  minThresh <- nrow(gene_expr[gene_expr$occurance <= bottom,])
  maxThresh <- nrow(gene_expr[gene_expr$occurance < quantile(gene_expr$occurance, top),])

  # updated to ignore bug created by too large 0 pool:
  gene_expr$group <- NA
  locs <- seq(minThresh, maxThresh+1, by=(as.integer(length(minThresh:maxThresh)/7)))
  locs <- c(1,locs,nrow(gene_expr))
  
  j <- 1
  for(i in 2:length(locs)){
    
    if(i == 2){
      gene_expr$group[locs[(i-1)]:locs[i]] <- "LOWLY EXPRESSED"
    } else if(i == length(locs)){
      gene_expr$group[locs[(i-1)]:locs[i]] <- "HIGHLY EXPRESSED"
    } else {
      gene_expr$group[locs[(i-1)]:locs[i]] <- paste0("g",(j-1))
    }
    
    
    j<-j+1
  }
  
  comparisons <- cbind(combn(unique(gene_expr$group[!gene_expr$group %in% c("LOWLY EXPRESSED", "HIGHLY EXPRESSED")]), 2), t(matrix(c(unique(gene_expr$group[!gene_expr$group %in% c("LOWLY EXPRESSED", "HIGHLY EXPRESSED")]), unique(gene_expr$group[!gene_expr$group %in% c("LOWLY EXPRESSED", "HIGHLY EXPRESSED")])),nrow=paritions)))
  compA <- comparisons[1,]
  compB <- comparisons[2,]
  
  message(paste0(paritions, " partitions become ", length(compA), " possible combinations..."))
  message(paste0("Getting partition combination backgrounds with ", permutations, " permutations..."))
  
  mocs <- list()
  for(c in 1:length(compA)){
    cat(".")
    bkgd_a <- gene_expr$gene[gene_expr$group==compA[c]]
    bkgd_b <- gene_expr$gene[gene_expr$group==compB[c]]
    
    bkgd_b_tmp <- sample(bkgd_b, replace = T, size = permutations)
    bkgd_a_tmp <- sample(bkgd_a, replace = T, size = permutations)
    
    # randomly re-order the vector until you have no identical gene comparisons, if that does happen (rare):
    j <- 0
    while(any(apply(data.frame(a=bkgd_a_tmp, b=bkgd_b_tmp), 1, function(x){identical(as.character(x)[1], as.character(x)[2])}))){
      bkgd_a_tmp <- sample(bkgd_a_tmp)
      j <- j+1
      if(j>50){
        warning(paste0("Background for a & b genes is highly similar! It is taking a long time to find unique background gene pairs for ", a, " and ", b, "..."))
        j <- 0
      }
    }
    
    k <- as.numeric(unlist(lapply(1:length(bkgd_b_tmp), function(x){
      calcMOC(SeuratObj, assay, slot, bkgd_a_tmp[as.numeric(x)], bkgd_b_tmp[as.numeric(x)])
    }) ))
    
    
    if(compA[c] == compB[c]){
      mocs <- append(mocs, list(k))
      names(mocs)[length(mocs)] <- c(paste0(compA[c],compB[c]))
    } else {
      mocs <- append(mocs, list(k))
      mocs <- append(mocs, list(k))
      names(mocs)[(length(mocs)-1):length(mocs)] <- c(paste0(compA[c],compB[c]), paste0(compB[c],compA[c]))
    }
    
    
  }
  cat("\n")
  return(list(list(gene_expr), mocs))
}
