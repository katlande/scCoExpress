#' @export
CoExpress <- function(obj, target_genes, gene2=NULL, seuratAssay=NULL, seuratSlot="data", 
                      nPartitions=7, nPermutations=50, BkgdGeneExpr=NULL, topExcl=0.98, bottomExcl=0.005, 
                      local.perms=6, skip.extremes=F, seurat=5, show.progress=FALSE){

 if(is.null(seuratAssay)){
    seuratAssay <- Seurat::DefaultAssay(obj)
  }
  
  # Background gene expression:
  if(is.null(BkgdGeneExpr)){
    BkgdGeneExpr <- getGeneExpr(obj, assay=seuratAssay, slot=seuratSlot, v=seurat)
  } 
  
  if(is.null(gene2)){
    
    # remove non-expressed genes:
    if(any(BkgdGeneExpr$occurance[BkgdGeneExpr$gene %in%  target_genes] == 0)){
      message(paste("Some genes have no expression, removing:", toString(BkgdGeneExpr$gene[BkgdGeneExpr$gene %in%  target_genes &  BkgdGeneExpr$occurance == 0]) ))
      target_genes <- target_genes[!target_genes %in% BkgdGeneExpr$gene[BkgdGeneExpr$gene %in%  target_genes &  BkgdGeneExpr$occurance == 0]]
    }
    
    mat <- utils::combn(unique(target_genes), 2)
    gene_A <- mat[1,]
    gene_B <- mat[2,]
    message(paste0(length(unique(target_genes)), " unique input genes becomes ", length(gene_B), " pairwise comparisons."))
    
  } else{
    if(length(gene2) > 1){
      warning("gene2 must be NULL or a single gene name; only the first gene in gene2 will be used.")
      gene2 <- gene2[1]
    }
    
    # remove non-expressed genes:
    if(BkgdGeneExpr$occurance[BkgdGeneExpr$gene == gene2] == 0){
      stop("Error: Gene2 has no expression, so it cannot be co-expressed!")
    }
    
    if(any(BkgdGeneExpr$occurance[BkgdGeneExpr$gene %in%  target_genes] == 0)){
      message(paste("Some genes have no expression, removing:", toString(BkgdGeneExpr$gene[BkgdGeneExpr$gene %in%  target_genes &  BkgdGeneExpr$occurance == 0]) ))
      target_genes <- target_genes[!target_genes %in% BkgdGeneExpr$gene[BkgdGeneExpr$gene %in%  target_genes &  BkgdGeneExpr$occurance == 0]]
    }
    
    gene_A <- target_genes[!target_genes %in% gene2]
    gene_B <- gene2
  }

  message("Calculating Co-Localization in Pairs...")
  lapply(1:length(gene_A), function(i){
    
    if(length(gene_B)==1){
      tmp <- suppressMessages(GetCoExpr(SeuratObj=obj, geneA=gene_A[i], geneB=gene_B, assay=seuratAssay, slot=seuratSlot)) 
    } else{
      tmp <- suppressMessages(GetCoExpr(SeuratObj=obj, geneA=gene_A[i], geneB=gene_B[i], assay=seuratAssay, slot=seuratSlot)) 
    }
    return(tmp)
  }) -> GetCoExpr_OUTPUTFILE
  do.call(rbind, GetCoExpr_OUTPUTFILE) -> GetCoExpr_OUTPUTFILE
  
  message("Getting MOC/Background Ratios...")
  if(is.null(nPartitions)){
    message("Using local background genes instead of partitions...")
    
    
    MOC_backgrounds <- bkgdMOCLocal(df=GetCoExpr_OUTPUTFILE, gene_expr=BkgdGeneExpr, SeuratObj=obj, bkgd=nPermutations,
                                    assay=seuratAssay, slot=seuratSlot)
    
    GetCoExpr_OUTPUTFILE <- cbind(GetCoExpr_OUTPUTFILE, data.frame(MOC_bkgd=unlist(lapply(MOC_backgrounds, mean))))
    
    GetCoExpr_OUTPUTFILE$MOC_Z <- NA
    for(i in 1:nrow(GetCoExpr_OUTPUTFILE)){
      GetCoExpr_OUTPUTFILE$MOC_Z[i] <- (GetCoExpr_OUTPUTFILE$MOC[i]-mean(MOC_backgrounds[[i]]))/stats::sd(MOC_backgrounds[[i]])
    }
    
    GetCoExpr_OUTPUTFILE$MOC_Ratio <- GetCoExpr_OUTPUTFILE$MOC/GetCoExpr_OUTPUTFILE$MOC_bkgd
    
  } else {
    
    MOC_backgrounds <- bkgdMOCParition(obj, gene_expr=BkgdGeneExpr, 
                                       permutations=nPermutations, 
                                       paritions=nPartitions, 
                                       assay=seuratAssay, 
                                       slot=seuratSlot, 
                                       top=topExcl, 
                                       bottom=bottomExcl)
    
    BkgdGeneExpr <- as.data.frame(MOC_backgrounds[[1]])
    expr_groups <- MOC_backgrounds[[2]]
    
    
    # very lowly and very highly expressed genes should be queried locally:
    if(any(BkgdGeneExpr$group[BkgdGeneExpr$gene %in% unique(c(gene_A, gene_B))] %in% c("HIGHLY EXPRESSED", "LOWLY EXPRESSED"))){
      
      message(paste0(as.numeric(table(BkgdGeneExpr$group[BkgdGeneExpr$gene %in% unique(c(gene_A, gene_B))] %in% c("HIGHLY EXPRESSED", "LOWLY EXPRESSED"))["TRUE"]), " Genes with extreme expression detected..."))
      if(skip.extremes == T){
        message("Comparisons containing these genes will be skipped...")
      } else{
        message("Will use local background for these comparisons...")
      }
    } 
    
    GetCoExpr_OUTPUTFILE$MOC_bkgd <- NA
    GetCoExpr_OUTPUTFILE$MOC_Z <- NA
    GetCoExpr_OUTPUTFILE$Comparison_Type <- NA
    for(i in 1:nrow(GetCoExpr_OUTPUTFILE)){
      
      a <- GetCoExpr_OUTPUTFILE$GeneA[i]
      b <- GetCoExpr_OUTPUTFILE$GeneB[i]
      GetCoExpr_OUTPUTFILE$Comparison_Type[i] <- paste0(BkgdGeneExpr$group[BkgdGeneExpr$gene==a], "-", BkgdGeneExpr$group[BkgdGeneExpr$gene==b])
      
      # if the comparison contains an edge gene:
      if(! BkgdGeneExpr$group[BkgdGeneExpr$gene==a] %in% c("HIGHLY EXPRESSED", "LOWLY EXPRESSED") & 
         ! BkgdGeneExpr$group[BkgdGeneExpr$gene==b] %in% c("HIGHLY EXPRESSED", "LOWLY EXPRESSED")){
        if(show.progress){
          cat(".")
        }
        g <- paste0(BkgdGeneExpr$group[BkgdGeneExpr$gene==a], BkgdGeneExpr$group[BkgdGeneExpr$gene==b])
        GetCoExpr_OUTPUTFILE$MOC_bkgd[i] <- mean(expr_groups[[g]])
        GetCoExpr_OUTPUTFILE$MOC_Z[i] <- (GetCoExpr_OUTPUTFILE$MOC[i]-mean(expr_groups[[g]]))/stats::sd(expr_groups[[g]])
      } else{
        
        
        if(skip.extremes){
          GetCoExpr_OUTPUTFILE$MOC_bkgd[i] <- NA
          GetCoExpr_OUTPUTFILE$MOC_Z[i] <- NA
        }else{
          if(show.progress){
            cat("*")
          }
          # get the background and Z score for local comparisons here...
          v <- Local_Once(SOBJ = obj, ass = seuratAssay, sl = seuratSlot, a = GetCoExpr_OUTPUTFILE$GeneA[i], b = GetCoExpr_OUTPUTFILE$GeneB[i], bkgd = local.perms, gene_expr = BkgdGeneExpr)
          GetCoExpr_OUTPUTFILE$MOC_bkgd[i] <- mean(v)
          GetCoExpr_OUTPUTFILE$MOC_Z[i] <- (GetCoExpr_OUTPUTFILE$MOC[i]-mean(v))/stats::sd(v)
          rm(v)
        }
      }
      
    }
    GetCoExpr_OUTPUTFILE$MOC_Ratio <- GetCoExpr_OUTPUTFILE$MOC/GetCoExpr_OUTPUTFILE$MOC_bkgd
  }
  
  
  if(nrow(GetCoExpr_OUTPUTFILE[!is.na(GetCoExpr_OUTPUTFILE$MOC_Z),]) > 5){
    
    if(nrow(GetCoExpr_OUTPUTFILE[!is.na(GetCoExpr_OUTPUTFILE$MOC_Z),]) < 20){
      warning("A low number of comparisons may cause inaccurate Z-score adjustment.")
    }
    
    tmp <- GetCoExpr_OUTPUTFILE[!is.na(GetCoExpr_OUTPUTFILE$MOC_Z),]
    message("Adjusting Z-scores...")
    reg <- stats::lm(tmp$MOC_Z~tmp$MOC_Ratio)
    rsq <- summary(reg)$r.squared
    
    if(rsq > 0.7){
      unlist(lapply(GetCoExpr_OUTPUTFILE$MOC_Ratio, function(x){
        lin(x, reg)
      })) -> GetCoExpr_OUTPUTFILE$Zadj
      if(rsq < 0.9){
        warning(paste0("R-squared value is below 0.9; Z-score adjustment may be inaccurate!"))
        if(is.null(nPartitions)){
          message("Consider using a greater number of permutations to increase Z-score adjustment accuracy.")
        } else {
          message("Consider using a greater number of permutations and/or partitions to increase Z-score adjustment accuracy.")
        }
      } 
    } else {
      
      if(is.null(nPartitions)){
        warning(paste0("R-squared too low to reliably adjust Z-score (", rsq, ")... Skipping.\nConsider using a greater number of permutations. Large Z-score residuals may indicate unreliable results!"))
      } else {
        warning(paste0("R-squared too low to reliably adjust Z-score (", rsq, ")... Skipping.\nConsider using a greater number of permutations or partitions. Large Z-score residuals may indicate unreliable results!"))
      }
      
    }
    
  } else{
    message("Too few comparisons for Z-score adjustment")
    GetCoExpr_OUTPUTFILE$Zadj <- NA
  }
  
  
  # calculate significance:
  if(any(colnames(GetCoExpr_OUTPUTFILE) == "Zadj")){
    GetCoExpr_OUTPUTFILE$p <- (1 - stats::pnorm(GetCoExpr_OUTPUTFILE$Zadj))
  } else {
    message("Calculating significance using unadjusted Z-score...")
    GetCoExpr_OUTPUTFILE$p <- (1 - stats::pnorm(GetCoExpr_OUTPUTFILE$MOC_Z))
  }
  
  GetCoExpr_OUTPUTFILE$FDR <- stats::p.adjust(GetCoExpr_OUTPUTFILE$p, method = "fdr")
  
  return(GetCoExpr_OUTPUTFILE)
}
