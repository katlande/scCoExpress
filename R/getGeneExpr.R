getGeneExpr <- function(SeuratObj, assay=NULL, slot="data", v=5){
  message("Getting background frequency information...")
  
  if(v >= 5){
    apply(SeuratObj@assays[[assay]][slot], 1, function(x){
      y <- table(as.numeric(x) > 0)
      if(length(y) == 2){
        return(y[["TRUE"]]/(y[["FALSE"]]+y[["TRUE"]]))
      } else if(! "TRUE" %in% names(y)){
        return(0)
      } else{
        return(1)
      }
    }) -> freqs
    gene_expr <- data.frame(gene=row.names(SeuratObj@assays[[assay]][slot]), occurance=freqs)
  } else {
    apply(SeuratObj@assays[[assay]][[slot]], 1, function(x){
      y <- table(as.numeric(x) > 0)
      if(length(y) == 2){
        return(y[["TRUE"]]/(y[["FALSE"]]+y[["TRUE"]]))
      } else if(! "TRUE" %in% names(y)){
        return(0)
      } else{
        return(1)
      }
    }) -> freqs
    gene_expr <- data.frame(gene=row.names(SeuratObj@assays[[assay]][[slot]]), occurance=freqs)
  }
  
  gene_expr <- gene_expr[order(gene_expr$occurance),]
  return(gene_expr)
}

