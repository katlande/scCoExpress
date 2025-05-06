plotCoExpr <- function(co){
  
  co <- co[!is.na(co$MOC_Z),]
  
  if(any(is.na(co$Zadj))){
    message("Adjusted Z-scores contain NAs! Using unadjusted Z-scores for plotting.")
    co <- setNames(co[colnames(co) %in% c("GeneA", "GeneB", "MOC_Z")], c("GeneA", "GeneB", "Z"))
  } else {
    message("Using adjusted Z-scores for plotting.")
    co <-setNames( co[colnames(co) %in% c("GeneA", "GeneB", "Zadj")], c("GeneA", "GeneB", "Z"))
  }
  
  df <- rbind(co, setNames(co[c(2,1,3)], c("GeneA", "GeneB", "Z")))
  df <- unique(df[order(df$Z, decreasing = T),])
  # make a matrix and fill in missing values:
  suppressMessages(reshape2::dcast(df, GeneA~GeneB, drop = F)) -> mat
  row.names(mat) <- mat$GeneA
  mat <- mat[!colnames(mat) == "GeneA"]
  mat[is.na(mat)] <- 0
  v <- max(mat)*5
  
  for(i in 1:nrow(mat)){
    mat[rownames(mat)[i],rownames(mat)[i]]<- v
  }
  
  # cluster the matrix by correlation:
  tree_row <- hclust(as.dist(1 - cor(t(mat))), method = "ward.D2") 
  row_order <- row.names(mat)[tree_row$order]
  df$GeneA <- factor(df$GeneA, levels=row_order)
  df$GeneB <- factor(df$GeneB, levels=row_order)
  
  # plot:
  ggplot2::ggplot(df, ggplot2::aes(x=GeneA, y=GeneB, fill=Z))+
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradient2("Co-Localization\nScore", midpoint=0, high="red", low="blue", mid="white", na.value = "grey")+
    ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill="black"),
                   legend.title = ggplot2::element_text(size=8, hjust=0.5),
                   plot.title = ggplot2::element_text(hjust=0.5, face = "bold"),
                   axis.text.x = ggplot2::element_text(colour = "black", angle=45, vjust=1, hjust=1))+
    ggplot2::geom_hline(yintercept = c(0:length(unique(df$GeneA)))+0.5, linewidth=0.5)+
    ggplot2::geom_vline(xintercept = c(0:length(unique(df$GeneA)))+0.5, linewidth=0.5)+
    ggplot2::scale_x_discrete("", expand=c(0,0))+
    ggplot2::scale_y_discrete("", expand=c(0,0)) -> outG
  
  return(outG)
}