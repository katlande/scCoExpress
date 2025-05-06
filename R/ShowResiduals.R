ShowResiduals <- function(CoExpr){
  
  reg <- lm(CoExpr$MOC_Z~CoExpr$MOC_Ratio)
  message(paste("Ratio-Z correlation has an Rsq of:", formatC(summary(reg)$r.squared, digits = 4)))
  
  CoExpr <- CoExpr[!is.na(CoExpr$MOC_Z),]
  ggplot2::ggplot(CoExpr)+
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey", linewidth = 0.25),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey", linewidth = 0.25),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black", face = "italic"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      legend.title = ggplot2::element_text(hjust = 0.5),
      strip.background = ggplot2::element_rect(fill="black"),
      strip.text = ggplot2::element_text(colour="white", face="bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.key = ggplot2::element_rect(fill = "white"))+
    ggplot2::geom_smooth(mapping=ggplot2::aes(x=MOC_Ratio, y=MOC_Z),method="lm", linetype="dashed", 
                         color="black", linewidth=0.25, formula = y ~ x)+
    ggplot2::ylab("Z-Score")+
    ggplot2::xlab("MOC Ratio") -> g
  
  if(any(!is.na(CoExpr$Zadj))){
    g+ggplot2::annotate(geom="segment", color="purple",
                        x=CoExpr$MOC_Ratio, xend=CoExpr$MOC_Ratio, 
                        y=CoExpr$Zadj, yend=CoExpr$MOC_Z)+
      ggplot2::geom_point(mapping=ggplot2::aes(x=MOC_Ratio, y=MOC_Z), color="darkblue", alpha=0.8)+
      ggplot2::geom_point(mapping=ggplot2::aes(x=MOC_Ratio, y=Zadj), color="darkred", alpha=0.8) -> g
  } else {
    g+ ggplot2::geom_point(mapping=ggplot2::aes(x=MOC_Ratio, y=MOC_Z), color="darkblue", alpha=0.8) -> g
  }
  
  return(g) 
}
