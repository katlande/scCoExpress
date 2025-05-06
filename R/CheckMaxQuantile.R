CheckMaxQuantile <- function(exprDF, quantiles=c(0.9, 0.96, 0.98)){
  
  ggplot2::ggplot(exprDF, ggplot2::aes(x=occurance))+
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
    
    ggplot2::annotate(geom="rect", xmin=-Inf, xmax=as.numeric(quantile(exprDF$occurance, quantiles[1])), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)+
    ggplot2::annotate(geom="rect", xmin=-Inf, xmax=as.numeric(quantile(exprDF$occurance, quantiles[2])), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)+
    ggplot2::annotate(geom="rect", xmin=-Inf, xmax=as.numeric(quantile(exprDF$occurance, quantiles[3])), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)+
    ggplot2::annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2)+
    ggplot2::geom_density(fill="lightblue")+
    ggplot2::geom_vline(xintercept = as.numeric(quantile(exprDF$occurance, quantiles[1])), linetype="dashed", linewidth=0.5, color="navy")+
    ggplot2::geom_vline(xintercept = as.numeric(quantile(exprDF$occurance, quantiles[2])), linetype="dashed", linewidth=0.5, color="navy")+
    ggplot2::geom_vline(xintercept = as.numeric(quantile(exprDF$occurance, quantiles[3])), linetype="dashed", linewidth=0.5, color="navy")+
    ggplot2::scale_x_continuous("Fraction of Cells Expressing"#,
                                #breaks=c(0, log10(1.1), log10(1.25), log10(1.5), log10(1.75), log10(2)),
                                #labels=c(0, 0.1, 0.25, 0.5, 0.75, 1), limits=c(0, log10(2))
    ) -> g
  
  ggplot2::ggplot_build(g)$layout$panel_scales_y[[1]]$range$range[[2]] -> ymax
  
  g+ggplot2::annotate(geom="text", hjust=1.2, x=as.numeric(quantile(exprDF$occurance, quantiles[1])), y=ymax*0.6, color="navy", label=quantiles[1], size=5)+
    ggplot2::annotate(geom="text", hjust=1.2, x=as.numeric(quantile(exprDF$occurance, quantiles[2])), y=ymax*0.4, color="navy", label=quantiles[2], size=5)+
    ggplot2::annotate(geom="text", hjust=1.2, x=as.numeric(quantile(exprDF$occurance, quantiles[3])), ymax*0.6, color="navy", label=quantiles[3], size=5)+
    ggplot2::scale_y_continuous("Gene Density", expand=c(0,0), limits=c((0-(ymax*0.0001)), ymax*1.1)) -> g
  
  return(g)
}