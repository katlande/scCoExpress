GetCoExpr <- function(SeuratObj, geneA, geneB, assay=NULL, slot="data"){
  
  # Manders' Overlap Coefficient (MOC)
  # sum the product of the intensities of each pixel in two channels (e.g., red and green), then divide that sum by the square root of the product of the sums of the squares of the intensities in each channel
  MOC <-  calcMOC(SeuratObj, assay, slot, geneA, geneB)
  return(data.frame(GeneA=geneA, 
                    GeneB=geneB,
                    MOC=MOC))
}
