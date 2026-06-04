# Internal function
#' @noRd
calcMOC <- function(SeuratObj, a, s, ga, gb){
  if(!Seurat::DefaultAssay(SeuratObj)==a){
    Seurat::DefaultAssay(SeuratObj) <- a
  }
  vec_a <- as.numeric(unlist(Seurat::FetchData(object = SeuratObj, vars = ga, layer = s)))
  vec_b <- as.numeric(unlist(Seurat::FetchData(object = SeuratObj, vars = gb, layer = s)))
  
  sum(vec_a * vec_b) / sqrt(sum(vec_a^2)*sum(vec_b^2)) -> m
  return(m)
}
