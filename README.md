## Single Cell Co-Expression
These functions identify the co-expression intensity of pairs of genes within individual cells or nuclei, using single cell or spatial data pre-processed in a Seurat Object. These functions use a modified version of Mander's Overlap Coefficient between any pair of genes, wherein an average background MOC is calculated using genes with a similar abundance to the target gene pair. All results are reported as a comparison of the target MOC to the local background MOC to account for the sparity of single cell data. 

#### The Single-Cell MOC Calculation
* c = The cells in the experiment
* A = The vector of expression of gene A in all cells k
* B = The vector of expression of gene B in all cells k

$$ \left( \sum_{k=1}^c A_k B_k \right) / \sqrt{ \left( \sum_{k=1}^c A_k^2 \right)  * \left( \sum_{k=1}^c B_k^2 \right) } $$

**All modes of scCoExpress** return a ratio of the MOC of the genes of interest against an average or distribution of MOCs from pairs of genes with similar expression to the target genes. 
**Partition mode** splits background genes into groups based on expression similarity. For any pair of two target genes, the background MOCs are generated from the partitions with the most similar expression to the targets. When querying a large number of genes, this will be fastest.
**Local mode** generates backgrounds using the genes that have the most similar expression to the two target genes. This mode is expected to be the most accurate, but for a large number of genes, this mode takes a *very* long time. However, if comparing a number of genes *smaller* than the number of partitions in partition mode, local mode will actually be faster. 


## Installation:
remotes::install_github("katlande/scCoExpress")

## Vignette:
[Instructions for Use](https://github.com/katlande/scCoExpress/blob/main/scCoExpress.md)
