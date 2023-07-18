library(Seurat)
library(irGSEA)

# seurat_obj: Seurat class of ST data
# gene_ls: list of gene sets (e.g. .rds files in gene_set directory)

seurat_obj <- irGSEA.score(object=seurat_obj,assay="RNA",slot="data",
                           custom=TRUE,geneset=gene_ls,msigdb=FALSE, 
                           method="AUCell", aucell.MaxRank=ceiling(0.2*nrow(seurat_obj)),
                           seeds=123,ncores=40)
