# deconvolution using spacexr (RCTD)

library(dplyr)
library(Seurat)
library(spacexr)

cancer_type=gsub('_[0-9]+','',dataset_1)
stRNA_sce=readRDS(file.path('../1_data_preprocessing',file.path(dataset_1,paste0(dataset_1,'.rds'))))
scRNA_sce=readRDS(file.path('scRNA_used',paste0(cancer_type,'.rds')))
scRNA_sce=scRNA_sce[,Idents(scRNA_sce)!='Undefined']
RCTD_result <- create.RCTD(SpatialRNA(data.frame(x=rep(1,dim(stRNA_sce)),y=rep(1,dim(stRNA_sce))),GetAssayData(stRNA_sce,assay='RNA',slot='counts'),use_fake_coords=TRUE), 
        Reference(GetAssayData(scRNA_sce,assay='RNA',slot='counts'),Idents(scRNA_sce)), 
        max_cores=1,test_mode=FALSE,CELL_MIN_INSTANCE=15,UMI_min=0) %>%
    run.RCTD(doublet_mode = 'full')
df=data.frame(RCTD_result@results$weights,check.names=FALSE)
df=normalize_weights(df)