library(dplyr)
library(Seurat)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/RNANormalization/RNANormalization.R')

dataset_used=readRDS('../dataset_used.rds')
l=sapply(dataset_used,function(dataset){
    seurat_obj=readRDS(paste0('../1_data_preprocessing/',dataset,'/',dataset,'.rds'))
    ind=intersect(rownames(LCP_cluster_result),colnames(seurat_obj))
    if (length(ind)<=50){
        return(NULL)
    }
    seurat_obj=seurat_obj[,ind]
    exp_mat=apply(seurat_obj@assays$RNA@counts,1,sum)
    return(exp_mat)
},USE.NAMES=1)

# l=l[!is.null(l)]
l=l[lapply(l,length)!=0]

df=dplyr::bind_rows(l)
df=data.frame(df)
df[is.na(df)]=0
rownames(df)=names(l)

# count to TPM
expr_norm=NormalizeCount(df,length.type='transcript',method='tpm')
# predict TCGA immune subtype
semble(X=expr_norm,geneids='symbol')