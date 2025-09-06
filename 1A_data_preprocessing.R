library(Seurat)
library(dplyr)
library(data.table)

# sample filtration
qc_table=parallel::mclapply(dataset,function(dataset_1){
    
    print(dataset_1)
    
    ### read in data to Seurat
    data.dir=file.path('../0_public',dataset_1)
    
    if (file.exists(file.path(data.dir,'spatial','tissue_lowres_image.png'))){
        img=Read10X_Image(file.path(data.dir,'spatial'),image.name = "tissue_lowres_image.png")
        res_type='lowres'
    }else{
        img=Read10X_Image(file.path(data.dir,'spatial'),image.name = "tissue_hires_image.png")
        res_type='hires'
    }
    Seurat::DefaultAssay(object = img)='RNA'
    if (file.exists(file.path(data.dir,'filtered_feature_bc_matrix.h5'))){
        data=Read10X_h5(file.path(data.dir,'filtered_feature_bc_matrix.h5'))
    } else if (dir.exists(file.path(data.dir,'filtered_feature_bc_matrix'))){
        df=fread(file.path(data.dir,'filtered_feature_bc_matrix','features.tsv.gz'))
        if (ncol(df)==1) {data=Read10X(file.path(data.dir,'filtered_feature_bc_matrix'),gene.column=1)}
        else{data=Read10X(file.path(data.dir,'filtered_feature_bc_matrix'))}
    }
    
    barcode_filtered=intersect(colnames(data),rownames(img@coordinates))
    data=data[,barcode_filtered]
    
    seurat_obj <- Seurat::CreateSeuratObject(counts=data,assay='RNA',
                                             min.cells=0,
                                             min.features=0)
    
    spot_num=dim(seurat_obj)[2]
    cell_umi=seurat_obj@meta.data[,'nCount_RNA']
    cell_feature=seurat_obj@meta.data[,'nFeature_RNA']
    mito_percent=PercentageFeatureSet(seurat_obj, pattern="^MT-") %>% unlist()
    ribo_percent=PercentageFeatureSet(seurat_obj, pattern="^RP[SL]") %>% unlist()
    
    return(data.frame(spot_num,median(cell_umi,na.rm=TRUE),median(cell_feature,na.rm=TRUE),
                      median(mito_percent,na.rm=TRUE)/100,median(ribo_percent,na.rm=TRUE)/100))
    
},mc.cores=10)

dataset_visium=filter(qc_table,spot_num>500,Median_UMI>500,Median_Mito<0.2) %>% rownames()

# data preprocessing
for (dataset_1 in dataset_visium){
    
    ### create directory for each dataset
    if (!dir.exists(dataset_1)){
        dir.create(dataset_1)
    }
    
    ### read in data to Seurat
    data.dir=file.path('../0_public',dataset_1)
    if (file.exists(file.path(data.dir,'spatial','tissue_lowres_image.png'))){
        img=Read10X_Image(file.path(data.dir,'spatial'),image.name = "tissue_lowres_image.png")
        res_type='lowres'
    }else{
        img=Read10X_Image(file.path(data.dir,'spatial'),image.name = "tissue_hires_image.png")
        res_type='hires'
        img@scale.factors$lowres = img@scale.factors$hires
        img@spot.radius=(img@scale.factors$fiducial*img@scale.factors$lowres)/max(dim(img))
    }
    Seurat::DefaultAssay(object = img)='RNA'
    if (file.exists(file.path(data.dir,'filtered_feature_bc_matrix.h5'))){
        data=Read10X_h5(file.path(data.dir,'filtered_feature_bc_matrix.h5'))
    } else if (dir.exists(file.path(data.dir,'filtered_feature_bc_matrix'))){
        df=fread(file.path(data.dir,'filtered_feature_bc_matrix','features.tsv.gz'))
        if (ncol(df)==1) {data=Read10X(file.path(data.dir,'filtered_feature_bc_matrix'),gene.column=1)}
        else{data=Read10X(file.path(data.dir,'filtered_feature_bc_matrix'))}
    }
    
    # return(res_type)
    
    barcode_filtered=intersect(colnames(data),rownames(img@coordinates))
    data=data[,barcode_filtered]
    img@coordinates=img@coordinates[barcode_filtered,]
    
    # data filtering (filter out spots with less than 200 genes and genes expressed in less than 3 spots)
    seurat_obj <- Seurat::CreateSeuratObject(counts=data,assay='RNA',
                                             min.cells=3,
                                             min.features=101)
    img=img[colnames(x=seurat_obj)]
    seurat_obj[['image']]=img
    seurat_obj[["percent.mt"]]=PercentageFeatureSet(seurat_obj, pattern = "^MT-")  # tags
    seurat_obj@images$image@assay=res_type
    
    DefaultAssay(seurat_obj) <- 'RNA'
    
    # normalization remove mitochondrial and ribosomal gens
    gene_removed <- setdiff(rownames(seurat_obj),c(grep("^MT-",rownames(seurat_obj),value=T)
                                                  ))
    seurat_obj <- seurat_obj[gene_removed,]
    
    seurat_obj <- NormalizeData(seurat_obj,assay='RNA',
                                normalization.method="LogNormalize")
    seurat_obj <- FindVariableFeatures(seurat_obj,assay='RNA', selection.method = "vst", nfeatures=3000) 
    # seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c('percent.mt', 'nFeature_RNA'),
    #                         features=rownames(seurat_obj),split.by='orig.ident')
    
    # library(glmGamPoi)
    # seurat_obj <- SCTransform(seurat_obj,assay = "RNA",new.assay.name = "SCT", # method='glmGamPoi',
    #                           verbose = FALSE,reference.SCT.model = NULL,
    #                           return.only.var.genes = FALSE, 
    #                           variable.features.n = 3500,
    #                           vars.to.regress = c('percent.mt', 'nFeature_RNA')
    #                           )
    
    seurat_obj <- RenameCells(seurat_obj,add.cell.id=dataset_1)
       
    saveRDS(seurat_obj,file.path(dataset_1,paste0(dataset_1,'.rds')))
    
    return(NULL)
}