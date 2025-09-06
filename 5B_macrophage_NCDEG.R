library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(spacexr)
library(Matrix)
library(doParallel)
library(Seurat)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/meta_analysis/custom_fun.R')

# CSIDE between Niche_4 and Niche_11
compare_group_control=c('Niche_4') # Macrophage without other immune cells
compare_group_experimental=c('Niche_11') # Macrophage with other immune cell

x=table(niche_cluster_result[,c('sample','Niche_combined')]) %>% as.data.frame.array()

cell_num_filter=50
dataset_used=which((x[,compare_group_control,drop=FALSE] %>% apply(1,sum)>=cell_num_filter) & 
                   (x[,compare_group_experimental,drop=FALSE] %>% apply(1,sum)>=cell_num_filter)) %>% names()

for (dataset in dataset_used
      ){
    tumor_type=gsub('_[0-9]+$','',dataset)

    dconv_result <- read.csv(paste0('../2_deconvolution/RCTD_new_v2/',dataset,'_RCTD_result.csv'),row.names =1,check.names=FALSE)

    scRNA_sce=readRDS(paste0('../2_deconvolution/scRNA_used/',tumor_type,'_new.rds'))

    if ('Undefined' %in% unique(Idents(scRNA_sce))){
        scRNA_sce=subset(scRNA_sce,idents='Undefined',invert=TRUE)
    }

    stRNA_sce=readRDS(paste0('../1_data_preprocessing/',dataset,'/',dataset,'.rds'))

    stRNA_sce@meta.data$Niche_combined=niche_cluster_result[colnames(stRNA_sce),'Niche_combined']
    stRNA_sce_subset <- subset(stRNA_sce,subset = Niche_combined %in% c(compare_group_control,compare_group_experimental))

    RCTD_result <- create.RCTD(SpatialRNA(data.frame(x=rep(1,dim(stRNA_sce_subset)),y=rep(1,dim(stRNA_sce_subset))),
                                          GetAssayData(stRNA_sce_subset,assay='RNA',slot='counts'),use_fake_coords=TRUE), 
                               Reference(GetAssayData(scRNA_sce,assay='RNA',slot='counts'),Idents(scRNA_sce)),
                               CELL_MIN_INSTANCE=10,
                               max_cores=1,test_mode=FALSE)
    RCTD_result <- import_weights(RCTD_result,dconv_result[colnames(stRNA_sce_subset),])

    RCTD_result@config$RCTDmode <- "full"
    RCTD_result@config$max_cores <- 15
    RCTD_result@config$doublet_mode <- "full"

    explanatory.variable <- rep(0,length.out=length(rownames(stRNA_sce_subset@meta.data)))
    names(explanatory.variable) <- rownames(stRNA_sce_subset@meta.data)
    explanatory.variable[stRNA_sce_subset@meta.data$Niche_combined %in% compare_group_experimental]=1

    cell_types='Macrophage'
    c=aggregate_cell_types(RCTD_result,barcodes=names(explanatory.variable),doublet_mode=FALSE)
    cell_types_present=c[c>=30] %>% names()
    cell_types_present=c(cell_types_present,'Macrophage') %>% unique()

    res=try({
        suppressWarnings(RCTD_result <-  run.CSIDE.single(RCTD_result,explanatory.variable,gene_threshold=5e-05,
                                                        cell_types=cell_types_present,
                                                        cell_types_present=NULL,
                                                        cell_type_threshold=0,fdr=1,
                                                        doublet_mode=FALSE,log_fc_thresh=0,
                                                        weight_threshold=0.05))
    })
    if(inherits(res,"try-error")) next
    
    df=RCTD_result@de_results$all_gene_list$Macrophage
    
} 

# meta analysis
l=lapply(dataset_ls,function(x){

    res=try({y <- read.csv(paste0(dir,'/',x,'.csv'),row.names=1,check.names=FALSE)},silent=TRUE)
    if (inherits(res,'try-error')) return(NULL)
    return(y)
})

l=l[!unlist(lapply(l,is.null))]

gene_num=lapply(l,rownames) %>% unlist() %>% table() %>% c()
gene_filtered=gene_num[gene_num>= 1/3 * length(l)] %>% names()
l=lapply(l,function(x) x[intersect(gene_filtered,rownames(x)),])

test_result=CombRank_DFLs(l,p_col='p_val',ES_col='log_fc',min_num=0,min_ratio=0,method='RankSum',quantile_est=1/4)

# CSIDE between Niche_10 and Niche_11
compare_group_control=c('Niche_10') # Macrophage without other immune cells
compare_group_experimental=c('Niche_11') # Macrophage with other immune cell

x=table(niche_cluster_result[,c('sample','Niche_combined')]) %>% as.data.frame.array()

cell_num_filter=50
dataset_used=which((x[,compare_group_control,drop=FALSE] %>% apply(1,sum)>=cell_num_filter) & 
                   (x[,compare_group_experimental,drop=FALSE] %>% apply(1,sum)>=cell_num_filter)) %>% names()

for (dataset in dataset_used
      ){
    tumor_type=gsub('_[0-9]+$','',dataset)

    dconv_result <- read.csv(paste0('../2_deconvolution/RCTD_new_v2/',dataset,'_RCTD_result.csv'),row.names =1,check.names=FALSE)

    scRNA_sce=readRDS(paste0('../2_deconvolution/scRNA_used/',tumor_type,'_new.rds'))

    if ('Undefined' %in% unique(Idents(scRNA_sce))){
        scRNA_sce=subset(scRNA_sce,idents='Undefined',invert=TRUE)
    }

    stRNA_sce=readRDS(paste0('../1_data_preprocessing/',dataset,'/',dataset,'.rds'))

    stRNA_sce@meta.data$Niche_combined=niche_cluster_result[colnames(stRNA_sce),'Niche_combined']
    stRNA_sce_subset <- subset(stRNA_sce,subset = Niche_combined %in% c(compare_group_control,compare_group_experimental))

    RCTD_result <- create.RCTD(SpatialRNA(data.frame(x=rep(1,dim(stRNA_sce_subset)),y=rep(1,dim(stRNA_sce_subset))),
                                          GetAssayData(stRNA_sce_subset,assay='RNA',slot='counts'),use_fake_coords=TRUE), 
                               Reference(GetAssayData(scRNA_sce,assay='RNA',slot='counts'),Idents(scRNA_sce)),
                               CELL_MIN_INSTANCE=10,
                               max_cores=1,test_mode=FALSE)
    RCTD_result <- import_weights(RCTD_result,dconv_result[colnames(stRNA_sce_subset),])

    RCTD_result@config$RCTDmode <- "full"
    RCTD_result@config$max_cores <- 15
    RCTD_result@config$doublet_mode <- "full"

    explanatory.variable <- rep(0,length.out=length(rownames(stRNA_sce_subset@meta.data)))
    names(explanatory.variable) <- rownames(stRNA_sce_subset@meta.data)
    explanatory.variable[stRNA_sce_subset@meta.data$Niche_combined %in% compare_group_experimental]=1

    cell_types='Macrophage'
    c=aggregate_cell_types(RCTD_result,barcodes=names(explanatory.variable),doublet_mode=FALSE)
    cell_types_present=c[c>=30] %>% names()
    cell_types_present=c(cell_types_present,'Macrophage') %>% unique()

    res=try({
        suppressWarnings(RCTD_result <-  run.CSIDE.single(RCTD_result,explanatory.variable,gene_threshold=5e-05,
                                                        cell_types=cell_types_present,
                                                        cell_types_present=NULL,
                                                        cell_type_threshold=0,fdr=1,
                                                        doublet_mode=FALSE,log_fc_thresh=0,
                                                        weight_threshold=0.05))
    })
    if(inherits(res,"try-error")) next
    
    df=RCTD_result@de_results$all_gene_list$Macrophage
    
}

# meta analysis
l=lapply(dataset_ls,function(x){

    res=try({y <- read.csv(paste0(dir,'/',x,'.csv'),row.names=1,check.names=FALSE)},silent=TRUE)
    if (inherits(res,'try-error')) return(NULL)
    return(y)
})

l=l[!unlist(lapply(l,is.null))]

gene_num=lapply(l,rownames) %>% unlist() %>% table() %>% c()
gene_filtered=gene_num[gene_num>= 1/3 * length(l)] %>% names()
l=lapply(l,function(x) x[intersect(gene_filtered,rownames(x)),])

test_result=CombRank_DFLs(l,p_col='p_val',ES_col='log_fc',min_num=0,min_ratio=0,method='RankSum',quantile_est=1/4)

