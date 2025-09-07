library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(spacexr)
library(Matrix)
library(doParallel)
library(Seurat)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/meta_analysis/custom_fun.R')

# CSIDE

exp_group=c('Niche_1','Niche_2','Niche_3','Niche_4','Niche_5','Niche_6')

for (compare_group_experimental in exp_group){

    dir.create(compare_group_experimental)
    compare_group_control=c('Niche_1')

    x=table(niche_cluster_result[,c('sample','Niche_combined')]) %>% as.data.frame.array()

    cell_num_filter=50
    dataset_used=which((x[,compare_group_control,drop=FALSE] %>% apply(1,sum)>=cell_num_filter) & 
                    (x[,compare_group_experimental,drop=FALSE] %>% apply(1,sum)>=cell_num_filter)) %>% names()


    for (dataset in dataset_used
        ){
        tumor_type=gsub('_[0-9]+$','',dataset)

        dconv_result <- read.csv(paste0('../2_deconvolution/RCTD_new_v2/',dataset,'_RCTD_result.csv'),row.names =1,check.names=FALSE)
        tumor_type=gsub('_[0-9]+','',dataset)
        tumor_cell=tumor_cell_type[[tumor_type]]    
        dconv_result=TumorCombn_v2(dconv_result,tumor_cell=tumor_cell)

        scRNA_sce=readRDS(paste0('../2_deconvolution/scRNA_used/',tumor_type,'_new.rds'))
        x=Idents(scRNA_sce) %>% as.character()
        x[x %in% tumor_cell]='Parachymal cell'
        Idents(scRNA_sce)=as.factor(x)

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

        c=aggregate_cell_types(RCTD_result,barcodes=names(explanatory.variable),doublet_mode=FALSE)
        cell_types_present=c[c>=30] %>% names()
        cell_types_present=c(cell_types_present,'Parachymal cell') %>% unique()

        res=try({
            suppressWarnings(RCTD_result <-  run.CSIDE.single(RCTD_result,explanatory.variable,gene_threshold=5e-05,
                                                            cell_types=cell_types_present,
                                                            cell_types_present=NULL,
                                                            cell_type_threshold=0,fdr=1,
                                                            doublet_mode=FALSE,log_fc_thresh=0,
                                                            weight_threshold=0.05))
        })
        
        if(inherits(res,"try-error")) next
        
        df=RCTD_result@de_results$all_gene_list$`Parachymal cell`
        
    }
}

# meta analysis
for (niche in paste0('Niche_',1:6)){
    
    dir=niche
    
    dataset_ls=gsub('.csv','',grep('BRCA|COAD|LUAD|GBM|PAAD|KIRC|LIHC|LUSC|OV|HNSC|CSCC|PRAD'
                               ,dir(dir),value=TRUE))

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
    write.csv(test_result,paste0(niche,'_meta_RankSum_result.csv'))
}

# enrichment analysis
up_enrich_result=lapply(unique(Para_gene_set[,'source']),function(source_1){
    Para_gene_set_sub=Para_gene_set %>% filter(source==source_1)
    
    lapply(paste0('Niche_',2:6),function(niche){

        test_result=read.csv(paste0(dir,'_meta/',niche,'_meta_RankSum_result.csv'),row.names=1,check.names=FALSE)
        test_result[,'rank']=rank(test_result[,'signed_combined_rank'])

        test_result[,'gene_name']=rownames(test_result)
        test_result[,'gene_type']=0
        test_result[ (test_result[,'signed_combined_rank']>=0.5) & (test_result[,'sig_study_ratio']>=sig_study_thres) ,'gene_type']=1
        test_result[ (test_result[,'signed_combined_rank']<=-0.5) & (test_result[,'sig_study_ratio']>=sig_study_thres) ,'gene_type']=-1
        test_result[,'gene_type']=as.character(test_result[,'gene_type'])

        up_gene_list=rownames(test_result[ test_result[,'gene_type']==1 ,])

        test_result_up=test_result[up_gene_list,] %>% arrange(desc(signed_combined_rank)) %>% head(10)

        up_enrich_result=enricher(up_gene_list,TERM2GENE=Para_gene_set_sub,
                                  pAdjustMethod='fdr',pvalueCutoff=0.05,qvalueCutoff=0.05)
        up_enrich_result=data.frame(up_enrich_result)

        if (nrow(up_enrich_result)==0){
            return(NULL)
        }

        up_enrich_result[,'Niche']=niche

        up_enrich_result[,'direct']='Higher_expressed'
        
        rownames(up_enrich_result)=NULL

        return(up_enrich_result)

    }) %>% dplyr::bind_rows()
    
}) %>% dplyr::bind_rows()

down_enrich_result=lapply(unique(Para_gene_set[,'source']),function(source_1){
    Para_gene_set_sub=Para_gene_set %>% filter(source==source_1)
    
    lapply(paste0('Niche_',2:6),function(niche){

        test_result=read.csv(paste0(dir,'_meta/',niche,'_meta_RankSum_result.csv'),row.names=1,check.names=FALSE)
        test_result[,'rank']=rank(test_result[,'signed_combined_rank'])

        test_result[,'gene_name']=rownames(test_result)
        test_result[,'gene_type']=0
        test_result[ (test_result[,'signed_combined_rank']>=0.5) & (test_result[,'sig_study_ratio']>=sig_study_thres) ,'gene_type']=1
        test_result[ (test_result[,'signed_combined_rank']<=-0.5) & (test_result[,'sig_study_ratio']>=sig_study_thres) ,'gene_type']=-1
        test_result[,'gene_type']=as.character(test_result[,'gene_type'])

        down_gene_list=rownames(test_result[ test_result[,'gene_type']==-1 ,])
        test_result_down=test_result[down_gene_list,] %>% arrange((signed_combined_rank)) %>% head(10)

        down_enrich_result=enricher(down_gene_list,TERM2GENE=Para_gene_set_sub,
                                    pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff=0.05)
        down_enrich_result=data.frame(down_enrich_result)

        if (nrow(down_enrich_result)==0){
            return(NULL)
        }

        down_enrich_result[,'Niche']=niche

        down_enrich_result[,'direct']='Lower_expressed'
        rownames(down_enrich_result)=NULL

        return(down_enrich_result)

    }) %>% dplyr::bind_rows()
    
}) %>% dplyr::bind_rows()

enrich_result=rbind(up_enrich_result,down_enrich_result)