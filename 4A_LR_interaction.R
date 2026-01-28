library(dplyr)
library(ggplot2)
library(aplot)
library(ComplexHeatmap)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')
source('https://github.com/KaWingLee9/in_house_tools/blob/main/meta_analysis/custom_fun.R')

# LR interaction score of spots
CalSpotLRScore_V2 <- function(exp_mat,lr_ls,min_spot=10){
    
    # LR score for spots
    lr_score_spot=sapply(1:length(lr_ls),function(i){

        lr_pair=lr_ls[[i]]
        
        # filter out LR with no genes detected
        if (length(setdiff(unlist(lr_pair),rownames(exp_mat)))!=0){
            return(NULL)
        }
        
        # geometric mean of all genes in LR pair
        L_genes=lr_pair[['L']]
        if (length(L_genes)==1){
            L_score=exp_mat[L_genes,,drop=FALSE] %>% as.matrix() %>% t() %>% data.frame()
        }else{
            L_score=apply(exp_mat[L_genes,,drop=FALSE],2,function(x){
                prod(x)^(1/length(x))
            }) %>% data.frame()
        }

        R_genes=lr_pair[['R']]
        if (length(R_genes)==1){
            R_score=exp_mat[R_genes,,drop=FALSE] %>% as.matrix() %>% t() %>% data.frame()
        }else{
            R_score=apply(exp_mat[R_genes,,drop=FALSE],2,function(x){
                prod(x)^(1/length(x))
            }) %>% data.frame()
        }
        
        df=sqrt(L_score*R_score)
        
        colnames(df)=names(lr_ls)[i]
        return(df)
    }) %>% dplyr::bind_cols() %>% data.frame(row.names=colnames(exp_mat),check.names=FALSE)

    # filter efficient L-R: number of spots with positive score not lower than `min_spot`
    LR_used=which(apply(lr_score_spot,2,function(x){sum(x>0)})>=min_spot) %>% names()
    lr_score_spot=lr_score_spot[,LR_used]
    
    return(as.sparse(t(lr_score_spot)))
    
}

seurat_obj=readRDS(file.path('../1_data_preprocessing',file.path(dataset_1,paste0(dataset_1,'.rds'))))
niche_cluster_result_sample=filter(niche_cluster_result,sample==dataset_1)

exp_mat=seurat_obj@assays$RNA@data
lr_score_spot=CalSpotLRScore_V2(exp_mat,lr_ls,min_spot=0)
lr_score_spot=lr_score_spot %>% as.matrix() %>% t()

# LR interaction  of each sample

ind=intersect(rownames(lr_score_spot),rownames(niche_cluster_result))
lr_score_spot=lr_score_spot[ind,]

seurat_obj=seurat_obj[,ind]
seurat_obj[['CCI']]=CreateAssayObject(t(lr_score_spot))

Idents(seurat_obj)=niche_cluster_result[ind,'Niche_combined']

future::plan("multicore",workers=10) 
options(future.globals.maxSize=10000*1024^2)
sce.markers=FindAllMarkers(object=seurat_obj,
                            only.pos=FALSE,min.pct=0,logfc.threshold=0,return.thresh=1)

# meta analysis
l=lapply(dataset_ls,function(x){

    res=try({y <- read.csv(paste0(dir,x,'.csv'),row.names=1,check.names=FALSE)},silent=TRUE)
    if (inherits(res,'try-error')) return(NULL)
    return(y)
})

l=l[!unlist(lapply(l,is.null))]

for (n in paste0('Niche_',1:13)){
    
    l_niche=parallel::mclapply(l,function(x){
        y=filter(x,cluster==n)
        rownames(y)=y[,'gene']
        return(y)
    },mc.cores=30) %>% .[lapply(.,nrow)!=0]

    test_result=CombRank_DFLs(l_niche,p_col='p_val',ES_col='avg_log2FC',
                              min_num=0,min_ratio=0.5,method='RankSum',quantile_est=1/4)

    write.csv(test_result,paste0('./FindMarkers_Cellchat_meta/',n,'.csv'))
    
}

meta_anlysis_result=lapply(dir('./FindMarkers_CellphoneDB_meta/'),function(x){
    y=read.csv(paste0('./FindMarkers_CellphoneDB_meta/',x),row.names=1,check.names=FALSE)
    y[,'Niche']=gsub('.csv','',x)
    y[,'lr_pair']=rownames(y)
    return(y)
}) %>% dplyr::bind_rows()

rownames(meta_anlysis_result)=1:nrow(meta_anlysis_result)

meta_anlysis_result_CellChat_filtered=meta_anlysis_result_CellChat %>% filter(available_study>=1/3*total_study,
                                                                              signed_combined_rank>=0.5,
                                                                              quantile_pval<=0.05)

# LR number of the niches
dd=meta_anlysis_result_CellChat_filtered %>% 
    group_by(lr_pair) %>% summarize(tumor_type=list(Niche))

options(repr.plot.width=25,repr.plot.height=8)
p=ggplot(data=dd,
       aes(x=tumor_type))+
    scale_x_upset(order_by='freq',sets=paste0('Niche_',1:13))+
    geom_bar()+
    xlab('')+
    ylab('LR number')+
    theme_classic()+
    scale_y_continuous(expand=expansion(mult=0))+
    theme_combmatrix(combmatrix.label.text=element_text(angle=180,vjust=0))

p

# enrichment anslysis of LR pathway
df=dplyr::left_join(meta_anlysis_result_CellChat,cellchat_table,by=c('lr_pair'='interaction_name_2'),multiple="all")
df=df %>% mutate(selected= 
                        (available_study>=1/3*total_study) & 
                        (signed_combined_rank>=0.5) & 
                        (sig_study_ratio>=0.25))
df=df %>% filter(! stringr::str_detect(lr_pair,'[...]'))

test_result=unique(df[,c('Niche','pathway_name')])
test_result[,c('OR','P')]=NA

x=apply(test_result,1,function(x){
    
    a=df %>% filter(Niche==x[[1]],pathway_name==x[[2]],selected==TRUE) %>% nrow()
    b=df %>% filter(Niche==x[[1]],pathway_name!=x[[2]],selected==TRUE) %>% nrow()
    c=df %>% filter(Niche==x[[1]],pathway_name==x[[2]],selected==FALSE) %>% nrow()
    d=df %>% filter(Niche==x[[1]],pathway_name!=x[[2]],selected==FALSE) %>% nrow()
    
    test_result=fisher.test(rbind(c(a,b),c(c,d)),alternative='two.sided')
    return(c(a,b,c,d,test_result[['estimate']],test_result[['p.value']]))
    
})

test_result[,c('a','b','c','d','OR','P')]=t(x)

x=table(cellchat_table$pathway_name) %>% c()
x=x[x>=5]

test_result_filtered=test_result %>% filter(OR>1,P<=0.05,a>=2
                                           )
test_result_filtered[,'OR1']=test_result_filtered[,'OR']
test_result_filtered[ test_result_filtered[,'OR1']>=5 ,'OR1']=5
test_result_filtered=test_result_filtered %>% filter(pathway_name %in% names(x))

options(repr.plot.height=3,repr.plot.width=7)
p=ggplot(test_result_filtered,aes(y=Niche,x=pathway_name))+
    geom_point(color='red')+
    # coord_flip()+
    xlab('')+
    ylab('')+
    theme_bw()+
    scale_y_discrete(position='right')+
    theme(axis.text.x=element_text(angle=45,hjust=1))
p

p=OrderedPlot(p,x='pathway_name',y='Niche',cluster_value='OR1',cluster_var=NULL,
              cluster_row=TRUE,show_row_dend=TRUE,row_dend_width=0.05,row_dend_direction='left',
              cluster_column=TRUE,show_column_dend=TRUE,
              method='ward.D2'
              )
p

# showing the expression of ligands and receptors respectively

# meta analysis
for (dataset in dataset_used){

    print(dataset)

    seurat_obj=readRDS( paste0('../1_data_preprocessing/',dataset,'/',dataset,'.rds') )
    
    ind=intersect(colnames(seurat_obj),rownames(niche_cluster_result))
    seurat_obj=seurat_obj[,ind]
    
    Idents(seurat_obj)=niche_cluster_result[colnames(seurat_obj),'Niche_combined']
    
    x=Idents(seurat_obj) %>% table()
    
    x=x[x>10] %>% names()
    
    # seurat_obj=subset(seurat_obj,idents=x)
    
    future::plan("multicore",workers=60) 
    options(future.globals.maxSize=10000*1024^2)
    sce.markers=FindAllMarkers(object=seurat_obj,
                                only.pos=FALSE,min.pct=0,logfc.threshold=0,return.thresh=1)
    
    sce.markers[,'dataset']=dataset

    sce.markers=filter(sce.markers,cluster %in% x)
    
    write.csv(sce.markers,paste0('./FindMarkers/',dataset,'.csv'))
    
}

dir='./FindMarkers/'
dataset_ls=gsub('.csv','',grep('BRCA|COAD|HNSC|CSCC|LUAD|GBM|OV|PRAD|PAAD|KIRC|LIHC|LUSC',dir(dir),value=TRUE))

l=lapply(dataset_ls,function(x){

    res=try({y <- read.csv(paste0(dir,x,'.csv'),row.names=1,check.names=FALSE)},silent=TRUE)
    if (inherits(res,'try-error')) return(NULL)
    return(y)
})

l=l[!unlist(lapply(l,is.null))]

for (n in paste0('Niche_',1:13)){
    
    l_niche=parallel::mclapply(l,function(x){
        y=filter(x,cluster==n)
        rownames(y)=y[,'gene']
        return(y)
    },mc.cores=20) %>% .[lapply(.,nrow)!=0]

    test_result=CombRank_DFLs(l_niche,p_col='p_val',ES_col='avg_log2FC',
                              min_num=0,min_ratio=0.5,method='RankSum')

    write.csv(test_result,paste0('./FindMarkers_meta/',n,'.csv'))
    
}

meta_anlysis_result=lapply(dir('./FindMarkers_meta/'),function(x){
    y=read.csv(paste0('./FindMarkers_meta/',x),row.names=1,check.names=FALSE)
    y[,'Niche']=gsub('.csv','',x)
    y[,'gene']=rownames(y)
    return(y)
}) %>% dplyr::bind_rows()

rownames(meta_anlysis_result)=1:nrow(meta_anlysis_result)

# showing the expression of ligands and receptors respectively
df=meta_anlysis_result
df[,'Niche']=factor(df[,'Niche'],c('Niche_1','Niche_2','Niche_3','Niche_4','Niche_5','Niche_6',
                                   'Niche_7','Niche_8','Niche_9','Niche_10','Niche_11','Niche_12','Niche_13'))

selected_chemokine_LR=c('CCL15-CCR1','CCL7-CCR1','CCL19-CCR7','CCL21-CCR7','CCL5-CCR1','CCL5-CCR3','CCL5-CCR4','CCL5-CCR5','CCL3-CCR1','CCL3-CCR5',
                        'CXCL12-ACKR3','CXCL12-CXCR4','CXCL13-CXCR5','CXCL16-CXCR6','CXCL9-CXCR3','CXCL10-CXCR3','CXCL10-CXCR3')

df_filter=df %>% filter(available_study>=1/3*total_study# ,sig_study_ratio>=0.3
                       )

LR_pairs=selected_chemokine_LR

L=lapply(LR_pairs,function(x) {strsplit(x,'-')[[1]][1]  %>% strsplit(.,'[+]') %>% unlist()})
R=lapply(LR_pairs,function(x) {strsplit(x,'-')[[1]][2]  %>% strsplit(.,'[+]') %>% unlist()})

link_df=lapply(1:length(L),function(x){
    expand.grid(L[[x]],R[[x]])
}) %>% dplyr::bind_rows()

link_df=link_df %>% filter((Var1 %in% df_filter$gene) & (Var2 %in% df_filter$gene))

# reset ligand and receptor order
df_order=link_df %>% reshape2::dcast(Var1~Var2) %>% data.frame(row.names=1)
df_order[!is.na(df_order)]=1
df_order[is.na(df_order)]=0
df_order=mutate_all(df_order,as.numeric)

p=LinkedPlot(df_filter,link_df,
             x_col='Niche',
             y_col='gene',
             fill_col='signed_combined_rank',
             size_col='sig_study_ratio',
             color_column=1,align='center')

p[[1]]=p[[1]]+ggtitle('Ligand')+
#     scale_fill_gradientn(colours=c(paletteer::paletteer_c("grDevices::RdBu",11)),
#                          values=scales::rescale(c(max(p[[1]]$data$signed_combined_rank), 0.8,0.7,0.6,0.5,0,
#                                   -0.5,-0.6,-0.7,-0.8,min(p[[1]]$data$signed_combined_rank) )))+
    scale_fill_gradientn(colours=c(paletteer::paletteer_c("grDevices::Reds",11)),
                         values=scales::rescale(c(max(p[[1]]$data$signed_combined_rank), 0.8,0.7,0.6,0.5,0.4,
                                  0,min(p[[1]]$data$signed_combined_rank) )))+
    theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))
p[[3]]=p[[3]]+ggtitle('Receptor')+
#     scale_fill_gradientn(colours=colorRampPalette(paletteer::paletteer_c("grDevices::RdBu",30))(11),
#                          values=scales::rescale(c(max(p[[3]]$data$signed_combined_rank), 0.8,0.7,0.6,0.5,0,
#                                   -0.5,-0.6,-0.7,-0.8,min(p[[3]]$data$signed_combined_rank) )))+
    scale_fill_gradientn(colours=c(paletteer::paletteer_c("grDevices::Blues",8)),
                         values=scales::rescale(c(max(p[[3]]$data$signed_combined_rank),0.8,0.7,0.6,0.5,0.4,
                                  0,min(p[[3]]$data$signed_combined_rank) )))+
    theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))

p