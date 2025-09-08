library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(Seurat)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

# MLCP definition
x=filter(CNV_result,CNVlabel=='Malignant') %>% rownames()
x=intersect(rownames(spot_type),x)
MLCP_spot_type=spot_type[x,]
MLCP_spot_type=MLCP_spot_type[apply(MLCP_spot_type,1,sum)>=0.7,]
MLCP_spot_type=MLCP_spot_type %>% apply(1,function(x){x/sum(x)})
MLCP_spot_type=t(MLCP_spot_type) %>% data.frame(check.names=FALSE)

sample_name=gsub('_[.A-Za-z0-9-]+$','',rownames(MLCP_spot_type))
tumor_type=gsub('_[0-9]+$','',sample_name)

MLCP_spot_type_scaled=tapply(rownames(MLCP_spot_type),tumor_type,function(x){
    MLCP_spot_type_x=MLCP_spot_type[x,]
    MLCP_spot_type_x=apply(MLCP_spot_type_x,2,function(y){
        y[quantile(y,0.99)]=quantile(y,0.99)
        return(y)
    })
    MLCP_spot_type_x=scale(MLCP_spot_type_x) %>% data.frame(check.names=FALSE) 
    return(MLCP_spot_type_x)
}) %>% dplyr::bind_rows()
MLCP_spot_type_scaled[MLCP_spot_type_scaled>=3]=3
MLCP_spot_type_scaled[MLCP_spot_type_scaled<=-1]=-1
MLCP_spot_type_scaled=asinh(MLCP_spot_type_scaled)

cluster_num=6:50
MLCP_cluster_df=sapply(cluster_num,function(x){
    MLCP_cluster_result <- kmeans(MLCP_spot_type_scaled,x,iter.max=100000,nstart=40)
    return(MLCP_cluster_result$cluster)
})
colnames(MLCP_cluster_df)=paste0('MLCP_',cluster_num)

MLCP_cluster_result[,'cluster_annotated']=recode(MLCP_cluster_result[,'cluster'],
              'MLCP_18'='MLCP_Tumor', # 
              'MLCP_26'='MLCP_Macro','MLCP_13'='MLCP_Macro', # 
              'MLCP_8'='MLCP_Macro/cDC','MLCP_9'='MLCP_Fibro/Endo/cDC',
              'MLCP_12'='MLCP_Fibro','MLCP_17'='MLCP_Fibro','MLCP_7'='MLCP_Fibro/Macro', # 
              'MLCP_23'='MLCP_Endo','MLCP_14'='MLCP_Endo','MLCP_16'='MLCP_Fibro/Endo', # 
              'MLCP_6'='MLCP_CD4 T','MLCP_11'='MLCP_cDC','MLCP_22'='MLCP_CD8 T', #
              'MLCP_20'='MLCP_Macro/cDC/CD8 T',
              'MLCP_10'='MLCP_B','MLCP_25'='MLCP_Macro/CD4 T',
              'MLCP_2'='MLCP_Macro/CD8 T',
              'MLCP_19'='MLCP_cDC/CD4 T',
              'MLCP_24'='MLCP_NK','MLCP_21'='MLCP_NK',
              'MLCP_1'='MLCP_Mast','MLCP_4'='MLCP_Mast',
              'MLCP_5'='MLCP_pDC','MLCP_3'='MLCP_panImmune','MLCP_15'='MLCP_panImmune')

spot_cluster <- group_by(data.frame(MLCP_spot_type_scaled,check.names=FALSE),as.factor(MLCP_cluster_result$cluster_annotated)) %>% 
    summarize_if(is.numeric,mean) %>% data.frame(row.names=1,check.names=FALSE)

spot_cluster=spot_cluster[c('MLCP_Tumor','MLCP_Endo','MLCP_Fibro/Endo','MLCP_Fibro','MLCP_Fibro/Macro','MLCP_Macro',
                            'MLCP_Macro/cDC','MLCP_cDC','MLCP_Fibro/Endo/cDC','MLCP_CD4 T','MLCP_CD8 T',
                            'MLCP_cDC/CD4 T','MLCP_Macro/CD4 T','MLCP_Macro/CD8 T','MLCP_Macro/cDC/CD8 T','MLCP_B','MLCP_NK',
                            'MLCP_pDC','MLCP_Mast','MLCP_panImmune'),
                          c('Parachymal cell','Endothelial cell','Fibroblast','Macrophage','cDC',
                             'CD4 T','CD8 T','B','NK','pDC','Mast cell')]

options(repr.plot.width=7,repr.plot.height=7)
ht_1=Heatmap(spot_cluster,name='Scaled cell\ntype comp',
             col=circlize::colorRamp2(c(-1,0.5,1),c('blue','white','red')),
             cluster_rows=FALSE,cluster_columns=FALSE,
             clustering_method_rows='ward.D2')

MLCP_tumor=reshape2::dcast(MLCP_comp,cancer_type~MLCP,value.var='Proportion',fun.aggregate=mean) %>% data.frame(row.names=1,check.names=FALSE)

MLCP_tumor=MLCP_tumor[,rownames(spot_cluster)]

s=unique(MLCP_comp$MLCP)
p_mat=sapply(s,function(x){
    filter(MLCP_comp,MLCP==x) %>% OneVsRes(group_col='cancer_type',value_col='Proportion')
},USE.NAMES=TRUE)
colnames(p_mat)=s

p_mat=p_mat[rownames(MLCP_tumor),colnames(MLCP_tumor)] %>% data.frame(check.names=FALSE)
p_mat[is.na(p_mat)]=1
p_mat=t(p_mat)

options(repr.plot.width=7,repr.plot.height=7)
ht_2=Heatmap(t(scale(MLCP_tumor)),name='Scaled\nLCP comp',col=circlize::colorRamp2(c(-2,0,4),c('#1a318b','#ffffff','#9a133d')),
             cluster_rows=FALSE,cluster_columns=FALSE,
             cell_fun=function(j, i, x, y, w, h, fill) {
                 if (p_mat[i,j]<0.001) {
                     grid.text("**", x, y)
                 } 
                 else if (p_mat[i,j]<0.01) {
                     grid.text("*", x, y)
                 }
             },
             right_annotation=rowAnnotation(' '=anno_barplot(table(MLCP_cluster_result$cluster_annotated) %>% 
                                                             as.data.frame.array() %>% .[rownames(spot_cluster),])))

options(repr.plot.width=12,repr.plot.height=7)
ht_1+ht_2

seurat_obj=CreateSeuratObject(counts=t(MLCP_spot_type_scaled))
seurat_obj@assays$RNA@scale.data=as.matrix(seurat_obj@assays$RNA@data)
seurat_obj@assays$RNA@var.features=rownames(seurat_obj)
seurat_obj=RunPCA(seurat_obj)
seurat_obj=FindNeighbors(seurat_obj,dims=1:10)
seurat_obj=RunTSNE(seurat_obj,dims=1:10)

# SLCP definition
x=filter(CNV_result,CNVlabel=='Stromal') %>% rownames()
x=intersect(rownames(spot_type),x)
SLCP_spot_type=spot_type[x,]
SLCP_spot_type=SLCP_spot_type[apply(SLCP_spot_type,1,sum)>=0.7,]
SLCP_spot_type=SLCP_spot_type %>%  select(-`Parachymal cell`) %>% apply(1,function(x){x/sum(x)})
SLCP_spot_type=t(SLCP_spot_type) %>% data.frame(check.names=FALSE)

sample_name=gsub('_[.A-Za-z0-9-]+$','',rownames(SLCP_spot_type))
tumor_type=gsub('_[0-9]+$','',sample_name)

SLCP_spot_type_scaled=tapply(rownames(SLCP_spot_type),tumor_type,function(x){
    SLCP_spot_type_x=SLCP_spot_type[x,]
    SLCP_spot_type_x=apply(SLCP_spot_type_x,2,function(y){
        y[quantile(y,0.99)]=quantile(y,0.99)
        return(y)
    })
    SLCP_spot_type_x=scale(SLCP_spot_type_x) %>% data.frame(check.names=FALSE) 
    return(SLCP_spot_type_x)
}) %>% dplyr::bind_rows()

SLCP_spot_type_scaled[SLCP_spot_type_scaled>=3]=3
SLCP_spot_type_scaled[SLCP_spot_type_scaled<=-1]=-1

cluster_num=6:50
SLCP_cluster_df=sapply(cluster_num,function(x){
    SLCP_cluster_result <- kmeans(SLCP_spot_type_scaled,x,iter.max=100000,nstart=40)
    return(SLCP_cluster_result$cluster)
})
colnames(SLCP_cluster_df)=paste0('SLCP_',cluster_num)

SLCP_cluster_result[,'cluster_annotated']=recode(SLCP_cluster_result[,'cluster'],
                'SLCP_4'='SLCP_Fibro','SLCP_6'='SLCP_Endo','SLCP_16'='SLCP_Fibro/Endo','SLCP_13'='SLCP_Macro',
                'SLCP_15'='SLCP_B','SLCP_14'='SLCP_cDC/CD4 T','SLCP_7'='SLCP_CD4 T','SLCP_9'='SLCP_Endo/cDC',
                'SLCP_10'='SLCP_Fibro/cDC','SLCP_8'='SLCP_pDC','SLCP_3'='SLCP_Mast','SLCP_5'='SLCP_NK',
                'SLCP_1'='SLCP_CD8 T','SLCP_12'='SLCP_Macro/CD8 T','SLCP_11'='SLCP_Fibro/Macro',
                'SLCP_2'='SLCP_cDC/Macro')

spot_cluster <- group_by(data.frame(SLCP_spot_type_scaled,check.names=FALSE),as.factor(SLCP_cluster_result$cluster_annotated)) %>% 
    summarize_if(is.numeric,mean) %>% data.frame(row.names=1,check.names=FALSE)

spot_cluster=spot_cluster[c('SLCP_Endo','SLCP_Fibro/Endo','SLCP_Fibro','SLCP_Macro','SLCP_Fibro/Macro','SLCP_Fibro/cDC','SLCP_Endo/cDC',
                            'SLCP_cDC/Macro','SLCP_CD4 T','SLCP_CD8 T','SLCP_cDC/CD4 T','SLCP_Macro/CD8 T',
                            'SLCP_B','SLCP_NK','SLCP_pDC','SLCP_Mast'),
                          c('Endothelial cell','Fibroblast','Macrophage','cDC',
                             'CD4 T','CD8 T','B','NK','pDC','Mast cell')]

options(repr.plot.width=7,repr.plot.height=7)
ht_3=Heatmap(spot_cluster,name='Percantage',clustering_method_rows='ward.D2', 
             cluster_rows=FALSE,cluster_columns=FALSE,
             col=circlize::colorRamp2(c(-0.5,0.3,1.5),c('blue','white','red')))

SLCP_tumor=reshape2::dcast(SLCP_comp,cancer_type~SLCP,value.var='Proportion',fun.aggregate=mean) %>% data.frame(row.names=1,check.names=FALSE)

SLCP_tumor=SLCP_tumor[,rownames(spot_cluster)]

s=unique(SLCP_comp$SLCP)
p_mat=sapply(s,function(x){
    filter(SLCP_comp,SLCP==x) %>% OneVsRes(group_col='cancer_type',value_col='Proportion')
},USE.NAMES=TRUE)
colnames(p_mat)=s

p_mat=p_mat[rownames(SLCP_tumor),colnames(SLCP_tumor)] %>% data.frame(check.names=FALSE)
p_mat[is.na(p_mat)]=1
p_mat=t(p_mat)

options(repr.plot.width=7,repr.plot.height=7)
ht_4=Heatmap(t(scale(SLCP_tumor)),name='Scaled\nLCP comp',col=circlize::colorRamp2(c(-2,0,4),c('#1a318b','#ffffff','#9a133d')),
             cluster_rows=FALSE,cluster_columns=FALSE,
             cell_fun=function(j, i, x, y, w, h, fill) {
                 if (p_mat[i, j] < 0.001) {
                     grid.text("**", x, y)
                 } else if (p_mat[i, j] < 0.01) {
                     grid.text("*", x, y)
                 }},
             right_annotation=rowAnnotation(' '=anno_barplot(table(SLCP_cluster_result$cluster_annotated) %>% 
                                                             as.data.frame.array() %>% .[rownames(spot_cluster),])))
seurat_obj=CreateSeuratObject(counts=t(SLCP_spot_type_scaled))
seurat_obj@assays$RNA@scale.data=as.matrix(seurat_obj@assays$RNA@data)
seurat_obj@assays$RNA@var.features=rownames(seurat_obj)
seurat_obj=RunPCA(seurat_obj)
seurat_obj=FindNeighbors(seurat_obj,dims=1:10)
seurat_obj=RunTSNE(seurat_obj,ims=1:10)

# methods for comparison
OneVsRes=function(df,group_col,value_col){
    sapply(unique(df[,group_col]),function(x){
        
        m=df[ which(df[,group_col]==x), value_col]
        n=df[ which(df[,group_col]!=x), value_col]
        
        df_test=data.frame(value=c(m,n),
          group=as.factor(c(rep('A',length.out=length(m)),rep('B',length.out=length(n)))))
        colnames(df_test)=c('value','group')
        df_test=data.frame(df_test)
        
        wilcox.test(m,n)$p.value
        
    },USE.NAMES=TRUE)
}