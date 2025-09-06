library(dplyr)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(DWLS)

# DWLS
seurat_obj=readRDS(paste0('../2_deconvolution/scRNA_used/',cancer_type,'_new.rds'))

if ('Undefined' %in% unique(Idents(seurat_obj))){
    seurat_obj=subset(seurat_obj,idents=c('Undefined'),invert=TRUE)
}

cell_type_label=Idents(seurat_obj) %>% as.character()
cell_type_label[which(cell_type_label %in% tumor_cell_type[[cancer_type]])]='Parachymal cell'
cell_type_label=sub(' ','.',cell_type_label)

sc_sig=buildSignatureMatrixUsingSeurat(scdata=seurat_obj@assays$RNA@data,
                                       id=cell_type_label,
                                       path=path.dir)

deconv_result=apply(bulk_RNA_mat,2,function(x){
  tr=trimData(sc_sig,x)
  solveDampenedWLS(tr$sig, tr$bulk)
})

df=lapply(dir(d),function(file){
    
    deconv_result=read.csv(paste0(d,file),check.names=FALSE,row.names=1) %>% t() %>% data.frame()
    deconv_result=apply(deconv_result,1,function(x) (x-min(x))/sum(x)) %>% t()
    deconv_result=apply(deconv_result,2,function(y){
        y[y>=quantile(y,quantile_cutoff)]=quantile(y,quantile_cutoff)
        return(y)
    }) %>% scale() %>% data.frame()
    
    return(deconv_result)
    
}) %>% dplyr::bind_rows()

# rownames(df)=gsub('.','-',rownames(df))

colnames(df)=gsub('\\.',' ',colnames(df))
df=df[,c('Parachymal cell','Fibroblast','Endothelial cell','CD4 T','CD8 T','B','Plasma cell',
         'Macrophage','cDC')]

rownames(df)=gsub('[.]','-',rownames(df))

df_scaled=df
df_scaled_unbin=df_scaled
df_scaled=scale(df_scaled)
df_scaled[df_scaled<0]=0
df_scaled[df_scaled!=0]=1
df_scaled=data.frame(df_scaled,check.names=FALSE)
# unique(df_scaled)
niche_type=rep(NA,length.out=nrow(df_scaled))
names(niche_type)=rownames(df_scaled)

# niche type assignment
# Niche 1 2 3 4
niche_type[df_scaled %>% filter(`Parachymal cell`==1,rowSums(.)==1) %>% rownames()]='Niche_1-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Endothelial cell`==1,rowSums(.)==2) %>% rownames()]='Niche_2-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Fibroblast`==1,`Endothelial cell`==1,rowSums(.)==3) %>% rownames()]='Niche_2-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Fibroblast`==1,rowSums(.)==2) %>% rownames()]='Niche_3-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Macrophage`==1,rowSums(.)==2) %>% rownames()]='Niche_4-like'
# Niche 7 8 10 12
niche_type[df_scaled %>% filter(`Endothelial cell`==1,rowSums(.)==1) %>% rownames()]='Niche_7-like'
niche_type[df_scaled %>% filter(`Fibroblast`==1,`Endothelial cell`==1,rowSums(.)==2) %>% rownames()]='Niche_7-like'
niche_type[df_scaled %>% filter(`Fibroblast`==1,rowSums(.)==1) %>% rownames()]='Niche_8-like'
niche_type[df_scaled %>% filter(`Macrophage`==1,rowSums(.)==1) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Fibroblast`==1,`Macrophage`==1,rowSums(.)==2) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Endothelial cell`==1,`Macrophage`==1,rowSums(.)==2) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Fibroblast`==1,`Endothelial cell`==1,`Macrophage`==1,rowSums(.)==3) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Fibroblast`==1,`Macrophage`==1,rowSums(.)==3) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Parachymal cell`==1,`Endothelial cell`==1,`Macrophage`==1,rowSums(.)==3) %>% rownames()]='Niche_10-like'
niche_type[df_scaled %>% filter(`Plasma cell`==1) %>% rownames()]='Niche_12-like'
# # Niche 5 6 9 11 13
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% filter(`Parachymal cell`==1) %>% rownames()]='Niche_5/Niche_6-like'
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% 
           filter(  `B`==1 ) %>% rownames()]='Niche_13-like'
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% 
           filter( (`cDC`==1) & (`CD4 T`==1)   ) %>% rownames()]='Niche_13-like'
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% filter(`Macrophage`==1,`Fibroblast`!=1,`Endothelial cell`!=1) %>% rownames()]='Niche_11-like'
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% filter(`Fibroblast`==1) %>% rownames()]='Niche_9-like'
niche_type[df_scaled[niche_type[is.na(niche_type)] %>% names(),] %>% filter(`Endothelial cell`==1) %>% rownames()]='Niche_9-like'

meta.data=lapply(c('KIRC','PRAD','NSCLC'),function(tumor_type){
    
    file_ls=grep(tumor_type,dir(d),value=TRUE)

    deconv_result=lapply(file_ls,function(file){
        paste0(d,file) %>% read.csv(check.names=FALSE,row.names=1)
    }) %>% dplyr::bind_cols()
    
    data.frame('sample'=colnames(deconv_result),'type'=tumor_type)
    
}) %>% dplyr::bind_rows()

response_data=lapply(file_ls,function(file){
    load(file)
    df=data.frame('sample'=colnames(tpm_mat),'response'=response)
    df[,'dataset']=file %>% gsub('[A-Za-z0-9_+]+/processedData/','',.) %>% gsub('\\.[A-Za-z]+$','',.)
    return(df)
}) %>% dplyr::bind_rows()
meta.data=dplyr::left_join(meta.data,response_data)
meta.data[,'sample']=gsub('\\.','-',meta.data[,'sample'])
rownames(meta.data)=meta.data[,'sample']
meta.data[,'niche_type']=niche_type[rownames(meta.data)]
meta.data['seperated']=meta.data['dataset']

# meta analysis
meta.independence=function(a,b,c,d,sm='OR',method='MH',...){
    
    meta_analysis_result=metabin(event.e=a,n.e=a+b,
                                 event.c=c,n.c=c+d,
                                 sm=sm,method=method,...)
    
    if (meta_analysis_result$pval.Q>0.05) {
        return(c(meta_analysis_result$pval.Q,meta_analysis_result$lower.common,
                 1/2*(meta_analysis_result$lower.common+meta_analysis_result$upper.common),
                 meta_analysis_result$upper.common,
                 meta_analysis_result$pval.common))
    } else {
        return(c(meta_analysis_result$pval.Q,meta_analysis_result$lower.random,
                 1/2*(meta_analysis_result$lower.random+meta_analysis_result$upper.random),
                 meta_analysis_result$upper.random,
                 meta_analysis_result$pval.random))
    }
}

meta_df=sapply(setdiff(meta.data.stat$niche_type,'Niche_11-like'),function(y){
    
    df=lapply(unique(meta.data.stat$dataset),function(x){
        meta.data.stat.tmp=meta.data.stat %>% filter(dataset==x)

        if (! y %in% meta.data.stat.tmp$niche_type){
            return(NULL)
        }

        a=meta.data.stat.tmp %>% filter(response=='R',niche_type==y) %>% .[,'Num'] %>% sum()
        b=meta.data.stat.tmp %>% filter(response=='R',niche_type!=y) %>% .[,'Num'] %>% sum()
        c=meta.data.stat.tmp %>% filter(response=='NR',niche_type==y) %>% .[,'Num'] %>% sum()
        d=meta.data.stat.tmp %>% filter(response=='NR',niche_type!=y) %>% .[,'Num'] %>% sum()

        return(data.frame(study=x,a=a,b=b,c=c,d=d))

    }) %>% dplyr::bind_rows()
    
    meta_analysis_result=meta.independence(df$a,df$b,df$c,df$d,
                                           sm='RD',method='MH')
    
    return(c(meta_analysis_result))
    
},USE.NAMES=TRUE) %>% t() %>% as.data.frame()

colnames(meta_df)=c('pval.Q','combine_lower','combine_mean','combine_upper','combine_p')
meta_df[,'niche_type']=rownames(meta_df)

meta_df[,'niche_type']=factor(meta_df[,'niche_type'],levels=rev( c('Niche_1-like','Niche_2-like','Niche_3-like','Niche_4-like','Niche_5/Niche_6-like','Niche_7-like',
                                                                   'Niche_8-like','Niche_9-like','Niche_10-like','Niche_11-like','Niche_12-like',
                                                                   'Niche_13-like')) )
meta_df=meta_df %>% arrange(niche_type)
meta_df[,'y']=1:nrow(meta_df)

y='Niche_4-like'
df=lapply(unique(meta.data.stat$dataset),function(x){
    meta.data.stat.tmp=meta.data.stat %>% filter(dataset==x)

    if (! y %in% meta.data.stat.tmp$niche_type){
        return(NULL)
    }

    a=meta.data.stat.tmp %>% filter(response=='R',niche_type==y) %>% .[,'Num'] %>% sum()
    b=meta.data.stat.tmp %>% filter(response=='R',niche_type!=y) %>% .[,'Num'] %>% sum()
    c=meta.data.stat.tmp %>% filter(response=='NR',niche_type==y) %>% .[,'Num'] %>% sum()
    d=meta.data.stat.tmp %>% filter(response=='NR',niche_type!=y) %>% .[,'Num'] %>% sum()

    return(data.frame(study=x,a=a,b=b,c=c,d=d))

}) %>% dplyr::bind_rows()

meta_analysis_result=metabin(event.e=df$a,n.e=df$a+df$b,
                             event.c=df$c,n.c=df$c+df$d,
                             studylab=df$study,
                             sm='ASD',method='MH')

options(repr.plot.width=10,repr.plot.height=3)
p=forest(meta_analysis_result,random=FALSE,file='./Niche_4-like_meta.pdf',func.gr=pdf,width=10,height=3)