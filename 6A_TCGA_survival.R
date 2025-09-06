library(dplyr)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(DWLS)

# DWLS
for (cancer_type in c('BRCA','COAD','HNSC','LUAD','LUSC','PRAD','OV','PAAD','KIRC','LIHC''KIRC','HNSC')){
    seurat_obj=readRDS(paste0('../2_deconvolution/scRNA_used/',cancer_type,'_new.rds'))

    cell_type_label=Idents(seurat_obj) %>% as.character()
    cell_type_label[which(cell_type_label %in% tumor_cell_type[[cancer_type]])]='Parachymal cell'
    cell_type_label=sub(' ','.',cell_type_label)

    sc_sig=buildSignatureMatrixUsingSeurat(scdata=seurat_obj@assays$RNA@data,
                                           id=cell_type_label,
                                           path=path.dir)
    
    bulk_RNA_mat=read.table(paste0('./TCGA_TPM/TCGA-',cancer_type,'_mRNA_TPM.txt'))

    deconv_result=apply(bulk_RNA_mat,2,function(x){
      tr=trimData(sc_sig,x)
      solveDampenedWLS(tr$sig, tr$bulk)
    })

}

meta.data=read.csv('./TCGA_clinical_data.csv',row.names=1,check.names=FALSE)
meta.data=data.frame(meta.data,row.names=1)

quantile_cutoff=1

deconv_result=lapply(c('BRCA','COAD','LUAD','HNSC','LUSC' ,'PRAD','OV',
             'PAAD','KIRC','LIHC'),function(tumor_type){ 
    x=read.csv(paste0('./TCGA_deconv_DWLS/',tumor_type,'.csv'),row.names=1,check.names=FALSE)
    x=t(x)

    x=apply(x,2,function(y){
        y[y>=quantile(y,quantile_cutoff)]=quantile(y,quantile_cutoff)
        return(y)
    })
    x=apply(x,1,function(y) {(y-min(y))/sum(y)}) %>% t() %>% as.data.frame()
    x=scale(x)
    x=data.frame(x,check.names=FALSE)
    return(x)
}) 

df=dplyr::bind_rows(deconv_result)
colnames(df)=gsub('\\.',' ',colnames(df))
df=df[,c('Parachymal cell','Fibroblast','Endothelial cell','CD4 T','CD8 T','B','Plasma cell',
         'Macrophage','cDC')]

rownames(df)=gsub('[.]','-',rownames(df))

selected_samples=strsplit(rownames(df),split='-') %>% sapply(function(x) {grepl('01',x[4])})
df=df[selected_samples,]

df_scaled=df
df_scaled[,'patient']=gsub('-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+$','',rownames(df_scaled))
df_scaled[,'sample']=rownames(df_scaled)
df_scaled[,'type']= meta.data[ df_scaled[,'patient'] ,'type']
df_scaled=df_scaled %>% group_by(type) %>% mutate_if(is.numeric,scale) %>% ungroup() %>% data.frame(check.names=FALSE)
rownames(df_scaled)=df_scaled[,'sample']
df_scaled=df_scaled[,c('Parachymal cell','Fibroblast','Endothelial cell','CD4 T','CD8 T','B','Plasma cell','Macrophage','cDC')]

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

# survival data
sample_info=data.frame(niche_type)
sample_info[,'patient']=gsub('-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+$','',rownames(sample_info))
sample_info=cbind(sample_info,meta.data[sample_info[,'patient'],])
cutoff.time=1000
sample_info_tmp[ sample_info_tmp %>% filter(OS.time>=cutoff.time,! is.na(OS.time)) %>% rownames() ,'OS']=0
sample_info_tmp[ sample_info_tmp %>% filter(OS.time>=cutoff.time,! is.na(OS.time)) %>% rownames() ,'OS.time']=cutoff.time

# adjusted survival curves
ggadjustedcurves(cox.model,
                 data=sample_info_tmp,
                 variable='niche_type',method='marginal',legend='right',
                 ylim=c(0.43,1),
                 size=0.3)

# comparsison of each pair
niche_test_result=lapply(unique(niche_type),function(x){
    lapply(niche_type,function(y){
        
        if (x!=y){
            
            niche=c(x,y)

            sample_info_tmp_niche=sample_info_tmp
            sample_info_tmp_niche[,'niche_type']=as.character(sample_info_tmp_niche[,'niche_type'])

            sample_info_tmp_niche=sample_info_tmp_niche %>% filter(niche_type %in% niche)

            cox.model.niche=coxph(Surv(OS.time,OS)~type+niche_type,
                            data=sample_info_tmp_niche,x=TRUE) %>% summary()

            return(data.frame(x,y,p=cox.model.niche$coefficients %>% .[nrow(.),ncol(.)]))
             
        } else {
            return(data.frame(x,y,p=2))
        }
        
    })
})  %>% unlist(recursive=FALSE) %>% .[!sapply(.,is.null)] %>% dplyr::bind_rows()

# survival curves of each tumor type
PlotSurvivalOneVsAll=function(surv_data,group_col='niche_type',
                              surv_time_col='OS.time',
                              surv_censoring_col='OS',
                              xlim=NULL,
                              groups=NULL,
                              cols=NULL){
    if (is.null(xlim)){
        xlim=c(0,max(surv_data[,surv_time_col],na.rm=TRUE))
    }
    if (is.null(groups)){
        groups=unique(surv_data[,group_col])
    }
    if (is.null(cols)){
        cols=c('#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
               '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
               '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
               '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
               '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF')[1:length(groups)]
    }
    
    cols=cols[groups %in% surv_data[,group_col]]
    groups=groups[groups %in% surv_data[,group_col]]
    
    f=as.formula(paste0('Surv(',surv_time_col,',',surv_censoring_col,')','~1'))
    
    for (i in 1:length(groups)){
        surv_data_subset=surv_data[ surv_data[,group_col]==groups[i] , ]
        km.model.group=survfit(f,data=surv_data_subset)
        plot(km.model.group,col=cols[i],conf.int=FALSE,mark.time=TRUE,xlim=xlim,ylim=c(0,1),las=1,
             xaxt='n',yaxt='n',bty='n',xlab='',ylab='')
        par(new=TRUE)
    }
    km.model.all=survfit(f,data=surv_data)
    plot(km.model.all,col='black',conf.int=FALSE,mark.time=TRUE,xlim=xlim,ylim=c(0,1),las=1)
}

tumor_type='HNSC'
sample_info_tmp_tumor=sample_info_tmp
sample_info_tmp_tumor=sample_info_tmp_tumor %>% filter(type==tumor_type)
surv_data=sample_info_tmp_tumor
km.model=survfit(Surv(OS.time,OS)~niche_type,data=sample_info_tmp_tumor)

PlotSurvivalOneVsAll(surv_data,
                     groups=g,cols=c,
                     xlim=c(0,1800))
