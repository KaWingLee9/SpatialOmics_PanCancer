# in-house pipeline to identify tissue compartments
CNVestimate_ST <- function(seurat_obj,deconv_result,MalType,
                           StrimmType=c('Fibroblast','Endothelial cell','CD4 T','CD8 T','B','Plasma cell',
                           'NK','Macrophage','cDC','pDC','Mast cell'),
                           out_dir='./', # output for saving files
                           ref_cell_ratio=0.1,# minimal percents of spots considered to be reference in inferCNV or CopyKAT
                           method='infercnv',runCNV_only=FALSE,
                           gene_ordering_table='',num_threads=6, # parameters for inferCNV
                           genome='hg20', # parameters for CopyKAT
                           cluster_method='exp',scale.max=10,npcs=50,k.param=20,dims=1:10,resolution=0.8, # parameters for clustering
                           exponent=3, cv_thres=0.1,# denoise parameters
                           StrImm_ratio=0.4,Para_cor=0.7, # parameters for benign/malignant identification
                           heatmap_cutoff=0.01, # pararmeters for heatmap visualization
                           marker=NULL,
                           ...){
    
    ### Key parameters of CNVestimate_ST:
    # seurat_obj: the SeuratV4 object
    # deconv_result: a data frame of the cell type fraction (the row number should be same as the column number of seurat object)
    # MalType: the malignant cell type (should be in the colnames of deconv_result)
    # StrimmType: the stromal/immune cell type (should be in the colnames of deconv_result)
    # out_dir: directory to store the output
    # gene_ordering_table: the data frame passed to gene_order_file in infercnv::CreateInfercnvObject
    # method: one of infercnv or infercnv_HMM
    # cluster_method: methods for cell type clustering, could be exp or CNV
    # StrImm_ratio: least propotion to consider spots with Parachymal parts
    # cv_thres, exponent: the lower bound of CV and the exponent of CV in the process of CNV denoising


    # threshold for benign or malignant tissue
    # least propotion to consider spots with Parachymal parts
    
    # Visualization
    # q: quantile for heatmap expression
    
    library(dplyr)
    library(Seurat)
    library(ComplexHeatmap)
    library(ggplot2)

    StrimmType=intersect(StrimmType,colnames(deconv_result))
    
    # library loading and create directory for output
    if (method=='infercnv'){
        library(infercnv)
        save_dir=file.path(out_dir,'infercnv')
        if (!dir.exists(save_dir)){
            dir.create(save_dir)
        }else{
            paste(save_dir,list.files(save_dir),sep="/") %>% file.remove()
        }
    } else if(method=='infercnv_HMM'){
        library(infercnv)
        save_dir=file.path(out_dir,'infercnv_HMM')
        if (!dir.exists(save_dir)){
            dir.create(save_dir)
        }else{
            paste(save_dir,list.files(save_dir),sep="/") %>% file.remove()
        }
    }else if (method=='copykat'){
        library(copykat)
        save_dir=file.path(out_dir,'copykat')
        if (!dir.exists(save_dir)){
            dir.create(save_dir)
        }else{
            paste(save_dir,list.files(save_dir),sep="/") %>% file.remove()
        }
    }else{
        stop(paste0('No such CNV estimation method: ',method))
    }
    
    # count matrix
    # if (class(seurat_obj))
    # if (class(deconv_result))
    # if (dim(seurat_obj)[1] != nrow(deconv_result))
    raw_count_mat <- as.matrix(seurat_obj@assays[["RNA"]]@counts)
    
    # write.csv(deconv_result,file.path(save_dir,'deconv_result.csv'))
    
    # extract spots with least malignant component as reference for CNV estimation
    NonMal_rate=deconv_result %>% select(-all_of(MalType)) %>% apply(1,sum)
    # NonMal_rate=deconv_result %>% select(all_of(StrimmType)) %>% apply(1,sum)
    NonMal_thres=quantile(NonMal_rate,probs=1-ref_cell_ratio)
    
    ref_spot=NonMal_rate<=NonMal_thres
    ref_spot=ref_spot+0
    cell_type_anno=factor(ref_spot,levels=c(0,1),labels=c('Ref_spot','Obs_spot')) %>% data.frame('Ref/Obs'=.,check.names=FALSE)
    # write.csv(NonMal_rate %>% data.frame('NonMal_rate'=.,check.names=FALSE),file.path(save_dir,'NonMal_rate.csv'))
    
    # remove sample with too less non-malignant components
    if (NonMal_thres<=1-StrImm_ratio){
        warning('There are too low non-malignant component in the spots')
        spot_df=data.frame('NonMal_rate'=NonMal_rate,CNVlabel='Malignant')
        rownames(spot_df)=colnames(seurat_obj)
        write.csv(spot_df,file.path(save_dir,'spot_result.csv'))
        return('')
    }
    
    # CNV estimation using infercnv, infercnv_HMM(not recommended), copykat, etc.
    if (method=='infercnv'){
        
        infercnv_obj=CreateInfercnvObject(raw_counts_matrix=raw_count_mat,
                                          annotations_file=cell_type_anno,
                                          gene_order_file=gene_ordering_table,
                                          ref_group_names='Ref_spot',
                                          min_max_counts_per_cell=c(-Inf,Inf),
                                          chr_exclude=c("chrX","chrY","chrM","KI270734.1"))
        infercnv_obj=infercnv::run(infercnv_obj,
                                   cutoff=0.1,
                                   cluster_by_groups=FALSE,
                                   sd_amplifier=1.5,
                                   out_dir=save_dir,
                                   denoise=TRUE,
                                   save_rds=FALSE,
                                   save_final_rds=TRUE,
                                   HMM=FALSE)
        cnv_mat=infercnv_obj@expr.data %>% t() %>% data.frame(check.names=FALSE)
    }
    
    if (method=='infercnv_HMM'){
        infercnv_obj=CreateInfercnvObject(raw_counts_matrix=raw_count_mat,
                                          annotations_file=cell_type_anno,
                                          gene_order_file=gene_ordering_table,
                                          ref_group_names='Ref_spot',
                                          min_max_counts_per_cell=c(-Inf,Inf),
                                          chr_exclude=c("chrX","chrY","chrM","KI270734.1"))
        infercnv_obj=infercnv::run(infercnv_obj,
                                   cutoff = 0.1,
                                   out_dir = save_dir,
                                   cluster_by_groups=TRUE,
                                   denoise=TRUE,
                                   HMM=TRUE,
                                   HMM_type="i6",
                                   analysis_mode='subcluster',
                                   BayesMaxPNormal=0,
                                   num_threads=num_threads,
                                   save_rds=FALSE,
                                   save_final_rds=TRUE,
                                   no_plot=FALSE)
        cnv_mat=infercnv_obj@expr.data %>% t() %>% data.frame(check.names=FALSE)
    }
    
    
    if (method=='copykat'){
        
        copykat_obj=copykat(rawmat=raw_count_mat,
                     id.type="S", # 'S': gene symbol; 'E': Ensembl ID
                     ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", 
                     distance="euclidean", 
                     norm.cell.names='Ref_spot',
                     output.seg="FLASE", 
                     plot.genes="TRUE", 
                     genome="hg20", # 'hg20' for human; 'mm10' for mouse
                     n.cores=10)
        # cnv_mat=
        # cnv_mat=exp(cnv_mat)
    }
    
    if (runCNV_only){
        return('')
    }

    # use provided cluster 
    if (is.null(cluster_method)){
        cluster_result=data.frame(cluster=Idents(seurat_obj))
    }

    if (cluster_method=='exp'){
        # remove cell cycling genes when clustering
        # gene_removed <- setdiff(rownames(seurat_obj),unlist(cc.genes.updated.2019),c(grep("^RP[SL]",rownames(seurat_obj),value=T)))
        # seurat_obj <- seurat_obj[gene_removed,]
        
        seurat_obj=seurat_obj %>% NormalizeData(verbose=FALSE) %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000,verbose=FALSE) %>%
            ScaleData(verbose=FALSE) %>%
            RunPCA(npcs=npcs,slot='data',verbose=FALSE) %>%
            FindNeighbors(dims=dims,k.param=k.param,annoy.metric="euclidean",verbose=FALSE) %>%
            FindClusters(algorithm = 1,resolution=resolution,verbose=FALSE)
        cluster_result=data.frame(cluster=Idents(seurat_obj))
    }
    
    # cluster spots using CNV profile in Seurat style
    if (cluster_method=='CNV'){
        # cluster_result=ClusterCNV(cnv_mat, 
        #                           scale.max=scale.max,
        #                           npcs=npcs,
        #                           k.param=k.param,
        #                           dims=dims,
        #                          resolution=resolution)
        cluster_result=ClusterCNV3(cnv_mat,scale.max=scale.max,
                                  npcs=npcs,
                                 k.param=k.param,
                                  dims=dims,
                                 resolution=resolution)
    }

    # denoise CNV profile using CV (coefficient of variation) of each cluster (should be parameterized)
    cnv_mat_denoise=DenoiseCNV(cnv_mat,
                               cluster_result,
                               cv_thres=cv_thres,
                               q=0.005,
                               p=exponent)
    
    # remove HLA genes
    # gene_removed <- setdiff(colnames(cnv_mat_denoise),c(grep('^HLA-',rownames(seurat_obj),value=T)))
    # cnv_mat_denoise <- cnv_mat_denoise[,gene_removed]
    
    # prepare for Benign/Malignant identification
    # cnv_score=apply(cnv_mat_denoise,1,function(x){sum(x)^2}) (Wrong!!!!!)
    cnv_score=apply(cnv_mat_denoise,1,function(x){sum(x^2)})
    # cnv_score=apply(cnv_mat_denoise,1,function(x){sum(abs(x))})
    NonMal_rate=deconv_result %>% select(all_of(StrimmType)) %>% apply(1,sum)
    spot_df=data.frame('Stromal/Immune propotion'=NonMal_rate,
                        'Ref/Obs'=cell_type_anno,
                        Cluster=cluster_result$cluster,
                        CNVscore=cnv_score,
                        check.names=FALSE)
    rownames(spot_df)=rownames(cluster_result)



    # benign/malignant prediction
    # output_ls=PredictMalignant(cnv_mat_denoise,
    #                            spot_df,
    #                            StrImm_ratio=StrImm_ratio,
    #                            Para_cor=Para_cor)
    output_ls=PredictMalignant_V2(cnv_mat_denoise,
                                spot_df)
    # return(output_ls)
    cluster_df=output_ls$cluster_df
    cluster_cnv=output_ls$cluster_cnv
    fct_label=output_ls$fct_label
    spot_df$CNVlabel=factor(spot_df$Cluster,labels=fct_label)

    write.table(t(cnv_mat_denoise),file=file.path(save_dir,'cnv.profile.denoise.txt'),quote=FALSE)
    write.csv(spot_df,file.path(save_dir,paste0('spot_result_res',resolution,'.csv')))
    write.csv(cluster_df,file.path(save_dir,paste0('cluster_result_res',resolution,'.csv')))
    
    
    # print figure
    pdf(file.path(save_dir,'infercnv_result_clustered.pdf'))
    par(mfrow=c(1,1))
    
    # show Stromal and Immune propotion spatially
    seurat_obj@meta.data$'Stromal and Immune ratio'=spot_df[,'Stromal/Immune propotion']
    print(SpatialFeaturePlot(seurat_obj,'Stromal and Immune ratio',image.alpha=0))
    # show Ref_spot and Obs_spot spatially
    Idents(seurat_obj)=cell_type_anno
    print(SpatialDimPlot(seurat_obj,image.alpha=0)+guides(fill=guide_legend(title='Reference/Observation spot')))

    # Heatmap with CNV-clustered data
    col_cnv=circlize::colorRamp2(unique(c(seq(quantile(unlist(cnv_mat_denoise),heatmap_cutoff),0,length.out=6),
                                        0,
                                        seq(0,quantile(unlist(cnv_mat_denoise),1-heatmap_cutoff),length.out=6))),
                                 c('#104680','#317CB7','#6DADD1','#B6D7E8','#E9F1F4','white',
                                   '#FBE3D5','#F6B293','#DC6D57','#B72230','#6D011F'))

    col_chr=rep(c('#B5B5B5','#1C1C1C'),11)
    names(col_chr)=1:22

    col_type=c('#5050FFFF','#CE3D32FF','#749B58FF','#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF','#924822FF',
      '#837B8DFF','#C75127FF','#D58F5CFF','#7A65A5FF','#E4AF69FF','#3B1B53FF','#CDDEB7FF','#612A79FF','#AE1F63FF','#E7C76FFF','#5A655EFF',
      '#CC9900FF','#99CC00FF','#A9A9A9FF','#CC9900FF','#99CC00FF','#33CC00FF','#00CC33FF','#00CC99FF','#0099CCFF','#0A47FFFF','#4775FFFF',
      '#FFC20AFF','#FFD147FF','#990033FF','#991A00FF','#996600FF','#809900FF','#339900FF','#00991AFF','#009966FF','#008099FF','#003399FF',
      '#1A0099FF','#660099FF','#990080FF','#D60047FF','#FF1463FF','#00D68FFF','#14FFB1FF')[1:nrow(cluster_cnv)]
    names(col_type)=rownames(cluster_cnv)

    Heatmap(cnv_mat_denoise %>% arrange(spot_df[,'Cluster']),name='CNV\nprofile',col=col_cnv,
            cluster_rows=FALSE,cluster_columns=FALSE,
            show_row_names=FALSE,show_column_names=FALSE,
            bottom_annotation=HeatmapAnnotation(Chr=as.numeric(infercnv_obj@gene_order[,'chr']),
                                                col=list(Chr=col_chr),show_legend=FALSE,annotation_name_side='left'),
            right_annotation=rowAnnotation('Cell\ntype'=arrange(spot_df,Cluster) %>% .[,'Cluster'],col=list('Cell\ntype'=col_type),
                                           show_annotation_name=FALSE)) %>% 
    print()
    
   Heatmap(cluster_cnv,name='CNV\nprofile',col=col_cnv,cluster_columns=FALSE,# clustering_distance_rows='spearman',
            show_column_names=FALSE,
            bottom_annotation=HeatmapAnnotation(Chr=as.numeric(infercnv_obj@gene_order[,'chr']),
                                                col=list(Chr=col_chr),show_legend=FALSE,annotation_name_side='left'),
            right_annotation=rowAnnotation('Cell\ntype'=rownames(cluster_cnv),col=list('Cell\ntype'=col_type),
                                           show_annotation_name=FALSE,show_legend=FALSE)) %>% 
    print() 
    
    options(scipen=7)

    # show cluster spatially
    Idents(seurat_obj)=spot_df$Cluster
    print(SpatialDimPlot(seurat_obj,cols=col_type,image.alpha=0)+guides(fill=guide_legend(title='Cluster')))
    
    # Line plot of mean CNV profile
    df=cluster_cnv %>% as.matrix
    colnames(df)=1:ncol(cluster_cnv)
    df=reshape2::melt(df,varnames=c('Cluster','Gene'),value.name='CNV')
    df$Cluster=as.factor(df$Cluster)

    print(ggplot(df,aes(x=Gene,y=CNV,color=Cluster))+
        geom_line(alpha=0.6,linewidth=1.1)+
        scale_color_manual(values=col_type)+
        geom_hline(yintercept=0,linetype='dashed',linewidth=1.2)+
        xlim(c(1,ncol(cluster_cnv)))+
        # ylim(c(-1,1))+
        theme_classic()+
        theme(axis.line.x=element_blank(),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()))

    # CNVscore and Stromal/Immune percent distribution
    
    print(ggplot(data=spot_df,aes(x=Cluster,y=CNVscore))+
        geom_boxplot()+
        theme_classic()+
        xlab('Cluster')+
        ylab('CNVscore'))

    print(ggplot(data=spot_df,aes(x=Cluster,y=`Stromal/Immune propotion`))+
        geom_boxplot()+
        theme_classic()+
        xlab('Cluster')+
        ylab('Stromal and Immune propotion/%'))

    seurat_obj@meta.data$'CNVscore'=spot_df[,'CNVscore']
    print(SpatialFeaturePlot(seurat_obj,'CNVscore',image.alpha=0))
    
    Heatmap(cor(t(cluster_cnv),method='pearson'),name='Correlation') %>% print()
    Heatmap(cor(t(cluster_cnv),method='spearman'),name='Correlation') %>% print()
    
    Idents(seurat_obj)=spot_df[,'CNVlabel']
    col_tissue=c('Stromal'='#501D8A','Benign'='#1C8041','Malignant'='#E55709')
    col_tissue=col_tissue[intersect(names(col_tissue),unique(Idents(seurat_obj)))]
    print(SpatialDimPlot(seurat_obj,cols=col_tissue,image.alpha=0)+guides(fill=guide_legend(title='Annotation')))

    if (! is.null(marker)){
        marker=intersect(marker,rownames(seurat_obj))
        if (! is.null(marker)){
            print(SpatialFeaturePlot(seurat_obj,marker,image.alpha=0))
    }}
    
    print(ggplot(data=spot_df,aes(x=CNVlabel,y=CNVscore))+
        geom_boxplot()+
        theme_classic()+
        xlab('CNVlabel')+
        ylab('CNVscore'))

    print(ggplot(data=spot_df,aes(x=`Stromal/Immune propotion`))+
        geom_line(stat='density')+
        theme_classic()+
        xlab('Stromal and Immune propotion/%')+
        ylab('Density'))

    print(ggplot(data=spot_df,aes(x=CNVscore))+
    geom_line(stat='density')+
    theme_classic()+
    xlab('CNVscore')+
    ylab('Density'))
    
    dev.off()
    
    return()
}

# cluster spots using CNV profile in Seurat style
ClusterCNV <- function(cnv_mat,
                       scale.max=scale.max,
                       npcs=npcs,
                       k.param=k.param,
                       dims=dims,
                       resolution=0.5){
    cnv_mat_scaled=ScaleData(cnv_mat,scale.max=scale.max,verbose=FALSE) %>% t()
    cell.embeddings=RunPCA(cnv_mat_scaled,
                           npcs=npcs,
                           verbose=FALSE) %>% Embeddings()
    g_SNN=FindNeighbors(cell.embeddings,k.param=k.param,dims=dims,verbose=FALSE)
    cluster_result=Seurat:::FindClusters(g_SNN[['snn']],resolution=resolution,verbose=FALSE)
    colnames(cluster_result)='cluster'
    cluster_result$cluster=as.factor(cluster_result$cluster)
    return(cluster_result)
}

# cluster spots using CNV profile in Seurat style
ClusterCNV3 <- function(cnv_mat,
                       scale.max=scale.max,
                       npcs=npcs,
                       k.param=k.param,
                       dims=dims,
                       resolution=0.5){
    
    # cluster CNV expression matrix in Seurat style
    cnv_mat_scaled=ScaleData(cnv_mat,scale.max=scale.max,verbose=FALSE) %>% t()
    cell.embeddings=RunPCA(cnv_mat_scaled,
                           npcs=npcs,
                           verbose=FALSE) %>% Embeddings()
    g_SNN=FindNeighbors(cell.embeddings,k.param=k.param,dims=dims,verbose=FALSE)
    cluster_result=Seurat:::FindClusters(g_SNN[['snn']],resolution=resolution,verbose=FALSE)

    dist_mat=parallelDist::parallelDist(cell.embeddings,threads =20)


    # dist_mat=parallelDist::parallelDist(as.matrix(cnv_mat),threads =20)
    x=hclust(dist_mat,method='ward.D2') %>% cutree(k=3)
    cluster_result=data.frame(cluster=x)

    # x=NbClust(t(cell.embeddings),distance=NULL,diss=dist_mat,min.nc=2,max.nc=8,index='silhouette',method='ward.D2')
    # cluster_result=data.frame(cluster=x$Best.partition)

    colnames(cluster_result)='cluster'
    cluster_result$cluster=as.factor(cluster_result$cluster)
    return(cluster_result)
}

# denoise CNV profile using coefficient of variation of each cluster
DenoiseCNV <- function(cnv_mat,
                       cluster_result,
                       cv_thres=0.1,
                       q=0.005,
                       p=3){
    cnv_mat_denoise=log2(cnv_mat) %>% data.frame(check.names=FALSE)
    # cnv_mat_denoise=cnv_mat-1
    # cnv_mat_denoise=data.frame(cnv_mat_denoise,check.names=FALSE)
    if (p>0){
        cv_mat=group_by(cnv_mat_denoise,cluster_result[,'cluster']) %>% summarize_if(is.numeric,function(x){sd(x)/abs(mean(x))}) %>% data.frame(row.names=1)
        cv_mat[cv_mat<=cv_thres]=cv_thres
        cv_mat_denoise=cv_mat[cluster_result[,'cluster'],]
        cnv_mat_denoise=cnv_mat_denoise/(cv_mat_denoise^p)

        q1=quantile(unlist(cnv_mat_denoise),1-q)
        q2=quantile(unlist(cnv_mat_denoise),q)
        cnv_mat_denoise[cnv_mat_denoise>=q1]=q1
        cnv_mat_denoise[cnv_mat_denoise<=q2]=q2

        # cnv_mat_denoise[cnv_mat_denoise>=1]=1
        # cnv_mat_denoise[cnv_mat_denoise<=-1]=-1
    }

    return(cnv_mat_denoise)
}

# benign/malignant prediction
PredictMalignant_V1 <- function(cnv_mat_denoise,
                             spot_df,
                             StrImm_ratio=0.25,
                             Para_cor=0.7){
    cluster_cnv=group_by(cnv_mat_denoise,spot_df[,'Cluster']) %>% summarise_if(is.numeric,mean) %>% data.frame(row.names=1)
    cnv_score=group_by(spot_df,Cluster) %>% summarize(CNVscore=mean(CNVscore)) %>% data.frame(row.names=1)

    # return(cnv_score)

    # Stromal/Parachymal identification
    NonMal_rate=group_by(spot_df,Cluster) %>% summarize('Stromal/Immune propotion'=mean(`Stromal/Immune propotion`)) %>% 
        data.frame(row.names=1,check.names=FALSE)
    x=spot_df[,'Stromal/Immune propotion']
    str_thres=EBImage::otsu(array(x,dim=c(length(x),1)))
    y=filter(NonMal_rate,`Stromal/Immune propotion`<=str_thres) %>% rownames()
    
    # malignant reference cluster by selecting cluster with highest CNVscore in parachymal cluster
    mal_label=filter(cnv_score[y,,drop=FALSE],CNVscore==max(cnv_score[y,,drop=FALSE])) %>% rownames()
    cor_score=cor(t(cluster_cnv),method='pearson') %>% data.frame(check.names=FALSE) %>% .[,mal_label]
    # return(cnv_score)
    # colnames(cor_score)='CORscore'

    cluster_df=data.frame(CNVscore=cnv_score,CORscore=cor_score,
                          NonMal_rate,check.names=FALSE)
    cluster_df[,'Cluster']=rownames(cluster_df)

    # malignant recognized as clusters highly correlated with reference cluster
    ben_label=filter(cluster_df,`Stromal/Immune propotion`<=str_thres,CORscore<=Para_cor) %>% rownames()
    str_label=filter(cluster_df,`Stromal/Immune propotion`>=str_thres) %>% rownames()
    mal_label=setdiff(cluster_df[,'Cluster'],c(ben_label,str_label))
    x=rep('Benign',length(ben_label))
    names(x)=ben_label
    y=rep('Malignant',length(mal_label))
    names(y)=mal_label
    z=rep('Stromal',length(str_label))
    names(z)=str_label
    fct_label=c(x,y,z) %>% .[cluster_df$Cluster]
    
    cluster_df$CNVlabel=factor(cluster_df$Cluster,labels=fct_label)
    
    output=list(cluster_df=cluster_df,cluster_cnv=cluster_cnv,fct_label=fct_label)
    return(output)
}


# benign/malignant prediction
PredictMalignant_V2 <- function(cnv_mat_denoise,spot_df){
    cluster_cnv=group_by(cnv_mat_denoise,spot_df[,'Cluster']) %>% summarise_if(is.numeric,mean) %>% data.frame(row.names=1)
    cnv_score=group_by(spot_df,Cluster) %>% summarize(CNVscore=mean(CNVscore)) %>% data.frame(row.names=1)

    # Stromal/Parachymal identification
    NonMal_rate=group_by(spot_df,Cluster) %>% summarize('Stromal/Immune propotion'=mean(`Stromal/Immune propotion`)) %>% 
        data.frame(row.names=1,check.names=FALSE)
    str_thres=EBImage::otsu(array(spot_df[,'Stromal/Immune propotion'],
                            dim=c(length(spot_df[,'Stromal/Immune propotion']),1)))

    str_label=filter(NonMal_rate,`Stromal/Immune propotion`>=str_thres) %>% rownames()

    # malignant reference identification
    # y=filter(NonMal_rate,`Stromal/Immune propotion`<=str_thres) %>% rownames()
    # mal_label_ref=filter(cnv_score,CNVscore==max(cnv_score)) %>% rownames()
    # mal_label_ref=filter(cnv_score,CNVscore==max(cnv_score)) %>% rownames()
    cluster_underestimate=setdiff(rownames(cnv_score),str_label) 

    if (length(cluster_underestimate)!=0) {
        cnv_score_ns=cnv_score[ setdiff(rownames(cnv_score),str_label) ,,drop=FALSE] %>% data.frame()
        mal_label_ref=filter(cnv_score_ns,CNVscore==max(cnv_score_ns)) %>% rownames()

        if (length(str_label)>1){
            c=cor(t(cluster_cnv[str_label,,drop=FALSE]),t(cluster_cnv[mal_label_ref,,drop=FALSE]))
            str_label_ref=which(c[,1]==min(c)) %>% names()
        }else{
            str_label_ref=str_label
        }
        # str_label_ref=filter(NonMal_rate,`Stromal/Immune propotion`==max(NonMal_rate)) %>% rownames()
        cluster_underestimate=setdiff(rownames(cnv_score),c(mal_label_ref,str_label))

        if (length(cluster_underestimate)==0){
            x=NULL
        } else {
            x=sapply(cluster_underestimate,function(i){

                # if (length(str_label)>=1){
                #     str_cnv=apply(cluster_cnv_denoise[str_label,],2,mean)
                # }
                # if (length(str_label)==1){
                #     str_cnv=unlist(cluster_cnv_denoise[str_label,])
                # }

                if (length(str_label)==0){
                    str_cnv=cnv_mat_denoise[filter(spot_df,`Ref/Obs`=='Ref_spot') %>% rownames(),] %>% apply(2,mean)
                }else{
                    str_cnv=unlist(cluster_cnv[str_label_ref,])
                }

                estimate_cnv=unlist(cluster_cnv[i,])
                mal_cnv=unlist(cluster_cnv[mal_label_ref,])
                # secondary denoising   
                # estimate_cnv=estimate_cnv/NonMal_rate[i,'Stromal/Immune propotion']
                # mal_cnv=mal_cnv/NonMal_rate[mal_label_ref,'Stromal/Immune propotion']
                estimate_cnv=estimate_cnv*((1-NonMal_rate[mal_label_ref,'Stromal/Immune propotion'])/(1-NonMal_rate[i,'Stromal/Immune propotion']))

                # distance comparison between Stromal_ref and Malignant_ref
                m=dist(rbind(estimate_cnv,str_cnv))
                n=dist(rbind(estimate_cnv,mal_cnv))
                if (m>=n){
                    return('Malignant')
                }
                if (m<n){
                    return('Benign')
                }

                # correlation comparison between Stromal_ref and Malignant_ref
                # m=cor(estimate_cnv,str_cnv,method='spearman')
                # n=cor(estimate_cnv,mal_cnv,method='spearman')
                # if (m<n){
                #     return('Malignant')
                # }
                # if (m>n){
                #     return('Benign')
                # }
                
            },USE.NAMES=TRUE)
        }

    } else {
        x=NULL
        mal_label_ref=NULL
    }




    
    # reassign clusters
    # if ('Benign' %in% x){
    #     ben_label_ref=which(x=='Benign') %>% names()
    #     if (length(ben_label_ref)>=1){
    #         ben_cnv=apply(cluster_cnv_denoise[ben_label_ref,],2,mean)
    #     }
    #     if (length(ben_label_ref)==1){
    #         ben_cnv=unlist(cluster_cnv_denoise[ben_label_ref,])
    #     }
    #     cluster_underestimate=setdiff(rownames(cnv_score),c(mal_label_ref,str_label,ben_label_ref))
    #     for (i in cluster_underestimate){
    #         m=dist(rbind(unlist(cluster_cnv_denoise[i,]),unlist(ben_cnv)))
    #         n=dist(rbind(unlist(cluster_cnv_denoise[i,]),unlist(cluster_cnv_denoise[mal_label_ref,])))
    #         if (m<n){
    #             x[i]='Benign'
    #         }
    #     }
    # }

    #     
    y=rep('Malignant',length(mal_label_ref))
    names(y)=mal_label_ref
    z=rep('Stromal',length(str_label))
    names(z)=str_label

    cluster_df=data.frame(CNVscore=cnv_score,
                          NonMal_rate,check.names=FALSE)
    cluster_df[,'Cluster']=rownames(cluster_df)

    fct_label=c(x,y,z) %>% .[cluster_df$Cluster]

    cluster_df$CNVlabel=factor(cluster_df$Cluster,labels=fct_label)

    output=list(cluster_df=cluster_df,cluster_cnv=cluster_cnv,fct_label=fct_label)
    return(output)
}



Other_information <- function(x){
    x=unlist(cnv_mat_denoise)
    x1=x[x>=0]
    x2=x[x<=0]
    thres1=EBImage::otsu(array(x1,dim=c(length(x1),1)))
    thres2=-EBImage::otsu(array(abs(x2),dim=c(length(x2),1)))
    thres1
    thres2

    cnv_mat_denoise_bin=cnv_mat_denoise
    cnv_mat_denoise_bin[cnv_mat_denoise_bin>=0 & cnv_mat_denoise_bin <= thres1]=0
    cnv_mat_denoise_bin[cnv_mat_denoise_bin<=0 & cnv_mat_denoise_bin >= thres2]=0
    # cnv_mat_denoise_bin[cnv_mat_denoise_bin>=1]=1
    # cnv_mat_denoise_bin[cnv_mat_denoise_bin<=-1]=-1
    cnv_mat_denoise_bin[cnv_mat_denoise_bin!=0]=1

    # cnv_mat_denoise_bin[cnv_mat_denoise_bin<=thres2 | cnv_mat_denoise_bin>=thres1]=1
    # cnv_mat_denoise_bin[cnv_mat_denoise_bin!=1]=0


    x=as.data.frame(log2(cnv_mat)) %>% unlist()
    x1=x[x>0]
    x2=x[x<0]
    thres1=EBImage::otsu(array(x1,dim=c(length(x1),1)))
    thres2=-EBImage::otsu(array(abs(x2),dim=c(length(x2),1)))
    thres1
    thres2

    cnv_mat_bin=log2(cnv_mat)
    cnv_mat_bin[cnv_mat_bin>=thres2 & cnv_mat_bin<=0]=0
    cnv_mat_bin[cnv_mat_bin<=thres1 & cnv_mat_bin>=0]=0
    cnv_mat_bin[cnv_mat_bin!=1]=0
}

RunInferCNV <- function(seurat_obj=NULL,
                        raw_count_mat=NULL,cell_type_anno=NULL,
                        gene_ordering_table=NULL,ref_group_names='Ref_spot',
                        chr_exclude=c("chrX","chrY","chrM"),return_mat=TRUE,
						method='infercnv',
						num_threads=6){
						
	library(infercnv)
	
	# if (is.null(raw_count_mat)){
		# if (is.null(cell_type_anno))
	
	# }
	
    infercnv_obj=CreateInfercnvObject(raw_counts_matrix=raw_count_mat,
                                        annotations_file=cell_type_anno,
                                        gene_order_file=gene_ordering_table,
                                        ref_group_names=ref_group_names,
                                        chr_exclude=chr_exclude)
    if (method=='infercnv'){
		infercnv_obj=infercnv::run(infercnv_obj,
                                cutoff=0.1,
                                cluster_by_groups=FALSE,
                                sd_amplifier=1.5,
                                out_dir=save_dir,
                                denoise=TRUE,
                                save_rds=FALSE,
                                save_final_rds=TRUE,
                                HMM=FALSE)
	
	}
    if (method=='infercnv_HMM'){
	        infercnv_obj=infercnv::run(infercnv_obj,
                                   cutoff = 0.1,
                                   out_dir = save_dir,
                                   cluster_by_groups=TRUE,
                                   denoise=TRUE,
                                   HMM=TRUE,
                                   HMM_type="i6",
                                   analysis_mode='subcluster',
                                   BayesMaxPNormal=0,
                                   num_threads=num_threads,
                                   save_rds=FALSE,
                                   save_final_rds=TRUE,
                                   no_plot=FALSE)
	
	}

    if (return_mat){
        cnv_mat=infercnv_obj@expr.data %>% t() %>% data.frame(check.names=FALSE)
        return(cnv_mat)
    }else{
        return(infercnv_obj)
    }
}

# Run CNVestimate_ST
gene_ordering_table <- read.table('./gene_ordering_file.txt',
                                  header=FALSE,sep=" ",check.names=FALSE) %>%
   .[! duplicated(.[,1],fromLast=TRUE,),] %>%
   data.frame(row.names=1,check.names=FALSE)

CNVestimate_ST(seurat_obj,deconv_result,MalType=malignant_cell_type,
                out_dir=d,gene_ordering_table=gene_ordering_table,
                method='infercnv',ref_cell_ratio=0.01,StrImm_ratio=0.5,
                exponent=0.5,
                cluster_method='exp',scale.max=10,npcs=50,k.param=20,dims=1:10,resolution=0.3,
                cv_thres=0.4,
                marker=c(marker_ls[[tumor_type]],'COL1A1','PTPRC'))

# inspect the result using benign/malignant markers
d='CNV_result/'

df=lapply(dataset,function(dataset_1){
    library(dplyr)
    tumor_type=gsub('_[0-9]+$','',dataset_1)
    d=paste0(d,dataset_1,'/infercnv/')
    if (file.exists(paste0(d,'/spot_result.csv'))){
        f='spot_result.csv'
    }else{
        f=grep('spot_result*_res0.3',list.files(d),value=TRUE)
    }
    spot_df=read.csv(paste0(d,f),row.names=1,header=TRUE,check.names=FALSE)
    df=table(spot_df[,'CNVlabel']) %>% as.data.frame.array() %>% data.frame()
    
    colnames(df)=paste0('num_',dataset_1)
    df[,'attr']=rownames(df)
    df=df[,c(2,1)]
    
    return(df)
    
}) %>% Reduce(dplyr::full_join,.) %>% data.frame(row.names=1) %>% t()

rownames(df)=gsub('num_','',rownames(df))
df=apply(df,1,function(x){x/sum(x,na.rm=TRUE)}) %>% t() %>% data.frame()

df_marker=sapply(gsub('num_','',rownames(df)),function(dataset_1){
    tumor_type=gsub('_[0-9]+$','',dataset_1)
    d=paste0(d,dataset_1,'/infercnv/')
    if (file.exists(paste0(d,'/spot_result.csv'))){
        f='spot_result.csv'
    }else{
        f=grep('spot_result*_res0.3',list.files(d),value=TRUE)
    }
    spot_df=read.csv(paste0(d,f),row.names=1,header=TRUE,check.names=FALSE)
    seurat_obj=readRDS(file.path('../1_data_preprocessing',file.path(dataset_1,paste0(dataset_1,'.rds'))))
    Idents(seurat_obj)=spot_df[,'CNVlabel']
    if ('Benign' %in% unique(spot_df[,'CNVlabel'])){
        p=marker_check(seurat_obj,
                   benign_marker=marker_ls[[tumor_type]]['Benign'],
                   malignant_marker=marker_ls[[tumor_type]]['Malignant'])
    }else{
        p=0
    }
    return(c(sum(p<=0.05)/length(p)))
},USE.NAMES=TRUE)

df[,'Marker_check']=df_marker

df[,'CNVscore_check']=sapply(dataset,function(dataset_1){
    d=paste0(d,dataset_1,'/infercnv/')
    if (file.exists(paste0(d,'/spot_result.csv'))){
        f='spot_result.csv'
    }else{
        f=grep('spot_result*_res0.3',list.files(d),value=TRUE)
    }
    spot_df=read.csv(paste0(d,f),row.names=1,header=TRUE,check.names=FALSE)
    if ('Benign' %in% unique(spot_df[,'CNVlabel'])){
        p=t.test(filter(spot_df,CNVlabel=='Malignant') %>% select(CNVscore),
                 filter(spot_df,CNVlabel=='Benign') %>% select(CNVscore),)$ p.value
    }else{
        p=0
    }
    return(p<0.05)
})
