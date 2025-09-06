library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(ggridges)
library(ggforce)
library(igraph)
library(tidygraph)
library(ggraph)
library(Seurat)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

# compartment connectivity analysis
CalCS=function(seurat_obj,n=1,ncores=20){
    
    spatial_graph=graph_constr_10X(seurat_obj,n=n,to_igraph=FALSE,ncores=ncores)
    spatial_graph=apply(spatial_graph,1,function(x) sort(c(x[1],x[2]))) %>% t() %>% data.frame() %>% unique()
    colnames(spatial_graph)=c('Node_1','Node_2')
    spatial_graph[,'LCP_1']=Idents(seurat_obj)[spatial_graph[,'Node_1']]
    spatial_graph[,'LCP_2']=Idents(seurat_obj)[spatial_graph[,'Node_2']]

    spatial_graph[,c('LCP_1','LCP_2')]=apply(spatial_graph,1,function(x) {sort(c(x['LCP_1'],x['LCP_2']))}) %>% t()
    spatial_graph[,'Edges']=1
    
    LCP_num=c(table(Idents(seurat_obj)))
#     LCP_num=LCP_num/sum(LCP_num)
    
    LCP_edges=sapply(Idents(seurat_obj) %>% unique() %>% as.character(),function(LCP_x){
        filter(spatial_graph,(LCP_1==LCP_x) | (LCP_2==LCP_x)) %>% nrow
    },USE.NAMES=TRUE)                 
    
    connectivity_tbl=reshape2::dcast(spatial_graph,LCP_1+LCP_2~Edges,value.var='Edges',fun.aggregate=sum)
    colnames(connectivity_tbl)[3]='Edges'
    
    connectivity_tbl[,'CS_score']=apply(connectivity_tbl,1,function(x){
        # (as.numeric(x[3])) /sqrt( LCP_num[x[1]] * LCP_num[x[2]]) # CS_score_old.csv
        (as.numeric(x[3])) /sqrt( LCP_edges[x[1]] * LCP_edges[x[2]]) # CS_score_by_edge.csv
        # (as.numeric(x[3])) / nrow(connectivity_tbl) # CS_score_by_all_edge.csv
    })
                        
    connectivity_tbl=data.frame(connectivity_tbl,check.names=FALSE)              
                        
    return(connectivity_tbl)
                        
}

connectivity_tbl=lapply(dataset_used,function(dataset){
    
    seurat_obj=readRDS(paste0('../1_data_preprocessing/',dataset,'/',dataset,'.rds'))
    ind=colnames(seurat_obj)
    Idents(seurat_obj)=CNV_cluster_result[ind,'CNVlabel']

    connectivity_tbl=CalCS(seurat_obj,n=1)
    connectivity_tbl[,'dataset']=dataset
    
    x=data.frame(true_label=Idents(seurat_obj))
    x=arrange(x,true_label)

    seurat_obj=seurat_obj[,rownames(x)]
    x$group=gsub('_[ /A-Za-z0-9]+$','',x$true_label)

    shuffle_labels=lapply(1:20,function(i){
        x=x %>% mutate('sample_label'=sample(true_label))
        return(x$sample_label)
    })

    shuffle_connectivity_tbl=lapply(shuffle_labels,function(y){

        Idents(seurat_obj)=y
        connectivity_tbl=CalCS(seurat_obj,n=1)

    }) %>% dplyr::bind_rows() %>% group_by(LCP_1,LCP_2) %>% summarise(m=mean(Edges),s=sd(Edges))
    
    connectivity_tbl[,'zscore']=apply(connectivity_tbl,1,function(x){
        (as.numeric(x['Edges'])-as.numeric(shuffle_connectivity_tbl[ (shuffle_connectivity_tbl[,'LCP_1']==x['LCP_1']) & (shuffle_connectivity_tbl[,'LCP_2']==x['LCP_1']) ,'m']))/as.numeric(shuffle_connectivity_tbl[ (shuffle_connectivity_tbl[,'LCP_1']==x['LCP_1']) & (shuffle_connectivity_tbl[,'LCP_2']==x['LCP_1']) ,'s'])
    })

    connectivity_tbl[is.na(connectivity_tbl[,'zscore']),'zscore']=0
    connectivity_tbl[,'dataset']=dataset
    
    return(connectivity_tbl)
    
}) %>% dplyr::bind_rows()

# LCP connectivity analysis
CalCS=function(seurat_obj,n=1,ncores=20){
    
    spatial_graph=graph_constr_10X(seurat_obj,n=n,to_igraph=FALSE,ncores=ncores)
    spatial_graph=apply(spatial_graph,1,function(x) sort(c(x[1],x[2]))) %>% t() %>% data.frame() %>% unique()
    colnames(spatial_graph)=c('Node_1','Node_2')
    spatial_graph[,'LCP_1']=Idents(seurat_obj)[spatial_graph[,'Node_1']]
    spatial_graph[,'LCP_2']=Idents(seurat_obj)[spatial_graph[,'Node_2']]

    spatial_graph[,c('LCP_1','LCP_2')]=apply(spatial_graph,1,function(x) {sort(c(x['LCP_1'],x['LCP_2']))}) %>% t()
    spatial_graph[,'Edges']=1
    
    LCP_num=c(table(Idents(seurat_obj)))
#     LCP_num=LCP_num/sum(LCP_num)
                        

                        
    LCP_edges=sapply(Idents(seurat_obj) %>% unique() %>% as.character(),function(LCP_x){
        filter(spatial_graph,(LCP_1==LCP_x) | (LCP_2==LCP_x)) %>% nrow
    },USE.NAMES=TRUE)             
                        
    connectivity_tbl=reshape2::dcast(spatial_graph,LCP_1+LCP_2~Edges,value.var='Edges',fun.aggregate=sum)
    colnames(connectivity_tbl)[3]='Edges'
                        
    connectivity_tbl$LCP_1_group=gsub('_[/ A-Za-z0-9]+$','',connectivity_tbl$LCP_1)
    connectivity_tbl$LCP_2_group=gsub('_[/ A-Za-z0-9]+$','',connectivity_tbl$LCP_2)  
                        
    connectivity_tbl=data.frame(connectivity_tbl,check.names=FALSE)
    c1=combn(sort(unique(c(connectivity_tbl$LCP_1,connectivity_tbl$LCP_2))),2) %>% t() %>% data.frame()
    colnames(c1)=c('LCP_1','LCP_2')
    c1$Edges=0
    # c1$CS_score=0
    c1$LCP_1_group=gsub('_[/ A-Za-z0-9]+$','',c1$LCP_1)
    c1$LCP_2_group=gsub('_[/ A-Za-z0-9]+$','',c1$LCP_2)                         

    
    c2=data.frame(LCP_1=unique(c(connectivity_tbl$LCP_1,connectivity_tbl$LCP_2)),LCP_2=unique(c(connectivity_tbl$LCP_1,connectivity_tbl$LCP_2)),
                  Edges=0# ,CS_score=0
                 )
    c2$LCP_1_group=gsub('_[/ A-Za-z0-9]+$','',c2$LCP_1)
    c2$LCP_2_group=gsub('_[/ A-Za-z0-9]+$','',c2$LCP_2)    
    connectivity_tbl=dplyr::bind_rows(list(connectivity_tbl,c1,c2))      
                        
    connectivity_tbl=connectivity_tbl %>% group_by(LCP_1,LCP_2,LCP_1_group,LCP_2_group) %>% summarise(Edges=sum(Edges))
                        
    connectivity_tbl$group=paste0(connectivity_tbl$LCP_1_group,'-',connectivity_tbl$LCP_2_group)                        
    total_num=reshape2::dcast(connectivity_tbl,group~1,value.var='Edges',fun.aggregate=sum)
    rownames(total_num)=total_num$group
                        
    connectivity_tbl[,'CS_score']=apply(connectivity_tbl,1,function(x){
        # v1
        (as.numeric(x[5]))/(total_num[x[6],2])
    })
                        
#     connectivity_tbl[,'CS_score']=apply(connectivity_tbl,1,function(x){
#         # v3
#         (as.numeric(x[5])) /sqrt( LCP_num[x[1]] * LCP_num[x[2]]) # can not be used across technologies
#         # (as.numeric(x[5])) /sqrt( LCP_edges[x[1]] * LCP_edges[x[2]]) # CS_score_by_edge.csv
#         # (as.numeric(x[5])) / nrow(connectivity_tbl) # CS_score_by_all_edge.csv
#     })
                        
    return(connectivity_tbl)
                        
}

connectivity_tbl=lapply(dataset_used,function(y){
    
    seurat_obj=readRDS(paste0('../1_data_preprocessing/',y,'/',y,'.rds'))
    ind=intersect(colnames(seurat_obj),rownames(LCP_cluster_result))

    seurat_obj=seurat_obj[,ind]
    Idents(seurat_obj)=LCP_cluster_result[ind,'cluster_annotated']
    connectivity_tbl=CalCS(seurat_obj,n=2)

    x=data.frame(true_label=Idents(seurat_obj))
    x=arrange(x,true_label)

    seurat_obj=seurat_obj[,rownames(x)]
    x$group=gsub('_[ /A-Za-z0-9]+$','',x$true_label)

    shuffle_labels=lapply(1:20,function(i){
        x=x %>% group_by(group) %>% mutate('sample_label'=sample(true_label))
        return(x$sample_label)
    })
    
    shuffle_connectivity_tbl=lapply(shuffle_labels,function(y){

        Idents(seurat_obj)=y
        connectivity_tbl=CalCS(seurat_obj,n=2)

    }) %>% dplyr::bind_rows() %>% group_by(LCP_1,LCP_2) %>% summarise(m=mean(Edges),s=sd(Edges))

    connectivity_tbl[,'zscore']=apply(connectivity_tbl,1,function(x){
        (as.numeric(x['Edges'])-as.numeric(shuffle_connectivity_tbl[ (shuffle_connectivity_tbl[,'LCP_1']==x['LCP_1']) & (shuffle_connectivity_tbl[,'LCP_2']==x['LCP_1']) ,'m']))/as.numeric(shuffle_connectivity_tbl[ (shuffle_connectivity_tbl[,'LCP_1']==x['LCP_1']) & (shuffle_connectivity_tbl[,'LCP_2']==x['LCP_1']) ,'s'])
    })

    connectivity_tbl[is.na(connectivity_tbl[,'zscore']),'zscore']=0
    connectivity_tbl[,'dataset']=y
    
    return(connectivity_tbl)
    
}) %>% dplyr::bind_rows()

# LCP self-connectivity

# LCP community
g=group_by(CS_score_filter,LCP_1,LCP_2) %>% summarise(mean_CS=mean(CS_score),pval=wilcox.test(zscore,mu=0,alternative='greater')$p.value
                                              ) %>%
    filter(pval<=0.05) %>% 
    .[,c('LCP_1','LCP_2','mean_CS')] # %>% filter(LCP_1!=LCP_2)
g=na.omit(g)

g=graph_from_data_frame(g,directed=FALSE)

vertex_attr(g,'node')=gsub('^[A-Z]+_','',vertex_attr(g,'name'))

g=delete_edges(g,E(g)[edge_attr(g,'mean_CS')<0.001])

x=cluster_louvain(g,weights=edge_attr(g,'mean_CS'),resolution=2.3)
vertex_attr(g,'Community')=as.factor(c(membership(x))[vertex_attr(g,'name')])
g=delete_edges(g,E(g)[edge_attr(g,'mean_CS')<0.005])

p=ggraph(g,layout='stress',weights=edge_attr(g,'mean_CS')
        )+
    geom_edge_fan2(aes(color=mean_CS,alpha=mean_CS))+
#     geom_edge_loop(aes(color=mean_CS,alpha=mean_CS,
#                            direction=0,strength=0.6,span=100))+
    geom_node_point(aes(color=Community))+ # ,size=3
    geom_node_text(aes(label=paste0(name,' (C_',Community,')'),color=Community),repel=TRUE)+
    theme_void()+
    scale_edge_color_gradientn(colours=RColorBrewer::brewer.pal(9,'Greys'),
                               values=scales::rescale(quantile(unlist(edge_attr(g,'mean_CS')),probs=(0:8)/8)))+
#     geom_mark_hull(aes(x,y,group=vertex_attr(g,'Community'),fill=vertex_attr(g,'Community')),
#                    alpha=0.25,show.legend=FALSE,expand=unit(3,'mm'))+
#     facet_nodes(.~tumor_type,scales='free')+
    theme(panel.border=element_rect(colour='black',fill=NA),strip.text=element_text(size=rel(1.5)))+
    guides(color='none')

