library(ggplot2)
library(ggalluvial)
library(ComplexHeatmap)
library(Seurat)

source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

graph_constr_10X <- function(seurat_obj,n=1,mat=NULL,to_binary=FALSE,to_igraph=TRUE,ncores=20){
    
    library(dplyr)
    library(igraph)
    
    if (class(seurat_obj@images$image)=='VisiumV1'){
        poi=seurat_obj@images$image@coordinates[,c('row','col')]
        poi[]=lapply(poi,as.numeric)
        # fix true coordinate of spots
        poi['row']=poi['row']*sqrt(3)
        d=dist(poi) %>% round(digits=5) %>% as.matrix()
        
    }

    if (class(seurat_obj@images$image)=='SlideSeq'){
        poi=seurat_obj@images$image@coordinates[,c('x','y')]
        poi[]=lapply(poi,as.numeric)
        d=dist(poi) %>% as.matrix()
        
    }
    
    # build graph based on distance (defined by how many rounds of spots centered by the target sopt)
    m=min(d[d!=0])*n
    bin_dist= d<=m & d!=0
    bin_dist[bin_dist==FALSE]=NA
    edge=reshape2::melt(bin_dist,na.rm=TRUE)
    colnames(edge)=c('Edge_1','Edge_2','Binary_weight')

    if (!to_binary){
        if (!is.null(mat)){
            edge[,'Binary_weight']=NULL
            # if (average){
            #     mat_adjut=sapply(rownames(mat),function(z){
            #         m=unlist(filter(edge,Edge_1==z) %>% select(Edge_2)) %>% as.character()
            #         apply(mat[c(z,m),],2,mean)
            #     }) %>% t()
            #     mat=mat_adjut
            # }

            # edge[,'weight']=apply(edge,1,function(y){
            #     1/dist(rbind(mat[y[1],],mat[y[2],]),method='euclidean')
            # })
            
            edge[,'weight']=parallel::mclapply(1:nrow(edge),function(z){
                y=edge[z,]
                c(1/dist(rbind(mat[y[,1],],mat[y[,2],]),method='euclidean'))
            },mc.cores=ncores) %>% unlist()

            # scale within whole network
            m=mean(edge[,'weight'])
            s=sd(edge[,'weight'])
            edge[,'weight']=(edge[,'weight']-m)/s
            edge[which(edge[,'weight']<=0),'weight']=0
            # edge[,'weight']=edge[,'weight']+abs(min(edge[,'weight']))

            # scele within node
            # edge[,'weight']=apply(edge,1,function(z){
            #     w1=filter(edge,Edge_1 %in% z[1])[,3]
            #     w2=filter(edge,Edge_1 %in% z[2])[,3]
            #     m1=mean(w1)
            #     s1=sd(w1)
            #     m2=mean(w2)
            #     s2=sd(w2)
            #     max((as.numeric(z[3])-m1)/s1,(as.numeric(z[3])-m2)/s2)
            # })
            # edge[which(edge[,'weight']<=0),'weight']=0

            # scele within node
            # edge[,'weight']=parallel::mclapply(edge,1,function(y){
            #     z=edge[y,]
            #     w1=filter(edge,Edge_1 %in% z[,1])[,3]
            #     w2=filter(edge,Edge_1 %in% z[,2])[,3]
            #     m1=mean(w1)
            #     s1=sd(w1)
            #     m2=mean(w2)
            #     s2=sd(w2)
            #     (z[,3]-m1)/s1
            #     max((as.numeric(z[3])-m1)/s1,(as.numeric(z[3])-m2)/s2)
            # },mc.cores=ncores)
            # edge[which(edge[,'weight']<=0),'weight']=0


            # edge[,'weight']=parallel::mclapply(1:nrow(edge),function(z){
            #     x=edge[z,]
            #     w1=filter(x,Edge_1 %in% x[,1])[,3]
            #     w2=filter(x,Edge_1 %in% x[,2])[,3]
            #     m1=mean(w1)
            #     s1=sd(w1)
            #     m2=mean(w2)
            #     s2=sd(w2)
            #     # max((as.numeric(x[,3])-m1)/s1,(as.numeric(x[,3])-m2)/s2)
            #     (as.numeric(x[,3])-m1)/s1
            # },mc.cores=50)
            # edge[which(edge[,'weight']<=0),'weight']=0

            # Banksy style
            # edge[,'weight']=apply(edge,1,function(z){
            #     w1=filter(edge,Edge_1 %in% z[1])[,3]
            #     w2=filter(edge,Edge_1 %in% z[1])[,3]
            #     as.numeric(z[3])/sum(min(w1,w2))
            # })

            # Banksy style
            # edge[,'weight']=parallel::mclapply(1:nrow(edge),function(z){
            #     x=edge[z,,drop=TRUE]
            #     w1=filter(edge,Edge_1 %in% x[1])[,3]
            #     w2=filter(edge,Edge_1 %in% x[1])[,3]
            #     as.numeric(x[3])/sum(min(w1,w2))
            # },mc.cores=20)
            }
        }
    
    if (to_igraph){
        return(graph.data.frame(edge,directed=FALSE))
    }else{
        return(edge)
    }
}

niche_type=lapply(dataset_visium,function(i){
    print(i)
    seurat_obj=readRDS(paste0('../../1_data_preprocessing/',i,'/',i,'.rds'))

    if (ncol(seurat_obj)<=100){
        return(NULL)
    }
    ind=intersect(colnames(seurat_obj),rownames(LCP_cluster_result))
    if (length(ind)==0){
        return(NULL)
    }
    seurat_obj=seurat_obj[,ind]
    spatial_graph=graph_constr_10X(seurat_obj,n=2,to_igraph=FALSE)
    
    LCP_cluster_result_tmp=LCP_cluster_result[ind,]
    m=parallel::mclapply(1:nrow(spatial_graph),function(i){
        y=spatial_graph[i,]
        return(data.frame(LCP_1=LCP_cluster_result_tmp[as.character(y[1,1]),'cluster_annotated'],
                          LCP_2=LCP_cluster_result_tmp[as.character(y[1,2]),'cluster_annotated']))
    },mc.cores=20) %>% dplyr::bind_rows()

    spatial_graph=cbind(spatial_graph,m)

    a=unique(spatial_graph[,c('Edge_1','LCP_1')])
    colnames(a)[2]='LCP_2'
    
    niche_tab=dplyr::bind_rows(spatial_graph[,c('Edge_1','LCP_2')],a) %>% 
        reshape2::dcast(Edge_1~LCP_2) %>% data.frame(row.names=1,check.names=FALSE) %>% apply(1,function(x){x/sum(x)}) %>% 
        t() %>% data.frame(check.names=FALSE)
    
}) %>% dplyr::bind_rows()

cluster_num=6:50
niche_cluster_df=sapply(cluster_num,function(x){
    niche_cluster_result <- kmeans(niche_spot,x,iter.max=100000,nstart=40)
    return(niche_cluster_result$cluster)
})
colnames(niche_cluster_df)=paste0('N_',cluster_num)


# LCP composition of niches
niche_cluster_comp <- group_by(data.frame(niche_spot,check.names=FALSE),as.factor(niche_cluster_result[,'Niche_combined'])) %>% 
    summarize_if(is.numeric,mean) %>% data.frame(row.names=1,check.names=FALSE)
Heatmap(niche_cluster_comp,name='Percantage',clustering_method_rows='ward.D2',cluster_columns = FALSE,cluster_rows=FALSE,
        # show_row_names=FALSE,show_column_names=FALSE,
        col=circlize::colorRamp2(c(seq(0,0.1,length.out=3),c(seq(0.2,0.5,length.out=4))),
                                 c('#F4FAED','#D6EFD0','#B2E1B9','#77CAC5','#42A6CB','#1373B2','#084384'))
)
# cell type composition of niches
x=read.csv('./Version 4/niche_cluster_result.csv',row.names=1,check.names=FALSE)
y=read.csv('../3_LCP_definition/Version 3/spot_type.csv',row.names=1,check.names=FALSE)
y=y[rownames(x),]
z=y %>% group_by(x[,'Niche_combined']) %>% summarise_if(is.numeric,mean) %>% data.frame(row.names=1)
z1=z[paste0('Niche_',1:13),c('Parachymal.cell','Endothelial.cell','Fibroblast','Macrophage',
                             'cDC','CD4.T','CD8.T','B','Plasma.cell','NK','pDC','Mast.cell')]
Heatmap(scale(z1),cluster_columns=FALSE,cluster_rows=FALSE)

# correspondance
LCP_cluster_result[,'Niche']=niche_cluster_result[rownames(LCP_cluster_result),'Niche_combined']

df=reshape2::dcast(LCP_cluster_result,Niche+cluster_annotated~1)
colnames(df)[3]='weight'

df[,'Niche_1']=df[,'Niche']
df_1=to_lodes_form(df,key='Demographic',weight='Niche',axes=1:2)
df_1=na.omit(df_1)
df_1$stratum=factor(df_1$stratum,levels=c(rownames(niche_cluster_comp),colnames(niche_cluster_comp)))

options(repr.plot.height=20)
p=ggplot(df_1,aes(x=Demographic,y=weight,stratum=stratum,
                alluvium=alluvium,fill=stratum,label=stratum))+
    geom_flow(aes(fill=Niche_1))+
    geom_stratum()+
    geom_text(stat='stratum',aes(label=after_stat(stratum)))+
    guides(fill='none')+
    theme_void()+
    # ggsci::scale_fill_igv()
    scale_fill_manual(values=c(MLCP_color,SLCP_color,niche_color))

# intra-tumor heterogeneity - entropy
E=apply(niche_sample_comp_,1,entropy::entropy)
tumor_type=gsub('_[0-9]+$','',names(E))
df=data.frame(E=E,tumor_type=tumor_type)
ggplot(data=df,aes( x=reorder(tumor_type,E,median),y=E,fill=tumor_type ) )+
    geom_boxplot()+
    scale_fill_manual(values=tumor_color)+
    guides(fill='none')+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1),axis.title.x=element_blank())

# inter-sample heterogeneity - sample clustering
niche_sample_comp=table(niche_cluster_result$sample,niche_cluster_result$Niche_combined)%>% as.data.frame.array() %>% 
    t() %>% apply(2,function(x){x/sum(x)}) %>% t() %>% data.frame()
# determine cluster number
SimilarityClustering(niche_sample_comp,mode='manual',select.cutoff=TRUE,
                    hc.method='ward.D2',similarity.method='pearson',min.nc=2,max.nc=20,show_index_result=NA,
                    # right_annotation=rowAnnotation(' '=gsub('_[0-9]+','',rownames(niche_sample_comp_clustering)),col=list(' '=tumor_color))
                 )
options(repr.plot.width=7,repr.plot.height=7)
c=SimilarityClustering(niche_sample_comp_clustering,mode='manual',select.cutoff=FALSE,
                    hc.method='ward.D2',similarity.method='pearson',
                    cluster_num=11,
                    right_annotation=rowAnnotation(' '=gsub('_[0-9]+','',rownames(niche_sample_comp_clustering)),col=list(' '=tumor_color))
                 )

# niche type ~ cancer type
df=data.frame(c)
tumor_type_order=reshape2::dcast(df,tumor_type~MetaType,value.var='sample_num',fun.aggregate=sum) %>% data.frame(row.names=1,check.names=FALSE) %>% 
    apply(1,function(x){x/sum(x)}) %>% t()%>% data.frame(check.names=FALSE) %>% arrange(ImmuneType) %>% rownames()

df_2=df %>% group_by(tumor_type,SpaType) %>% summarise(Num=sum(Num=sample_num)) %>% data.frame()
df_2=df_2 %>% group_by(tumor_type) %>% mutate(Ratio=Num/sum(Num)) %>% data.frame()