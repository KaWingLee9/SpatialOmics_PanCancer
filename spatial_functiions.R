# construct spatial network from seurat_obj
graph_constr_10X <- function(seurat_obj,n=1,to_binary=FALSE,to_igraph=TRUE,ncores=20){

    ### Usage for graph_constr_10X:
    # seurat_obj: SeuratV4 object
    # n: the round number of the spots considered
    # to_binary: whether the weight of the edge is binary
    # to_igraph: return igraph object; otherwise return a data frame
    # ncores: the CPU used 
    
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
        # if (!is.null(mat)){
            mat=poi
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
            # }
        }
    
    if (to_igraph){
        return(graph.data.frame(edge,directed=FALSE))
    }else{
        return(edge)
    }
}

# sptially weighted expression
spatialweighted_exp=function(seurat_obj,n=1,
                             assay='RNA',new_assay='weighted',
                             ncores=20,...){
    coord_df=graph_constr_10X(seurat_obj,n=n,to_binary=TRUE,to_igraph=FALSE,ncores=ncores)
    coord_df[,'Edge']=1
    coord_df=coord_df %>% reshape2::dcast(Edge_1~Edge_2,value.var='Edge') %>% data.frame(row.names=1,check.names=FALSE)
    coord_df[is.na(coord_df)]=0
    x=setdiff(x,rownames(coord_df))
    coord_df[,x]=0
    coord_df[x,]=0

    x=intersect(rownames(coord_df),colnames(seurat_obj))

    exp_mat=seurat_obj@assays$RNA@data %>% as.matrix()
    exp_mat=exp_mat[,x]
    coord_df=coord_df[x,x] %>% as.matrix()

    exp_lag=(coord_df %*% t(exp_mat))/apply(coord_df,1,sum)
    mat=t((t(exp_mat)+exp_lag)/2)
    mat=mat[rownames(seurat_obj),colnames(seurat_obj)]
    seurat_obj[['weighted']]=CreateAssayObject( mat )
    return(seurat_obj)
}

# Double Layer visualization: Background (Discrete), Foreground(Discrete); return a ggplot object
DL_DimPlot=function(seurat_obj,foreground_label='LCP',foreground_selected_label='',reverse=FALSE,
                    background_label='LCP_m1',background_size=2,foreground_size=0.75
                    ){
    metadata=seurat_obj@meta.data
    if (class(seurat_obj@images$image)=='VisiumV1'){
        metadata['imagerow']=seurat_obj@images$image@coordinates$imagecol
        metadata['imagecol']=-seurat_obj@images$image@coordinates$imagerow
    }
    if (class(seurat_obj@images$image)=='SlideSeq'){
        metadata['imagerow']=seurat_obj@images$image@coordinates$x
        metadata['imagecol']=-seurat_obj@images$image@coordinates$y
    }

    
    if (foreground_selected_label!=''){
        ind= ! metadata[,foreground_label] %in% foreground_selected_label
        if (reverse){
            ind=! ind
        }
        metadata[ind,'LCP']=NA
    }
    
    foreground_size=1-foreground_size
    
    p=ggplot(metadata,aes_string(x='imagerow',y='imagecol'))+
      geom_point(aes_string(color=background_label),size=background_size)+
      geom_point(data=na.omit(metadata),
                 aes_string(color=foreground_label),size=foreground_size*background_size)+
      theme_classic()+
      theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
    
    return(p)
}

# Double Layer visualization: Background (Discrete), Foreground(Continuous); return a ggplot object
DL_FeaturePlot=function(seurat_obj,foreground_label='KRT8',foreground_assay='RNA',foreground_slot='data',
                        scale=TRUE,min.cutoff=NULL,max.cutoff=NULL,
                        background_label='LCP_m1',
                        show_background_value=NA,highlight_background=NA,
                        background_size=2,foreground_size=0.75,
                        scale_size=FALSE,scale_alpha=FALSE,rm.min=TRUE
                        ){
    
    metadata=seurat_obj@meta.data
    if (class(seurat_obj@images$image)=='VisiumV1'){
        metadata['imagerow']=seurat_obj@images$image@coordinates$imagecol
        metadata['imagecol']=-seurat_obj@images$image@coordinates$imagerow
    }
    if (class(seurat_obj@images$image)=='SlideSeq'){
        metadata['imagerow']=seurat_obj@images$image@coordinates$x
        metadata['imagecol']=-seurat_obj@images$image@coordinates$y
    }
    
    if (foreground_label %in% rownames(seurat_obj@assays[[foreground_assay]])){
        metadata['value']=GetAssayData(seurat_obj,assay=foreground_assay,slot=foreground_slot)[foreground_label,]
    } else {
        metadata['value']=0
    }
    

    if (scale){
        metadata['value']=scale(metadata['value'])
    }
    
    if (! is.null(min.cutoff)){
        if (is.character(min.cutoff)){
            q=gsub('q','',min.cutoff)
            q=as.numeric(q)
            min.cutoff=quantile(metadata[,'value'],q*0.01)
            metadata[metadata['value']<=min.cutoff,'value']=min.cutoff
        }
        if (min.cutoff>=min(metadata['value'])){
            metadata[metadata['value']<=min.cutoff,'value']=min.cutoff
    }}

    if (! is.null(max.cutoff)){
        if (is.character(max.cutoff)){
            q=gsub('q','',max.cutoff)
            q=as.numeric(q)
            min.cutoff=quantile(metadata[,'value'],q*0.01)
            metadata[metadata['value']>=max.cutoff,'value']=max.cutoff
        }
        if (max.cutoff<=max(metadata['value'])){
            metadata[metadata['value']>=max.cutoff,'value']=max.cutoff
    }}
    
    if (rm.min){
        metadata[metadata['value']==min(metadata['value']),'value']=NA
    }
    
    if (!is.na(show_background_value)){
        metadata[ ! metadata[,background_label] %in% show_background_value,'value']=NA
    }

    if (!is.na(highlight_background)){
        metadata[ ! metadata[,background_label] %in% highlight_background , background_label]='Other'
    }

    if (scale_size){
        scale_size='value'
    } else if (! scale_size){
        scale_size=NULL
    }

    if (scale_alpha){
        scale_alpha='value'
    } else if (! scale_alpha){
        scale_alpha=NULL
    }

    # if (scale_size_alpha){
    #     p=ggplot(metadata,aes_string(x='imagerow',y='imagecol'))+
    #       geom_point(aes_string(color=background_label),size=background_size)+
    #       geom_point(data=na.omit(metadata),
    #                  aes_string(fill='value',size=scale_size,alpha=scale_alpha
    #                             ),stroke=0,color='#FFFFFF00',shape=21)+
    #       scale_size_area(max_size=foreground_size*background_size)+
    #       theme_classic()+
    #       # scale_color_manual(values=c('Benign'='#CFFFC6','MLCP'='#FFFCC6','SLCP'='#C6E2FF'))+
    #       theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
    #       scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9,'Reds'))+
    #       labs(fill=foreground_label)+
    #       guides(alpha='none',size='none')
    # }else{
    #     p=ggplot(metadata,aes_string(x='imagerow',y='imagecol'))+
    #       geom_point(aes_string(color=background_label),size=background_size)+
    #       geom_point(data=na.omit(metadata),
    #                  aes_string(fill='value'),color='#FFFFFF00',shape=21,stroke=0,size=foreground_size*background_size)+
    #       theme_classic()+
    #       # scale_color_manual(values=c('Benign'='#CFFFC6','MLCP'='#FFFCC6','SLCP'='#C6E2FF'))+
    #       theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
    #       scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9,'Reds'))+
    #       labs(fill=foreground_label)
    # }

    p=ggplot(metadata,aes_string(x='imagerow',y='imagecol'))+
        geom_point(aes_string(color=background_label),size=background_size)+
        geom_point(data=na.omit(metadata),
                    aes_string(fill='value',size=scale_size,alpha=scale_alpha
                            ),stroke=NA,color='#FFFFFF00',shape=21)+
        # scale_size_area(max_size=foreground_size*background_size)+
        theme_classic()+
        # scale_color_manual(values=c('Benign'='#CFFFC6','MLCP'='#FFFCC6','SLCP'='#C6E2FF'))+
        theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
        scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9,'Reds'))+
        labs(fill=foreground_label)+
        guides(alpha='none',size='none')

    if (! is.null(scale_size)){
        p=p+scale_size_area(max_size=foreground_size*background_size)
    }

    return(p)
}