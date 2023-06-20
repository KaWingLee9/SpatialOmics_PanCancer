graph_constr_10X <- function(seurat_obj,n=1,mat=NULL,to_binary=FALSE,to_igraph=TRUE,ncores=20){
    
    library(dplyr)
    library(igraph)
    
    poi=seurat_obj@images$image@coordinates[,c('row','col')]
    # fix true coordinate of spots
    poi['row']=poi['row']*sqrt(3)
    d=dist(poi) %>% round(digits=5) %>% as.matrix()
    
    # build graph based on distance (defined by how many rounds of spots centered by the target sopt)
    m=2*n
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

neighbour_average <- function(seurat_obj,mat,n=1,lambda=0.5,similarity_weighted=TRUE){
    
    edge=graph_constr_10X(seurat_obj,n=n,mat=NULL,to_binary=TRUE,to_igraph=FALSE)
    
    # Banksy style
    x=lapply(rownames(mat),function(y){
        n1=y
        n2=filter(edge,Edge_1==n1) %>% .[,'Edge_2']
        mat1=mat[n1,]
        mat2=mat[n2,]
        if (similarity_weighted){
            w=apply(mat2,1,function(x){1/dist(rbind(x,mat1))})
            w=w/sum(w)
            mat2=w*mat2
        }
        mat1_adjust=mat1*lambda+apply(mat2*(1-lambda),2,sum)
        return(mat1_adjust)
    }) %>% dplyr::bind_rows()
    return(x)
}
