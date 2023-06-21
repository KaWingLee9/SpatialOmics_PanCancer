Binarization <- function(x,method='',binarization=FALSE,Plot=FALSE){
    if (method=='otsu'){
        thres=EBImage::otsu(array(x,dim=c(length(x),1)))
    }
    if (method=='2sd'){
        thres=mean(x)+2*sd(x)
    }
    if (method=='GMM'){
        
        thres
    }
    if (Plot){
        plot(density(x),main=paste0('thres=',thres))
        abline(v=thres,col='red')
    }

    x[x<=thres]=0
    if (binarization==TRUE){
        x[x>=thres]=1
    }
    return(x)
}
