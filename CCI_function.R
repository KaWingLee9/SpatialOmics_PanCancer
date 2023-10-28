# LR score for all spots
CalSpotLRScore <- function(exp_mat,lr_ls,min_spot=10){
    
    library(Seurat)
    
    # LR score for spots
    lr_score_spot=sapply(1:length(lr_ls),function(i){

        lr_pair=lr_ls[[i]]
        
        # filter out LR with no genes detected
        if (length(setdiff(lr_pair,rownames(exp_mat)))!=0){
            return(NULL)
        }
        
        # geometric mean of all genes in LR pair
        df=apply(exp_mat[lr_pair,],2,function(x){
            prod(x)^(1/length(x))
        }) %>% data.frame()

        colnames(df)=names(lr_ls)[i]
        return(df)
    }) %>% dplyr::bind_cols() %>% data.frame(row.names=colnames(exp_mat),check.names=FALSE)

    # filter efficient L-R: number of spots with positive score not lower than `min_spot`
    LR_used=which(apply(lr_score_spot,2,function(x){sum(x>0)})>=min_spot) %>% names()
    lr_score_spot=lr_score_spot[,LR_used]
    
    return(as.sparse(t(lr_score_spot)))
    
}

# LR coexpression between cell types of spots
CalSampleCCI <- function(lr_score_spot,cell_type_score,cor_method='pearson',
                         p.adjust.method='bonferroni'){
    
    pearson_result=Hmisc::rcorr(as.matrix(cell_type_score),lr_score_spot,type=cor_method)
    r=pearson_result$r[colnames(cell_type_score),colnames(lr_score_spot)] %>% reshape2::melt(varnames=c('cell_type','lr_pair'),value.name='pearson_r')
    p=pearson_result$P[colnames(cell_type_score),colnames(lr_score_spot)] %>% reshape2::melt(varnames=c('cell_type','lr_pair'),value.name='pearson_p')
    cor_result=dplyr::left_join(r,p,by=c('cell_type','lr_pair'))
    cor_result[which(cor_result[,'pearson_r']<=0),'pearson_p']=1
    cor_result[,'pearson_bonferroni']=p.adjust(cor_result[,'pearson_p'],method=p.adjust.method)
    return(cor_result)
    
}

# LCP and cell type related LR pairs
CalLCPCCI <- function(lr_score_spot,cell_type_score,MLCP_result,cor_result,
                      ref_LCP,obs_LCP,LCP_cutoff=20){
    test_result=sapply(obs_LCP,function(obs_LCP){
        obs_spot=filter(MLCP_result,cluster_result==obs_LCP) %>% rownames()
        ref_spot=filter(MLCP_result,cluster_result==ref_LCP) %>% rownames()
        if (length(obs_spot)<=LCP_cutoff | length(ref_spot)<=LCP_cutoff){
            return(NULL)
        }
        cell_type_score_tmp=cell_type_score[c(ref_spot,obs_spot),]
        x=sapply(colnames(cell_type_score_tmp),function(i){
            y=cell_type_score_tmp[,i]
            df_lm=data.frame(character=y)
            df_lm[,'LCP']=rep(c(ref_LCP,obs_LCP),times=c(length(ref_spot),length(obs_spot))) %>% 
            factor(levels=c(ref_LCP,obs_LCP))
            lm_model=lm(character~LCP,data=df_lm) %>% summary()
            LCP_test_result=lm_model$coefficients %>% .[nrow(.),]
            if (LCP_test_result['Estimate']<=0 | LCP_test_result['Pr(>|t|)']>=0.05){
                return(NULL)
            }
            lr_score_spot_tmp=lr_score_spot[c(ref_spot,obs_spot),] %>% data.frame(check.names=FALSE)
            lr_score_spot_tmp=lr_score_spot_tmp[,filter(cor_result,cell_type==i) %>% .[,'lr_pair'],drop=FALSE]
            lr=sapply(lr_score_spot_tmp,function(z){
                df_lm=data.frame(character=z)
                df_lm[,'LCP']=rep(c(ref_LCP,obs_LCP),times=c(length(ref_spot),length(obs_spot))) %>% 
                    factor(levels=c(ref_LCP,obs_LCP))
                lm_model=lm(character~LCP,data=df_lm) %>% summary()
                LCP_test_result=lm_model$coefficients %>% .[nrow(.),]
                if (LCP_test_result['Estimate']<=0 | LCP_test_result['Pr(>|t|)']>=0.05){
                    return(NULL)
                }
                return(1)
            })
            lr=lr[!sapply(lr,is.null)]
            lr=names(lr)
            return(lr)
        })
        x=x[!sapply(x,is.null)]
        return(x)
    }) %>% .[!sapply(.,is.null)]
    
    if (length(test_result)==0){
        return(NULL)
    }
    
    test_result=test_result %>% .[sapply(.,function(x){length(x)!=0})]
    
    test_result=reshape2::melt(test_result)
    test_result=data.frame(test_result)
    
    if (nrow(test_result)==0){
        return(NULL)
    }
    
    colnames(test_result)=c('sig_lr_pair','cell_type','LCP')
    return(test_result)
}

# LCP-related LR pairs
CalLCPCCI_V2 <- function(lr_score_spot,LCP_result,
                      ref_LCP,obs_LCP,LCP_cutoff=20){
    test_result=lapply(obs_LCP,function(obs_LCP){
        obs_spot=filter(LCP_result,cluster_result==obs_LCP) %>% rownames()
        ref_spot=filter(LCP_result,cluster_result==ref_LCP) %>% rownames()
        if (length(obs_spot)<=LCP_cutoff | length(ref_spot)<=LCP_cutoff){
            return(NULL)
        }

        lr_score_spot_tmp=lr_score_spot[c(ref_spot,obs_spot),] %>% data.frame(check.names=FALSE)

        lr_test=sapply(lr_score_spot_tmp,function(z){
            df_lm=data.frame(character=z)
            df_lm[,'LCP']=rep(c(ref_LCP,obs_LCP),times=c(length(ref_spot),length(obs_spot))) %>% 
                factor(levels=c(ref_LCP,obs_LCP))
            if (length(unique(z))==1){
                test_result=c('Estimate'=0.0000000001,'Std. Error'=0,'Std. Error'=0,'Pr(>|t|)'=1)
            }else{
                lm_model=lm(character~LCP,data=df_lm) %>% summary()
                test_result=lm_model$coefficients %>% .[nrow(.),]
            }
            return(test_result)
        })  %>% t() %>% data.frame(check.names=FALSE)
        
        lr_test[,'lr_pair']=rownames(lr_test)
        lr_test[,'LCP']=obs_LCP

        return(lr_test)

    }) %>% .[!sapply(.,is.null)] %>% dplyr::bind_rows()
    
    # test_result=test_result[!sapply(test_result,is.null)]
    # test_result=dplyr::bind_rows(test_result)
    
    if (nrow(test_result)==0){
        return(NULL)
    }
    
    return(test_result)
    
}
