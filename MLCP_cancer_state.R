library(Seurat)
library(irGSEA)
library(dplyr)

# gene set enrichment for spots
# seurat_obj: Seurat class of ST data
# gene_ls: list of gene sets (e.g. .rds files in `gene_set directory`)
seurat_obj <- irGSEA.score(object=seurat_obj,assay="RNA",slot="data",
                           custom=TRUE,geneset=gene_ls,msigdb=FALSE, 
                           method="AUCell", aucell.MaxRank=ceiling(0.2*nrow(seurat_obj)),
                           seeds=123,ncores=40)


# differential cancer state of tumor cells between MLCP
# sum of the covariate_cell_type could not be 1 to avoid multicollinearity
LCPModuleComparison <- function(character_mat,LCP_label,cell_type_comp,
         ref_LCP,obs_LCP,
         covariate_cell_type=NULL){
    
    ref_spot=filter(LCP_label,cluster_result==ref_LCP) %>% rownames()
    obs_spot=filter(LCP_label,cluster_result==obs_LCP) %>% rownames()

    ref_character=character_mat[ref_spot,]
    obs_character=character_mat[obs_spot,]

    ref_cell_type_comp=cell_type_comp[ref_spot,]
    obs_cell_type_comp=cell_type_comp[obs_spot,]

    character_mat_subset=rbind(ref_character,obs_character)
    cell_type_comp_subset=rbind(ref_cell_type_comp,obs_cell_type_comp)
    
    test_result=sapply(character_mat_subset,function(x){
        
        df_lm=data.frame(character=x)
        df_lm[,'LCP']=rep(c(ref_LCP,obs_LCP),times=c(length(ref_spot),length(obs_spot))) %>% 
            factor(levels=c(ref_LCP,obs_LCP))
        if (!is.null(covariate_cell_type)){
            df_lm=data.frame(df_lm,cell_type_comp_subset[,covariate_cell_type],check.names=FALSE)
            colnames(df_lm)[3:(ncol(df_lm))]=covariate_cell_type
            covariate_cell_type=sprintf(ifelse(grepl(" ", covariate_cell_type), "`%s`", "%s"), covariate_cell_type)
        }
        lm_model=lm(paste('character~',paste(c(covariate_cell_type,'LCP'),collapse='+')),data=df_lm) %>% summary()

        lm_model$coefficients %>% .[nrow(.),]
    }) %>% t() %>% data.frame(check.names=FALSE)
    
    test_result[,'p.adjust']=p.adjust(test_result[,4])
    return(test_result)
}
