CNVestimate_ST <- function(seurat_obj,
                           deconv_result,MalType, # parameters used to estimate non-malignant cell type propotion
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
