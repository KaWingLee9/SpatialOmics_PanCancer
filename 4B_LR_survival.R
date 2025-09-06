library(dplyr)
library(survival)
library(maxstat)
library(survminer)

lr_ls=readRDS('./CCI_CellChatDB.rds')
lr_names=names(lr_ls)

clinical_data=read.csv('../8_bulk_seq/TCGA_clinical_data.csv',row.names=1,check.names=FALSE)

test_result=lapply(c('BLCA','BRCA','CESC','CHOL','COAD','ESCA','GBM','HNSC','KICH',
       'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','PRAD','READ',
       'SKCM','STAD','TGCT','THCA','UCEC','UVM'),function(tumor_type){
    
    tpm_mat=read.table(paste0('../8_bulk_seq/TCGA_TPM/TCGA-',tumor_type,'_mRNA_TPM.txt'),check.names=FALSE) %>% t() %>% data.frame(check.names=FALSE)
    selected_samples=strsplit(rownames(tpm_mat),split='-') %>% sapply(function(x) {grepl('01',x[4])})
    tpm_mat=tpm_mat[selected_samples,]
    tpm_mat=tpm_mat %>% .[,colSums(.)!=0]
    tpm_mat=log(tpm_mat+1)
    
    df_tumor=parallel::mclapply(lr_names,function(lr){

        lr_genes=lr_ls[[lr]] %>% unlist()
        
        if (! lr_genes %in% colnames(tpm_mat) %>% all()){
            return(NULL)
        }

        df=tpm_mat[,lr_genes]
#         df=df %>% apply( 2,function(x) {(x-min(x))/(max(x)-min(x))} )

        df=data.frame(df)

        colnames(df)=paste0('g',1:length(lr_genes))

        df[,'product']=df %>% apply(1,prod)
#         df[,'sum']=df %>% apply(1,sum)

        df['bcr_patient_barcode']=gsub('-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+$','',rownames(df))
        df=dplyr::left_join(df,clinical_data[,c('bcr_patient_barcode','OS','OS.time',
                                                'age_at_initial_pathologic_diagnosis','clinical_stage')],by='bcr_patient_barcode')
        
        res=try({df_grouped=surv_cutpoint(df,time='OS.time',event='OS',variables=c('product')) %>% surv_categorize()},silent=TRUE)
        if (inherits(res,'try-error')) {return(NULL)}
        group=rep('group2',times=nrow(df_grouped))
        group[(df_grouped[,'product']=='high')]='group1'
        df_grouped['group']=factor(group)
        res=try({cox.model=coxph(Surv(OS.time,OS)~relevel(group,'group2'),data=df_grouped)},silent=TRUE)
        if (inherits(res,'try-error')) {return(NULL)}
        cox.model.result=summary(cox.model)
        cox.model.result.df=cox.model.result$coefficients
        df.result=data.frame(HR=cox.model.result.df['relevel(group, "group2")group1','exp(coef)'],
                             Pval=cox.model.result.df['relevel(group, "group2")group1','Pr(>|z|)'])
        df.result[,'lr_pair']=lr
        
        return(df.result)
        
    },mc.cores=15) %>% .[! sapply(.,is.null) ] %>% dplyr::bind_rows()
    
    df_tumor[,'tumor_type']=tumor_type
    print(tumor_type)
    return(df_tumor)
    
}) %>% dplyr::bind_rows()