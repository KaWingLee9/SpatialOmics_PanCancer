import numpy as np
import pandas as pd
import scanpy as sc
import bbknn
import SOAPy_st as sp

def Read_MERSCOPE(exp_file,meta_file,image_file=None,
                  fov_filter=None,gene_filter='Blank-'):
    
    exp_mat=pd.read_csv(exp_file,index_col=0)
    if gene_filter:
        exp_mat = exp_mat.drop(columns=[col for col in exp_mat.columns if gene_filter in col])
    
    meta_data=pd.read_csv(meta_file,index_col=0)
    meta_data=meta_data.loc[exp_mat.index,:]
    meta_data=meta_data.dropna(axis=1)
    
    # if not fov_filter:
        
    # if not image_file:
    # if not fov_filter:
    
    adata=sc.AnnData(exp_mat)
    adata.obs=meta_data
    adata.obs["total_counts"]=exp_mat.apply(np.sum,axis=1)
    adata.obs["total_features"]=exp_mat.apply(lambda x:sum(x!=0),axis=1)
    adata.obs.index=adata.obs.index.astype(str)
    adata.obsm['spatial']=np.array(meta_data[['center_x','center_y']])
    
    return(adata)

# cell type clustering of all cells
# within-sample nomalization and scale -> batch removal by bbknn
data_list=['HumanBreastCancerPatient1','HumanColonCancerPatient1','HumanColonCancerPatient2','HumanLiverCancerPatient1','HumanLiverCancerPatient2',
           'HumanLungCancerPatient1','HumanLungCancerPatient2','HumanOvarianCancerPatient1','HumanOvarianCancerPatient2Slice1',
           'HumanOvarianCancerPatient2Slice2','HumanOvarianCancerPatient2Slice3','HumanProstateCancerPatient1','HumanProstateCancerPatient2']

adata_list=[]
for x in data_list:
    
    exp_file=r'./MERSCOPE/'+x+'/'+x+'_cell_by_gene.csv'
    meta_file=r'./MERSCOPE/'+x+'/'+x+'_cell_metadata.csv'

    # read in data
    adata=Read_MERSCOPE(exp_file,meta_file)
    adata.obs['dataset']=x
    # qc
    sc.pp.calculate_qc_metrics(adata,inplace=True)
    sc.pp.filter_cells(adata,min_genes=100)
    sc.pp.filter_genes(adata,min_cells=5)
    # normalize
    sc.pp.normalize_total(adata)
    sc.pp.scale(adata,max_value=10)
    adata_list.append(adata)

adata=sc.concat(adata_list,label="batch",keys=data_list)
sc.tl.pca(adata)
bbknn.bbknn(adata,batch_key='batch')
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1)
adata.obs['louvain']=[ int(i) for i in adata.obs['louvain'] ]

# cell type clustering of TME cells
adata_TME=adata[adata.obs['CellType_Lv1']=='TME cell',:]
sc.tl.pca(adata_TME)
bbknn.bbknn(adata_TME,batch_key='batch')

sc.tl.umap(adata_TME)
sc.tl.louvain(adata_TME,resolution=1)

# combine meta data
df_1=adata_TME.obs[['dataset','CellType_Lv2']]
df_1=df_1.reset_index()
df_1['cell']=[ df_1.loc[i,'dataset']+'_'+df_1.loc[i,'cell'] for i in range(df_1.shape[0]) ]
df_1=df_1[['cell','CellType_Lv2']]
df_1.index=df_1['cell']

df_2=adata.obs
df_2=df_2.reset_index()
df_2['cell']=[ df_2.loc[i,'dataset']+'_'+df_2.loc[i,'cell'] for i in range(df_2.shape[0]) ]
df_2.index=df_2['cell']

x=[i in df_1.index for i in df_2.index]
df_2.loc[x,'CellType_Lv2']=df_1.loc[:,'CellType_Lv2']

y=list(pd.isna(df_2['CellType_Lv2']))
df_2['CellType_Lv2']=df_2['CellType_Lv2'].astype(str)
df_2.loc[y,'CellType_Lv2']='Tumor cell'

adata.obs=df_2

# dotplot
marker_genes=['EPCAM','CDH1','ERBB2','EGFR', # Epithelial cell
             'COL1A1','FN1','PDGFRB', # Fibroblast
             'PECAM1','CDH5','VWF', # Endothelial cell
             'CD3E','CD2','TRAC','CD4','CD8A','NKG7','GNLY', # T/NK
             'CD19','MS4A1','CD79A','XBP1','MZB1','FCRL5', # B/Plasma
             'LYZ','C1QC','HLA-DRA','XCR1','FCER1A','CD1C','LAMP3', # Macro/cDC
             'KIT','IRF4','CLEC4C','TCL1A' # Mast/pDC
             ]
sc.pl.dotplot(adata,marker_genes,groupby="CellType_Lv2",standard_scale='var',
              categories_order=['Tumor cell','Fibroblast','Endothelial cell','CD4 T cell','CD8 T cell','NK',
                                'B cell','Plasma cell','Macrophage','cDC','pDC','Mast cell','Undefined'])

# niche identification
sp.pp.make_network(adata,sample_key='dataset',method='knn',
                   cutoff=15,scale=1.0,cluster_key='CellType_Lv2')
sp.tl.get_c_niche(adata,k_max=30,celltype_key='CellType_Lv2',sample_key='dataset')

