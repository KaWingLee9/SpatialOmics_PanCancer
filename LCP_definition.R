# definition of MLCP and SLCP

## MLCP
### select malignant spots
malignant_spot <- filter(infercnv_label,CNVlabel=='Malignant') %>% rownames()
malignant_type <- spot_type[malignant_spot,]
sample_name <- gsub('_[ACTG]+-1','',rownames(malignant_type))
### remvoe NK, pDC, Mast cell for too low proportion
malignant_type=malignant_type[,c('Parachymal cell','Fibroblast','Endothelial cell','CD4 T','CD8 T','B','Macrophage','cDC')]
malignant_type <- apply(malignant_type,1,function(x){x/sum(x)}) %>% t()
### cutoff set to limit parachymal cell percentage
spot_type_scaled <- malignant_type
spot_type_scaled[spot_type_scaled[,'Parachymal cell']>=0.5,'Parachymal cell'] <- 0.5
### scale data
spot_type_scaled <- scale(spot_type_scaled)
spot_type_scaled[spot_type_scaled>=2]=2


## SLCP
