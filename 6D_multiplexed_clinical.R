library(dplyr)
library(ggplot2)

library(survival)
library(maxstat)
library(survminer)

# survival analysis
spatial_score_df=read.csv(paste0(path,'RN_shuffle_except1_df.csv'),check.names=FALSE)
colnames(spatial_score_df)[1]='ImageID'

meta.data=read.csv(paste0(path,'meta_data.csv'),row.names=1,check.names=FALSE)
spatial_score_df=dplyr::left_join(spatial_score_df,meta.data,by=c('ImageID'='RoiID'))

spatial_score_df=spatial_score_df %>% filter(! is.na(OS))

OS.cutoff=2800

spatial_score_df[spatial_score_df[,'OS']>=OS.cutoff,'Ev.O']=0
spatial_score_df[spatial_score_df[,'OS']>=OS.cutoff,'OS']=OS.cutoff
spatial_score_df=spatial_score_df %>% group_by(Patient_ID,Age,Gender,Typ,Grade,Ev.O,OS,DFS,DX.name) %>% 
    summarise_at(vars(Cell_num,CT1_num,CT2_num,CT3_num,CT1_CT2,CT1_CT3),sum,na.rm=TRUE) %>% data.frame()

spatial_score_df['CT1_CT2_ratio']=spatial_score_df['CT1_CT2']/spatial_score_df['Cell_num']
spatial_score_df['CT1_CT3_ratio']=spatial_score_df['CT1_CT3']/spatial_score_df['Cell_num']

survfit(Surv(OS,Ev.O)~CT1_CT2_group+CT1_CT3_group,data=spatial_score_df) %>%
  ggsurvplot(data=spatial_score_df,pval=TRUE,risk.table=TRUE)

spatial_score_df[,'CT1_CT2_group']=NA
spatial_score_df$CT1_CT2_group[spatial_score_df$CT1_CT2_ratio<= median(spatial_score_df$CT1_CT2_ratio,na.rm=TRUE)]="Low"
spatial_score_df$CT1_CT2_group[spatial_score_df$CT1_CT2_ratio>= median(spatial_score_df$CT1_CT2_ratio,na.rm=TRUE)]="High"

spatial_score_df[,'CT1_CT3_group']=NA
spatial_score_df$CT1_CT3_group[spatial_score_df$CT1_CT3_ratio<= median(spatial_score_df$CT1_CT3_ratio,na.rm=TRUE)]="Low"
spatial_score_df$CT1_CT3_group[spatial_score_df$CT1_CT3_ratio>= median(spatial_score_df$CT1_CT3_ratio,na.rm=TRUE)]="High"
p=survfit(Surv(OS,Ev.O)~CT1_CT2_group+CT1_CT3_group,data=spatial_score_df) %>%
  ggsurvplot(data=spatial_score_df,pval=TRUE,risk.table=FALSE,surv.median.line='none',break.x.by=700)
p

# response to immunotherapy

meta.data=read.csv(paste0(path,'ROI_meta_afterQC.csv'),row.names=1,
                   check.names=FALSE)
colnames(meta.data)[1]='ImageID'
meta.data$PatientID=gsub('_ROI[0-9]+$','',meta.data$ImageID)

spatial_score_df=read.csv(paste0(path,'RN_shuffle_except1_df.csv'),check.names=FALSE)
colnames(spatial_score_df)[1]='ImageID'

df=dplyr::left_join(spatial_score_df,meta.data,by=c('ImageID'))

df=df %>% filter(Cell_num>=100,CT1_num>=0) %>% 
    group_by(PatientID,response) %>%
    summarise(Cell_num=sum(Cell_num),CT1_num=sum(CT1_num),CT2_num=sum(CT2_num),CT3_num=sum(CT3_num),
              CT1_CT2=sum(CT1_CT2),CT1_CT3=sum(CT1_CT3))

df['CT1_CT2_ratio']=df['CT1_CT2']/df['Cell_num']
df['CT1_CT3_ratio']=df['CT1_CT3']/df['Cell_num']
df['ratio']=df['CT1_CT3']/df['CT1_CT2']
df[ df['ratio']>=30 ,'ratio']=30

df=reshape2::melt(df,id.vars=c('PatientID','response'),measure.vars=c('CT1_CT2_ratio','CT1_CT3_ratio'),
               variable.name='indicator',value.name='ratio')

df=dplyr::left_join(spatial_score_df,meta.data,by=c('ImageID'))

df=df %>% filter(Cell_num>=100,CT1_num>=0) %>% 
    group_by(PatientID,response) %>%
    summarise(Cell_num=sum(Cell_num),CT1_num=sum(CT1_num),CT2_num=sum(CT2_num),CT3_num=sum(CT3_num),
              CT1_CT2=sum(CT1_CT2),CT1_CT3=sum(CT1_CT3))
# df['CT1_CT2_ratio']=(df['CT1_num']/df['Cell_num'])*(df['CT1_CT2']/df['Cell_num'])
# df['CT1_CT3_ratio']=(df['CT1_num']/df['Cell_num'])*(df['CT1_CT3']/df['Cell_num'])
df['CT1_CT2_ratio']=df['CT1_CT2']/df['Cell_num']
df['CT1_CT3_ratio']=df['CT1_CT3']/df['Cell_num']
df['ratio']=df['CT1_CT3']/df['CT1_CT2']
df[ df['ratio']>=30 ,'ratio']=30

df=as.data.frame(df)

df[ df[,'CT1_CT2_ratio']>=median(df[,'CT1_CT2_ratio']) ,'CT1_CT2_group']='high'
df[ df[,'CT1_CT2_ratio']<median(df[,'CT1_CT2_ratio']) ,'CT1_CT2_group']='low'

df[ df[,'CT1_CT3_ratio']>=median(df[,'CT1_CT3_ratio']) ,'CT1_CT3_group']='high'
df[ df[,'CT1_CT3_ratio']<median(df[,'CT1_CT3_ratio']) ,'CT1_CT3_group']='low'

df[,'group']=paste0(df[,'CT1_CT2_group'],'-',df[,'CT1_CT3_group'])

df_type=table(df[,'group'],df[,'response']) %>% as.data.frame.array()
df=df_type %>% apply(2,function(x) {x/sum(x)}) %>% reshape2::melt(varnames=c('sample_type','pCR'))

