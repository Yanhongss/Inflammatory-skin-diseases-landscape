library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(openxlsx)
library(dplyr)


infla<-getGmt("gene_infla.gmt")
infla_gene2<-read.gmt("gene_infla.gmt")
infla_gene2<-infla_gene2%>% 
  dplyr::filter(term!="edc")


##GSE85034 
GSE85034_sample_exp <- read.delim("/data/syh/science/pso/geo/GSE85034_GPL10558_sample_exp.txt", row.names=1)
GSE85034_sample <- read.delim("/data/syh/science/pso/geo/GSE85034_sample.txt", row.names=1)
dat <- as.matrix(GSE85034_sample_exp)
gsva_mat <- gsva(expr=dat, 
                 gset.idx.list=infla, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核
gsva_mat<-data.frame(t(gsva_mat))

write.xlsx(gsva_mat, file = "gsva_GSE85034.xlsx",rowNames=TRUE)

gsva_GSE85034<-cbind(rownames(gsva_mat),gsva_mat)
names(gsva_GSE85034)[1]<-c("GEO")

GSE85034<-GSE85034_sample %>% 
  dplyr::select(Group,treatment,timepoint)
GSE85034<-subset(GSE85034,timepoint!=c("WK 1NL"))
GSE85034<-subset(GSE85034,timepoint!=c("WK 16"))

table(GSE85034$timepoint)
time_id<-GSE85034$timepoint
table(time_id)
time_id<-gsub("NL", "NL_Time0", time_id)
time_id<-gsub("LS", "LS_Time0", time_id)
time_id<-gsub("WK 1", "LS_Wk1", time_id)
time_id<-gsub("WK 1LS", "LS_Wk1", time_id)
time_id<-gsub("WK2", "LS_Wk2", time_id)
time_id<-gsub("WK 2", "LS_Wk2", time_id)
time_id<-gsub("WK 4", "LS_Wk4", time_id)

GSE85034$time_id<-time_id
GSE85034 <- cbind(rownames(GSE85034),GSE85034)
names(GSE85034)[1]<-c("GEO")
GSE85034<-merge(GSE85034,gsva_GSE85034,by="GEO")
write.xlsx(GSE85034, file = "GSE85034.xlsx",rowNames=TRUE)

GSE85034_subset<-subset(GSE85034,time_id=="LS_Time0")
write.xlsx(GSE85034_subset, file = "GSE85034_subset.xlsx",rowNames=TRUE)

##GSE130588 
GSE130588_sample_exp <- read.delim("/data/syh/AA_new/GSE130588_GPL570_sample_rawexp.txt", row.names=1)
GSE130588_sample <- read.delim("/data/syh/AA_new/GSE130588_sample_R_NR.txt", row.names=1)
dat <- as.matrix(GSE130588_sample_exp)
gsva_mat <- gsva(expr=dat, 
                 gset.idx.list=infla, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核
gsva_mat<-data.frame(t(gsva_mat))

write.xlsx(gsva_mat, file = "gsva_GSE130588.xlsx",rowNames=TRUE)

gsva_GSE130588<-cbind(rownames(gsva_mat),gsva_mat)
names(gsva_GSE130588)[1]<-c("GEO")

GSE130588<-GSE130588_sample %>% 
  dplyr::select(time,tissue,responder)
GSE130588<-subset(GSE130588,time==c("Week0"))
GSE130588<-subset(GSE130588,tissue==c("LS"))
GSE130588<-subset(GSE130588,responder!=c("Placebo"))
GSE130588<-subset(GSE130588,responder!=c("Delete"))

GSE130588 <- cbind(rownames(GSE130588),GSE130588)
names(GSE130588)[1]<-c("GEO")

GSE130588<-merge(GSE130588,gsva_GSE130588,by="GEO")
write.xlsx(GSE130588, file = "GSE130588_subset.xlsx",rowNames=TRUE)

GSE130588<-merge(GSE130588,gsva_GSE130588,by="GEO")
write.xlsx(GSE130588, file = "GSE130588.xlsx",rowNames=TRUE)
