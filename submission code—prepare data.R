library(ggbeeswarm)
library(Seurat)
library(harmony)
library(Rcpp)
library(clusterProfiler)
library(GSEABase)
library(AUCell)
library("FactoMineR")
library("factoextra")
library(dplyr)
require(devtools)
library(ggplot2)


palette1<-c('#8dd3c7','#fee391','#bebada','#fb8072','#d9d9d9','#fdb462','#b3de69','#80b1d3',
            '#bc80bd','#ccebc5','#ffed6f','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
            '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
palette2<-c("#173271",	"#5c6f9b","#b9c1d4","#e7eaf0","#f3f6f4","#2878b5","#9ac9db") 
palette3<-c("#99bba7","#7ea441","#3a562b","#a87431",
            "#3d579f","#8fb2e1","#bfbdaa","#25143f","#856ba2",
            "#c1809f","#e3d1df","#bf8436","#efdb52","#2878b5","#9ac9db")
palette4<-c('#8dd3c7','#fee391','#bebada','#b3de69','#d9d9d9','#fdb462',
            "#99bba7","#7ea441","#3a562b","#a87431",
            "#3d579f","#8fb2e1","#bfbdaa","#25143f","#856ba2",
            "#c1809f","#e3d1df","#bf8436","#efdb52","#2878b5","#9ac9db")
palette5<-c('#8dd3c7','#fee391','#bebada',
            "#99bba7","#7ea441","#3a562b","#a87431",
            "#3d579f","#8fb2e1","#bfbdaa","#25143f","#856ba2",
            "#c1809f","#e3d1df","#bf8436","#efdb52","#2878b5","#9ac9db")


##GSE222840 PN and AD
PN1 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN1")
PN1<-CreateSeuratObject(counts =PN1,project = "GSE222840",min.cells = 3, min.features = 200)

PN2 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN2")
PN2<-CreateSeuratObject(counts =PN2,project = "GSE222840",min.cells = 3, min.features = 200)

PN3 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN3")
PN3<-CreateSeuratObject(counts =PN3,project = "GSE222840",min.cells = 3, min.features = 200)

PN4 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN4")
PN4<-CreateSeuratObject(counts =PN4,project = "GSE222840",min.cells = 3, min.features = 200)

PN5 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN5")
PN5<-CreateSeuratObject(counts =PN5,project = "GSE222840",min.cells = 3, min.features = 200)

PN6 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN6")
PN6<-CreateSeuratObject(counts =PN6,project = "GSE222840",min.cells = 3, min.features = 200)

PN7 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/PN7")
PN7<-CreateSeuratObject(counts =PN7,project = "GSE222840",min.cells = 3, min.features = 200)

AD1 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AD1")
AD1<-CreateSeuratObject(counts =AD1,project = "GSE222840",min.cells = 3, min.features = 200)

AD2 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AD2")
AD2<-CreateSeuratObject(counts =AD2,project = "GSE222840",min.cells = 3, min.features = 200)

AD3 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AD3")
AD3<-CreateSeuratObject(counts =AD3,project = "GSE222840",min.cells = 3, min.features = 200)

AD4 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AD4")
AD4<-CreateSeuratObject(counts =AD4,project = "GSE222840",min.cells = 3, min.features = 200)

AD5 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AD5")
AD5<-CreateSeuratObject(counts =AD5,project = "GSE222840",min.cells = 3, min.features = 200)

AD6 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE222840/AP")
AD6<-CreateSeuratObject(counts =AD6,project = "GSE222840",min.cells = 3, min.features = 200)


PN1@meta.data$group <- "PN"
PN2@meta.data$group <- "PN"
PN3@meta.data$group <- "PN"
PN4@meta.data$group <- "PN"
PN5@meta.data$group <- "PN"
PN6@meta.data$group <- "PN"
PN7@meta.data$group <- "PN"
AD1@meta.data$group <- "AD"
AD2@meta.data$group <- "AD"
AD3@meta.data$group <- "AD"
AD4@meta.data$group <- "AD"
AD5@meta.data$group <- "AD"
AD6@meta.data$group <- "AD"

PN1@meta.data$group2 <- "GSE222840"
PN2@meta.data$group2 <- "GSE222840"
PN3@meta.data$group2 <- "GSE222840"
PN4@meta.data$group2 <- "GSE222840"
PN5@meta.data$group2 <- "GSE222840"
PN6@meta.data$group2 <- "GSE222840"
PN7@meta.data$group2 <- "GSE222840"
AD1@meta.data$group2 <- "GSE222840"
AD2@meta.data$group2 <- "GSE222840"
AD3@meta.data$group2 <- "GSE222840"
AD4@meta.data$group2 <- "GSE222840"
AD5@meta.data$group2 <- "GSE222840"
AD6@meta.data$group2 <- "GSE222840"

PN1@meta.data$group3 <- "GSE222840_PN1"
PN2@meta.data$group3 <- "GSE222840_PN2"
PN3@meta.data$group3 <- "GSE222840_PN3"
PN4@meta.data$group3 <- "GSE222840_PN4"
PN5@meta.data$group3 <- "GSE222840_PN5"
PN6@meta.data$group3 <- "GSE222840_PN6"
PN7@meta.data$group3 <- "GSE222840_PN7"
AD1@meta.data$group3 <- "GSE222840_AD1"
AD2@meta.data$group3 <- "GSE222840_AD2"
AD3@meta.data$group3 <- "GSE222840_AD3"
AD4@meta.data$group3 <- "GSE222840_AD4"
AD5@meta.data$group3 <- "GSE222840_AD5"
AD6@meta.data$group3 <- "GSE222840_AD6"

PN1@meta.data$group4 <- "PN"
PN2@meta.data$group4 <- "PN"
PN3@meta.data$group4 <- "PN"
PN4@meta.data$group4 <- "PN"
PN5@meta.data$group4 <- "PN"
PN6@meta.data$group4 <- "PN"
PN7@meta.data$group4 <- "PN"
AD1@meta.data$group4 <- "AD"
AD2@meta.data$group4 <- "AD"
AD3@meta.data$group4 <- "AD"
AD4@meta.data$group4 <- "AD"
AD5@meta.data$group4 <- "AD"
AD6@meta.data$group4 <- "AD"

PN1@meta.data$group5 <- "PN1"
PN2@meta.data$group5 <- "PN2"
PN3@meta.data$group5 <- "PN3"
PN4@meta.data$group5 <- "PN4"
PN5@meta.data$group5 <- "PN5"
PN6@meta.data$group5 <- "PN6"
PN7@meta.data$group5 <- "PN7"
AD1@meta.data$group5 <- "AD1"
AD2@meta.data$group5 <- "AD2"
AD3@meta.data$group5 <- "AD3"
AD4@meta.data$group5 <- "AD4"
AD5@meta.data$group5 <- "AD5"
AD6@meta.data$group5 <- "AD6"

PN1[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN1[["RNA"]]))
PN2[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN2[["RNA"]]))
PN3[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN3[["RNA"]]))
PN4[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN4[["RNA"]]))
PN5[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN5[["RNA"]]))
PN6[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN6[["RNA"]]))
PN7[["RNA"]]@meta.features<- data.frame(row.names = rownames(PN7[["RNA"]]))
AD1[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD1[["RNA"]]))
AD2[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD2[["RNA"]]))
AD3[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD3[["RNA"]]))
AD4[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD4[["RNA"]]))
AD5[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD5[["RNA"]]))
AD6[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD6[["RNA"]]))

GSE222840<-merge(PN1,y=c(PN2,PN3,PN4,PN5,PN6,PN7,AD1,AD2,AD3,AD4,AD5,AD6))
GSE222840_nFeature<-as.numeric(GSE222840@meta.data$nFeature_RNA)
class(GSE222840_nFeature)
GSE222840_nCount<-GSE222840@meta.data$nCount_RNA
class(GSE222840_nCount)
min(GSE222840_nFeature)
max(GSE222840_nFeature)
VlnPlot(object = GSE222840, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1197.5-9177.5
GSE222840<-subset(GSE222840,subset = nFeature_RNA > 1197.5 & nFeature_RNA < 9177.5)
rm(PN1,PN2,PN3,PN4,PN5,PN6,PN7,AD1,AD2,AD3,AD4,AD5,AD6)
saveRDS(GSE222840, "GSE222840.rds")
table(GSE222840$group2)

##GSE212450 AA
AA2 <- Read10X(data.dir = "/data/syh/AA/AA/AA2")
AA2<-CreateSeuratObject(counts =AA2,project = "GSE212450",min.cells = 3, min.features = 200)

AA4 <- Read10X(data.dir = "/data/syh/AA/AA/AA4")
AA4<-CreateSeuratObject(counts =AA4,project = "GSE212450",min.cells = 3, min.features = 200)

AA7 <- Read10X(data.dir = "/data/syh/AA/AA/AA7")
AA7<-CreateSeuratObject(counts =AA7,project = "GSE212450",min.cells = 3, min.features = 200)

AA8 <- Read10X(data.dir = "/data/syh/AA/AA/AA8")
AA8<-CreateSeuratObject(counts =AA8,project = "GSE212450",min.cells = 3, min.features = 200)

PB1 <- Read10X(data.dir = "/data/syh/AA/AA/PB1")
PB1<-CreateSeuratObject(counts =PB1,project = "GSE212450",min.cells = 3, min.features = 200)

PB2 <- Read10X(data.dir = "/data/syh/AA/AA/PB2")
PB2<-CreateSeuratObject(counts =PB2,project = "GSE212450",min.cells = 3, min.features = 200)

PB3 <- Read10X(data.dir = "/data/syh/AA/AA/PB3")
PB3<-CreateSeuratObject(counts =PB3,project = "GSE212450",min.cells = 3, min.features = 200)

SD1 <- Read10X(data.dir = "/data/syh/AA/AA/SD1")
SD1<-CreateSeuratObject(counts =SD1,project = "GSE212450",min.cells = 3, min.features = 200)

SD2 <- Read10X(data.dir = "/data/syh/AA/AA/SD2")
SD2<-CreateSeuratObject(counts =SD2,project = "GSE212450",min.cells = 3, min.features = 200)

SD3 <- Read10X(data.dir = "/data/syh/AA/AA/SD3")
SD3<-CreateSeuratObject(counts =SD3,project = "GSE212450",min.cells = 3, min.features = 200)

AA2@meta.data$group <- "AA"
AA4@meta.data$group <- "AA"
AA7@meta.data$group <- "AA"
AA8@meta.data$group <- "AA"
PB1@meta.data$group <- "Normal"
PB2@meta.data$group <- "Normal"
PB3@meta.data$group <- "Normal"
SD1@meta.data$group <- "AA"
SD2@meta.data$group <- "AA"
SD3@meta.data$group <- "AA"

AA2@meta.data$group2 <- "GSE212450"
AA4@meta.data$group2 <- "GSE212450"
AA7@meta.data$group2 <- "GSE212450"
AA8@meta.data$group2 <- "GSE212450"
PB1@meta.data$group2 <- "GSE212450"
PB2@meta.data$group2 <- "GSE212450"
PB3@meta.data$group2 <- "GSE212450"
SD1@meta.data$group2 <- "GSE212450"
SD2@meta.data$group2 <- "GSE212450"
SD3@meta.data$group2 <- "GSE212450"

AA2@meta.data$group3 <- "GSE212450_AA2"
AA4@meta.data$group3 <- "GSE212450_AA4"
AA7@meta.data$group3 <- "GSE212450_AA7"
AA8@meta.data$group3 <- "GSE212450_AA8"
PB1@meta.data$group3 <- "GSE212450_PB1"
PB2@meta.data$group3 <- "GSE212450_PB2"
PB3@meta.data$group3 <- "GSE212450_PB3"
SD1@meta.data$group3 <- "GSE212450_SD1"
SD2@meta.data$group3 <- "GSE212450_SD2"
SD3@meta.data$group3 <- "GSE212450_SD3"

AA2@meta.data$group4 <- "AA"
AA4@meta.data$group4 <- "AA"
AA7@meta.data$group4 <- "AA"
AA8@meta.data$group4 <- "AA"
PB1@meta.data$group4 <- "Normal"
PB2@meta.data$group4 <- "Normal"
PB3@meta.data$group4 <- "Normal"
SD1@meta.data$group4 <- "AA_NL"
SD2@meta.data$group4 <- "AA_NL"
SD3@meta.data$group4 <- "AA_NL"

AA2@meta.data$group5 <- "AA1"
AA4@meta.data$group5 <- "AA2"
AA7@meta.data$group5 <- "AA3"
AA8@meta.data$group5 <- "AA4"
PB1@meta.data$group5 <- "Normal1"
PB2@meta.data$group5 <- "Normal2"
PB3@meta.data$group5 <- "Normal3"
SD1@meta.data$group5 <- "AA_NL1"
SD2@meta.data$group5 <- "AA_NL2"
SD3@meta.data$group5 <- "AA_NL3"

AA2[["RNA"]]@meta.features<- data.frame(row.names = rownames(AA2[["RNA"]]))
AA4[["RNA"]]@meta.features<- data.frame(row.names = rownames(AA4[["RNA"]]))
AA7[["RNA"]]@meta.features<- data.frame(row.names = rownames(AA7[["RNA"]]))
AA8[["RNA"]]@meta.features<- data.frame(row.names = rownames(AA8[["RNA"]]))
PB1[["RNA"]]@meta.features<- data.frame(row.names = rownames(PB1[["RNA"]]))
PB2[["RNA"]]@meta.features<- data.frame(row.names = rownames(PB2[["RNA"]]))
PB3[["RNA"]]@meta.features<- data.frame(row.names = rownames(PB3[["RNA"]]))
SD1[["RNA"]]@meta.features<- data.frame(row.names = rownames(SD1[["RNA"]]))
SD2[["RNA"]]@meta.features<- data.frame(row.names = rownames(SD2[["RNA"]]))
SD3[["RNA"]]@meta.features<- data.frame(row.names = rownames(SD3[["RNA"]]))

AA<-merge(AA2,y=c(AA4,AA7,AA8,PB1,PB2,PB3,SD1,SD2,SD3))
AA_nFeature<-as.numeric(AA@meta.data$nFeature_RNA)
class(AA_nFeature)
AA_nCount<-AA@meta.data$nCount_RNA
class(AA_nCount)
min(AA_nFeature)
max(AA_nFeature)
VlnPlot(object = AA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1012-7572
AA<-subset(AA,subset = nFeature_RNA > 1012 & nFeature_RNA < 7572)
rm(AA2,AA4,AA7,AA8,PB1,PB2,PB3,SD1,SD2,SD3)
saveRDS(AA, "GSE212450_AA.rds")

##GSE175990 HS
HS_1_1 <- Read10X(data.dir = "/data/syh/AA/HS/1_GSM5352392_HS_1")
HS_1_1<-CreateSeuratObject(counts =HS_1_1,project = "GSE175990",min.cells = 3, min.features = 200)

HS_1_2 <- Read10X(data.dir = "/data/syh/AA/HS/1_GSM5352393_HS_2")
HS_1_2<-CreateSeuratObject(counts =HS_1_2,project = "GSE175990",min.cells = 3, min.features = 200)

HS_1_3 <- Read10X(data.dir = "/data/syh/AA/HS/1_GSM5352394_HS_3")
HS_1_3<-CreateSeuratObject(counts =HS_1_3,project = "GSE175990",min.cells = 3, min.features = 200)

HC_1 <- Read10X(data.dir = "/data/syh/AA/HS/1_GSM5352395_HC")
HC_1<-CreateSeuratObject(counts =HC_1,project = "GSE175990",min.cells = 3, min.features = 200)

HS_1_1@meta.data$group <- "HS"
HS_1_2@meta.data$group <- "HS"
HS_1_3@meta.data$group <- "HS"
HC_1@meta.data$group <- "normal"

HS_1_1@meta.data$group2 <- "project_GSE175990"
HS_1_2@meta.data$group2 <- "project_GSE175990"
HS_1_3@meta.data$group2 <- "project_GSE175990"
HC_1@meta.data$group2 <- "project_GSE175990"

HS_1_1@meta.data$group3 <- "project_GSE175990-HS_1_1"
HS_1_2@meta.data$group3 <- "project_GSE175990-HS_1_2"
HS_1_3@meta.data$group3 <- "project_GSE175990-HS_1_3"
HC_1@meta.data$group3 <- "project_GSE175990-normal"

HS_1_1@meta.data$group4 <- "HS"
HS_1_2@meta.data$group4 <- "HS"
HS_1_3@meta.data$group4 <- "HS"
HC_1@meta.data$group4 <- "Normal"

HS_1_1@meta.data$group5 <- "HS1"
HS_1_2@meta.data$group5 <- "HS2"
HS_1_3@meta.data$group5 <- "HS3"
HC_1@meta.data$group5 <- "Normal4"

HS_1_1[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_1_1[["RNA"]]))
HS_1_2[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_1_2[["RNA"]]))
HS_1_3[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_1_3[["RNA"]]))
HC_1[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC_1[["RNA"]]))

HS_1<-merge(HS_1_1,y=c(HS_1_2,HS_1_3,HC_1))
HS_1_nFeature<-as.numeric(HS_1@meta.data$nFeature_RNA)
class(HS_1_nFeature)
HS_1_nCount<-HS_1@meta.data$nCount_RNA
class(HS_1_nCount)
min(HS_1_nFeature)
max(HS_1_nFeature)
VlnPlot(object = HS_1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1105.5-8373.5
HS_1<-subset(HS_1,subset = nFeature_RNA > 1105.5 & nFeature_RNA < 8373.5)
rm(HS_1_1,HS_1_2,HS_1_3,HC_1)
saveRDS(HS_1, "GSE175990_HS.rds")

##GSE154775 HS
HS_2_1<- Read10X(data.dir = "/data/syh/AA/HS/GSM4679492_HS1")
HS_2_1<-CreateSeuratObject(counts =HS_2_1,project = "GSE154775",min.cells = 3, min.features = 200)

HS_2_2<- Read10X(data.dir = "/data/syh/AA/HS/GSM4679493_HS2")
HS_2_2<-CreateSeuratObject(counts =HS_2_2,project = "GSE154775",min.cells = 3, min.features = 200)

HS_2_3<- Read10X(data.dir = "/data/syh/AA/HS/GSM4679500_HS_07_2018-06-29")
HS_2_3<-CreateSeuratObject(counts =HS_2_3,project = "GSE154775",min.cells = 3, min.features = 200)

HS_2_1@meta.data$group <- "HS"
HS_2_2@meta.data$group <- "HS"
HS_2_3@meta.data$group <- "HS"

HS_2_1@meta.data$group2 <- "GSE154775"
HS_2_2@meta.data$group2 <- "GSE154775"
HS_2_3@meta.data$group2 <- "GSE154775"

HS_2_1@meta.data$group3 <- "GSE154775_HS_2_1"
HS_2_2@meta.data$group3 <- "GSE154775_HS_2_2"
HS_2_3@meta.data$group3 <- "GSE154775_HS_2_3"

HS_2_1@meta.data$group4 <- "HS"
HS_2_2@meta.data$group4 <- "HS"
HS_2_3@meta.data$group4 <- "HS"

HS_2_1@meta.data$group5 <- "HS4"
HS_2_2@meta.data$group5 <- "HS5"
HS_2_3@meta.data$group5 <- "HS6"

HS_2_1[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_2_1[["RNA"]]))
HS_2_2[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_2_2[["RNA"]]))
HS_2_3[["RNA"]]@meta.features<- data.frame(row.names = rownames(HS_2_3[["RNA"]]))

HS_2<-merge(HS_2_1,y=c(HS_2_2,HS_2_3))
HS_2_nFeature<-as.numeric(HS_2@meta.data$nFeature_RNA)
class(HS_2_nFeature)
HS_2_nCount<-HS_2@meta.data$nCount_RNA
class(HS_2_nCount)
min(HS_2_nFeature)
max(HS_2_nFeature)
VlnPlot(object = HS_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1177.3-6499.7
HS_2<-subset(HS_2,subset = nFeature_RNA > 1177.3 & nFeature_RNA < 6499.7)
rm(HS_2_1,HS_2_2,HS_2_3)
saveRDS(HS_2, "GSE154775_HS.rds")

#GSE162183 PSO
GSE162183 <- read.table('/data/syh/AA/pso_AD_supplement/GSE162183/GSE162183_Raw_gene_counts_matrix.tab',sep='\t')
rownames(GSE162183)<-GSE162183$V1
ROWNAMES<-GSE162183$V1
GSE162183<-GSE162183[,-1]
COLNAMES<- GSE162183[1,]
colnames(GSE162183)<-COLNAMES
colnames(GSE162183)
rownames(GSE162183)
GSE162183<-GSE162183[-1,]

Ctrl1<- grep(pattern ="^Ctrl1",colnames
             (GSE162183))
Ctrl1_count<-GSE162183[,Ctrl1]
rownames(Ctrl1_count)<-rownames(GSE162183)
colnames(Ctrl1_count)
Ctrl1<-CreateSeuratObject(counts =Ctrl1_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Ctrl1_count)

Ctrl2<- grep(pattern ="^Ctrl2",colnames
             (GSE162183))
Ctrl2_count<-GSE162183[,Ctrl2]
Ctrl2<-CreateSeuratObject(counts =Ctrl2_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Ctrl2_count)

Ctrl3<- grep(pattern ="^Ctrl3",colnames
             (GSE162183))
Ctrl3_count<-GSE162183[,Ctrl3]
Ctrl3<-CreateSeuratObject(counts =Ctrl3_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Ctrl3_count)

Psor1<- grep(pattern ="^Psor1",colnames
             (GSE162183))
Psor1_count<-GSE162183[,Psor1]
Psor1<-CreateSeuratObject(counts =Psor1_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Psor1_count)

Psor2<- grep(pattern ="^Psor2",colnames
             (GSE162183))
Psor2_count<-GSE162183[,Psor2]
Psor2<-CreateSeuratObject(counts =Psor2_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Psor2_count)

Psor3<- grep(pattern ="^Psor3",colnames
             (GSE162183))
Psor3_count<-GSE162183[,Psor3]
Psor3<-CreateSeuratObject(counts =Psor3_count,project = "GSE162183",min.cells = 3, min.features = 200)
rm(Psor3_count)

Ctrl1@meta.data$group <- "Normal"
Ctrl2@meta.data$group <- "Normal"
Ctrl3@meta.data$group <- "Normal"
Psor1@meta.data$group <- "Psoriasis"
Psor2@meta.data$group <- "Psoriasis"
Psor3@meta.data$group <- "Psoriasis"

Ctrl1@meta.data$group2 <- "GSE162183"
Ctrl2@meta.data$group2 <- "GSE162183"
Ctrl3@meta.data$group2 <- "GSE162183"
Psor1@meta.data$group2 <- "GSE162183"
Psor2@meta.data$group2 <- "GSE162183"
Psor3@meta.data$group2 <- "GSE162183"

Ctrl1@meta.data$group3 <- "GSE162183_Ctrl1"
Ctrl2@meta.data$group3 <- "GSE162183_Ctrl2"
Ctrl3@meta.data$group3 <- "GSE162183_Ctrl3"
Psor1@meta.data$group3 <- "GSE162183_Psor1"
Psor2@meta.data$group3 <- "GSE162183_Psor2"
Psor3@meta.data$group3 <- "GSE162183_Psor3"

Ctrl1@meta.data$group4 <- "Normal"
Ctrl2@meta.data$group4 <- "Normal"
Ctrl3@meta.data$group4 <- "Normal"
Psor1@meta.data$group4 <- "Psoriasis"
Psor2@meta.data$group4 <- "Psoriasis"
Psor3@meta.data$group4 <- "Psoriasis"

Ctrl1@meta.data$group5 <- "Normal5"
Ctrl2@meta.data$group5 <- "Normal6"
Ctrl3@meta.data$group5 <- "Normal7"
Psor1@meta.data$group5 <- "Pso1"
Psor2@meta.data$group5 <- "Pso2"
Psor3@meta.data$group5 <- "Pso3"

Ctrl1[["RNA"]]@meta.features<- data.frame(row.names = rownames(Ctrl1[["RNA"]]))
Ctrl2[["RNA"]]@meta.features<- data.frame(row.names = rownames(Ctrl2[["RNA"]]))
Ctrl3[["RNA"]]@meta.features<- data.frame(row.names = rownames(Ctrl3[["RNA"]]))
Psor1[["RNA"]]@meta.features<- data.frame(row.names = rownames(Psor1[["RNA"]]))
Psor2[["RNA"]]@meta.features<- data.frame(row.names = rownames(Psor2[["RNA"]]))
Psor3[["RNA"]]@meta.features<- data.frame(row.names = rownames(Psor3[["RNA"]]))

rm(GSE162183)
GSE162183<-merge(Ctrl1,y=c(Ctrl2,Ctrl3,Psor1,Psor2,Psor3))
GSE162183_nFeature<-as.numeric(GSE162183@meta.data$nFeature_RNA)
class(GSE162183_nFeature)
GSE162183_nCount<-GSE162183@meta.data$nCount_RNA
class(GSE162183_nCount)
min(GSE162183_nFeature)
max(GSE162183_nFeature)
VlnPlot(object = HS_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1387.5-9039.5
GSE162183<-subset(GSE162183,subset = nFeature_RNA > 1387.5 & nFeature_RNA < 9039.5)
rm(Ctrl1,Ctrl2,Ctrl3,Psor1,Psor2,Psor3)
saveRDS(GSE162183, "GSE162183_PSO.rds")

##GSE153760 AD
HC1 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC1")
HC1<-CreateSeuratObject(counts =HC1,project = "GSE153760",min.cells = 3, min.features = 200)

HC2 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC2")
HC2<-CreateSeuratObject(counts =HC2,project = "GSE153760",min.cells = 3, min.features = 200)

HC3 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC3")
HC3<-CreateSeuratObject(counts =HC3,project = "GSE153760",min.cells = 3, min.features = 200)

HC4 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC4")
HC4<-CreateSeuratObject(counts =HC4,project = "GSE153760",min.cells = 3, min.features = 200)

HC5 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC5")
HC5<-CreateSeuratObject(counts =HC5,project = "GSE153760",min.cells = 3, min.features = 200)

HC6 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC6")
HC6<-CreateSeuratObject(counts =HC6,project = "GSE153760",min.cells = 3, min.features = 200)

HC7 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/HC7")
HC7<-CreateSeuratObject(counts =HC7,project = "GSE153760",min.cells = 3, min.features = 200)

AD1 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD1")
AD1<-CreateSeuratObject(counts =AD1,project = "GSE153760",min.cells = 3, min.features = 200)

AD2 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD2")
AD2<-CreateSeuratObject(counts =AD2,project = "GSE153760",min.cells = 3, min.features = 200)

AD3 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD3")
AD3<-CreateSeuratObject(counts =AD3,project = "GSE153760",min.cells = 3, min.features = 200)

AD4 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD4")
AD4<-CreateSeuratObject(counts =AD4,project = "GSE153760",min.cells = 3, min.features = 200)

AD5 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD5")
AD5<-CreateSeuratObject(counts =AD5,project = "GSE153760",min.cells = 3, min.features = 200)

AD6 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD6")
AD6<-CreateSeuratObject(counts =AD6,project = "GSE153760",min.cells = 3, min.features = 200)

AD7 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD7")
AD7<-CreateSeuratObject(counts =AD7,project = "GSE153760",min.cells = 3, min.features = 200)

AD8 <- Read10X(data.dir = "/data/syh/AA/pso_AD_supplement/GSE153760/AD8")
AD8<-CreateSeuratObject(counts =AD8,project = "GSE153760",min.cells = 3, min.features = 200)

HC1@meta.data$group <- "Normal"
HC2@meta.data$group <- "Normal"
HC3@meta.data$group <- "Normal"
HC4@meta.data$group <- "Normal"
HC5@meta.data$group <- "Normal"
HC6@meta.data$group <- "Normal"
HC7@meta.data$group <- "Normal"
AD1@meta.data$group <- "AD"
AD2@meta.data$group <- "AD"
AD3@meta.data$group <- "AD"
AD4@meta.data$group <- "AD"
AD5@meta.data$group <- "AD"
AD6@meta.data$group <- "AD"
AD7@meta.data$group <- "AD"
AD8@meta.data$group <- "AD"

HC1@meta.data$group2 <- "GSE153760"
HC2@meta.data$group2 <- "GSE153760"
HC3@meta.data$group2 <- "GSE153760"
HC4@meta.data$group2 <- "GSE153760"
HC5@meta.data$group2 <- "GSE153760"
HC6@meta.data$group2 <- "GSE153760"
HC7@meta.data$group2 <- "GSE153760"
AD1@meta.data$group2 <- "GSE153760"
AD2@meta.data$group2 <- "GSE153760"
AD3@meta.data$group2 <- "GSE153760"
AD4@meta.data$group2 <- "GSE153760"
AD5@meta.data$group2 <- "GSE153760"
AD6@meta.data$group2 <- "GSE153760"
AD7@meta.data$group2 <- "GSE153760"
AD8@meta.data$group2 <- "GSE153760"

HC1@meta.data$group3 <- "GSE153760_HC1"
HC2@meta.data$group3 <- "GSE153760_HC2"
HC3@meta.data$group3 <- "GSE153760_HC3"
HC4@meta.data$group3 <- "GSE153760_HC4"
HC5@meta.data$group3 <- "GSE153760_HC5"
HC6@meta.data$group3 <- "GSE153760_HC6"
HC7@meta.data$group3 <- "GSE153760_HC7"
AD1@meta.data$group3 <- "GSE153760_AD1"
AD2@meta.data$group3 <- "GSE153760_AD2"
AD3@meta.data$group3 <- "GSE153760_AD3"
AD4@meta.data$group3 <- "GSE153760_AD4"
AD5@meta.data$group3 <- "GSE153760_AD5"
AD6@meta.data$group3 <- "GSE153760_AD6"
AD7@meta.data$group3 <- "GSE153760_AD7"
AD8@meta.data$group3 <- "GSE153760_AD8"

HC1@meta.data$group4 <- "Normal"
HC2@meta.data$group4 <- "Normal"
HC3@meta.data$group4 <- "Normal"
HC4@meta.data$group4 <- "Normal"
HC5@meta.data$group4 <- "Normal"
HC6@meta.data$group4 <- "Normal"
HC7@meta.data$group4 <- "Normal"
AD1@meta.data$group4 <- "AD"
AD2@meta.data$group4 <- "AD"
AD3@meta.data$group4 <- "AD"
AD4@meta.data$group4 <- "AD"
AD5@meta.data$group4 <- "AD"
AD6@meta.data$group4 <- "AD"
AD7@meta.data$group4 <- "AD"
AD8@meta.data$group4 <- "AD"

HC1@meta.data$group5 <- "Normal8"
HC2@meta.data$group5 <- "Normal9"
HC3@meta.data$group5 <- "Normal10"
HC4@meta.data$group5 <- "Normal11"
HC5@meta.data$group5 <- "Normal12"
HC6@meta.data$group5 <- "Normal13"
HC7@meta.data$group5 <- "Normal14"
AD1@meta.data$group5 <- "AD7"
AD2@meta.data$group5 <- "AD8"
AD3@meta.data$group5 <- "AD9"
AD4@meta.data$group5 <- "AD10"
AD5@meta.data$group5 <- "AD11"
AD6@meta.data$group5 <- "AD12"
AD7@meta.data$group5 <- "AD13"
AD8@meta.data$group5 <- "AD14"

HC1[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC1[["RNA"]]))
HC2[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC2[["RNA"]]))
HC3[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC3[["RNA"]]))
HC4[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC4[["RNA"]]))
HC5[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC5[["RNA"]]))
HC6[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC6[["RNA"]]))
HC7[["RNA"]]@meta.features<- data.frame(row.names = rownames(HC7[["RNA"]]))

AD1[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD1[["RNA"]]))
AD2[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD2[["RNA"]]))
AD3[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD3[["RNA"]]))
AD4[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD4[["RNA"]]))
AD5[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD5[["RNA"]]))
AD6[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD6[["RNA"]]))
AD7[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD7[["RNA"]]))
AD8[["RNA"]]@meta.features<- data.frame(row.names = rownames(AD8[["RNA"]]))

GSE153760<-merge(HC1,y=c(HC2,HC3,HC4,HC5,HC6,HC7,AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8))
GSE153760_nFeature<-as.numeric(GSE153760@meta.data$nFeature_RNA)
class(GSE153760_nFeature)
GSE153760_nCount<-GSE153760@meta.data$nCount_RNA
class(GSE153760_nCount)
min(GSE153760_nFeature)
max(GSE153760_nFeature)
VlnPlot(object = GSE153760, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1062.8-8173.2
GSE153760<-subset(GSE153760,subset = nFeature_RNA > 1062.8 & nFeature_RNA < 8173.2)
rm(HC1,HC2,HC3,HC4,HC5,HC6,HC7,AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
rm(COLNAMES)

saveRDS(GSE153760, "GSE153760_AD.rds")

##GSE150672 PSO
load("/data/scRNA_mix/zrh_mix/rawdata_GSE150672_new.Rda")

##GSE147424 AD
load("/data/scRNA_mix/zrh_mix/rawdata_GSE147424_new.Rda")

##filter
##GSE150672
Normal1@meta.data$group4 <- "Normal"
Normal2@meta.data$group4 <- "Normal"
Normal3@meta.data$group4 <- "Normal"
psoriasis1@meta.data$group4 <- "Psoriasis"
psoriasis2@meta.data$group4 <- "Psoriasis"
psoriasis3@meta.data$group4 <- "Psoriasis"
psoriasis4@meta.data$group4 <- "Psoriasis"
psoriasis5@meta.data$group4 <- "Psoriasis"

Normal1@meta.data$group5 <- "Normal15"
Normal2@meta.data$group5 <- "Normal16"
Normal3@meta.data$group5 <- "Normal17"
psoriasis1@meta.data$group5 <- "Pso4"
psoriasis2@meta.data$group5 <- "Pso5"
psoriasis3@meta.data$group5 <- "Pso6"
psoriasis4@meta.data$group5 <- "Pso7"
psoriasis5@meta.data$group5 <- "Pso8"

GSE150672<-merge(Normal1,y=c(Normal2,Normal3,psoriasis1,psoriasis2,psoriasis3,psoriasis4,psoriasis5))
VlnPlot(object = GSE150672, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
GSE150672_nFeature<-as.numeric(GSE150672@meta.data$nFeature_RNA)
class(GSE150672_nFeature)
GSE150672_nCount<-GSE150672@meta.data$nCount_RNA
class(GSE150672_nCount)
mean_GSE150672_nFeature<-mean(GSE150672_nFeature)
sd_GSE150672_nFeature<-sd(GSE150672_nFeature)
min(GSE150672_nFeature)
max(GSE150672_nFeature)
## 10%-90% 969.8-6592.2
GSE150672<-subset(GSE150672,subset = nFeature_RNA > 969.8 & nFeature_RNA < 6592.2)
rm(Normal1,Normal2,Normal3,psoriasis1,psoriasis2,psoriasis3,psoriasis4,psoriasis5)
saveRDS(GSE150672, "GSE150672_PSO.rds")


##GSE147424
## 10%-90% 1173.8-3188.2
fibad<-subset(fibad,subset = nFeature_RNA > 1146.3 & nFeature_RNA < 8748.7)#这步完成于上一次保存前
fibad <- readRDS("/data/syh/AA/fibad.rds")
table(fibad$group3)
fibad.list <- SplitObject(fibad, split.by = "group3")
fib_ad1<-fibad.list$GSE147424_ad1
fib_ad2<-fibad.list$GSE147424_ad2
fib_ad3<-fibad.list$GSE147424_ad3
fib_ad4<-fibad.list$GSE147424_ad4
fib_n1<-fibad.list$GSE147424_n1
fib_n2<-fibad.list$GSE147424_n2
fib_n3<-fibad.list$GSE147424_n3
fib_n4<-fibad.list$GSE147424_n4
fib_n5<-fibad.list$GSE147424_n5
fib_n6<-fibad.list$GSE147424_n6
fib_n7<-fibad.list$GSE147424_n7
fib_n8<-fibad.list$GSE147424_n8
fib_NL1<-fibad.list$GSE147424_NL1
fib_NL2<-fibad.list$GSE147424_NL2
fib_NL3<-fibad.list$GSE147424_NL3
fib_NL4<-fibad.list$GSE147424_NL4
rm(fibad.list)

fib_ad1@meta.data$group4 <- "AD"
fib_ad2@meta.data$group4 <- "AD"
fib_ad3@meta.data$group4 <- "AD"
fib_ad4@meta.data$group4 <- "AD"
fib_n1@meta.data$group4 <- "Normal"
fib_n2@meta.data$group4 <- "Normal"
fib_n3@meta.data$group4 <- "Normal"
fib_n4@meta.data$group4 <- "Normal"
fib_n5@meta.data$group4 <- "Normal"
fib_n6@meta.data$group4 <- "Normal"
fib_n7@meta.data$group4 <- "Normal"
fib_n8@meta.data$group4 <- "Normal"
fib_NL1@meta.data$group4 <- "AD_NL"
fib_NL2@meta.data$group4 <- "AD_NL"
fib_NL3@meta.data$group4 <- "AD_NL"
fib_NL4@meta.data$group4 <- "AD_NL"

fib_ad1@meta.data$group5 <- "AD15"
fib_ad2@meta.data$group5 <- "AD16"
fib_ad3@meta.data$group5 <- "AD17"
fib_ad4@meta.data$group5 <- "AD18"
fib_n1@meta.data$group5 <- "Normal18"
fib_n2@meta.data$group5 <- "Normal19"
fib_n3@meta.data$group5 <- "Normal20"
fib_n4@meta.data$group5 <- "Normal21"
fib_n5@meta.data$group5 <- "Normal22"
fib_n6@meta.data$group5 <- "Normal23"
fib_n6@meta.data$group5 <- "Normal24"
fib_n8@meta.data$group5 <- "Normal25"
fib_NL1@meta.data$group5 <- "AD_NL1"
fib_NL2@meta.data$group5 <- "AD_NL2"
fib_NL3@meta.data$group5 <- "AD_NL3"
fib_NL4@meta.data$group5 <- "AD_NL4"

fibad<-merge(fib_ad1,y=c(fib_ad2,fib_ad3,fib_ad4,fib_n1,fib_n2,fib_n3,fib_n4,fib_n5,fib_n6,fib_n7,fib_n8,fib_NL1,fib_NL2,fib_NL3,fib_NL4))
fibad_nFeature<-as.numeric(fibad@meta.data$nFeature_RNA)
class(fibad_nFeature)
fibad_nCount<-fibad@meta.data$nCount_RNA
class(fibad_nCount)
min(fibad_nFeature)
max(fibad_nFeature)
VlnPlot(object = fibad, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 1173.8-3188.2
fibad<-subset(fibad,subset = nFeature_RNA > 1146.3 & nFeature_RNA < 8748.7)
rm(fib_ad1,fib_ad2,fib_ad3,fib_ad4,fib_n1,fib_n2,fib_n3,fib_n4,fib_n5,fib_n6,fib_n7,fib_n8,fib_NL1,fib_NL2,fib_NL3,fib_NL4)
rm(fib_NL5)
rm(epidermis_fore12,epidermis_fore8,epidermis_fore9,epidermis_psoriasis1,epidermis_psoriasis2,epidermis_psoriasis3,epidermis_scalp1,
   epidermis_scalp2,epidermis_scalp3,epidermis_trunk1,epidermis_trunk2,epidermis_trunk3)

saveRDS(fibad, "GSE147424_AD.rds")

##sciencedata E-MTAB-8142
load("/data/syh/science/sciencedata.Rda")
seuratObject$group3 <- paste(seuratObject$sample_id,seuratObject$Status, sep = "_")
seuratObject$group<-seuratObject$Status
seuratObject$group2<-c("sciencedata")
seuratObject$group4<-paste(seuratObject$Status,seuratObject$Site, sep = "_")

seuratObject[["RNA"]]@meta.features<- data.frame(row.names = rownames(seuratObject[["RNA"]]))

seuratObject.list <- SplitObject(seuratObject, split.by = "sample_id")
rm(seuratObject)

#sample delete(<500)
seuratObject.list$SKN8090526<-NULL
seuratObject.list$SKN8090541<-NULL
seuratObject.list$SKN8090562<-NULL
seuratObject.list$SKN8090566<-NULL
seuratObject.list$SKN8090536<-NULL

SKN8090524<-seuratObject.list$SKN8090524
SKN8090524@meta.data$group5 <- "AD_NL5"

SKN8090525<-seuratObject.list$SKN8090525
table(SKN8090525$group4)
SKN8090525@meta.data$group5 <- "AD_NL6"

SKN8090527<-seuratObject.list$SKN8090527
table(SKN8090527$group4)
SKN8090527@meta.data$group5 <- "AD_NL7"

SKN8090528<-seuratObject.list$SKN8090528
table(SKN8090528$group4)
SKN8090528@meta.data$group5 <- "AD19"

SKN8090529<-seuratObject.list$SKN8090529
table(SKN8090529$group4)
SKN8090529@meta.data$group5 <- "AD20"

SKN8090530<-seuratObject.list$SKN8090530
table(SKN8090530$group4)
SKN8090530@meta.data$group5 <- "AD21"

SKN8090531<-seuratObject.list$SKN8090531
table(SKN8090531$group4)
SKN8090531@meta.data$group5 <- "AD22"

SKN8090537<-seuratObject.list$SKN8090537
table(SKN8090537$group4)
SKN8090537@meta.data$group5 <- "AD_NL8"

SKN8090538<-seuratObject.list$SKN8090538
table(SKN8090538$group4)
SKN8090538@meta.data$group5 <- "AD_NL9"

SKN8090539<-seuratObject.list$SKN8090539
table(SKN8090539$group4)
SKN8090539@meta.data$group5 <- "AD_NL10"

SKN8090540<-seuratObject.list$SKN8090540
table(SKN8090540$group4)
SKN8090540@meta.data$group5 <- "AD23"

SKN8090542<-seuratObject.list$SKN8090542
table(SKN8090542$group4)
SKN8090542@meta.data$group5 <- "AD24"

SKN8090543<-seuratObject.list$SKN8090543
table(SKN8090543$group4)
SKN8090543@meta.data$group5 <- "AD25"

SKN8090548<-seuratObject.list$SKN8090548
table(SKN8090548$group4)
SKN8090548@meta.data$group5 <- "AD_NL11"

SKN8090549<-seuratObject.list$SKN8090549
table(SKN8090549$group4)
SKN8090549@meta.data$group5 <- "AD_NL12"

SKN8090550<-seuratObject.list$SKN8090550
table(SKN8090550$group4)
SKN8090550@meta.data$group5 <- "AD_NL13"

SKN8090551<-seuratObject.list$SKN8090551
table(SKN8090551$group4)
SKN8090551@meta.data$group5 <- "AD_NL14"

SKN8090552<-seuratObject.list$SKN8090552
table(SKN8090552$group4)
SKN8090552@meta.data$group5 <- "AD26"

SKN8090553<-seuratObject.list$SKN8090553
table(SKN8090553$group4)
SKN8090553@meta.data$group5 <- "AD27"

SKN8090554<-seuratObject.list$SKN8090554
table(SKN8090554$group4)
SKN8090554@meta.data$group5 <- "AD28"

SKN8090555<-seuratObject.list$SKN8090555
table(SKN8090555$group4)
SKN8090555@meta.data$group5 <- "AD29"

SKN8090560<-seuratObject.list$SKN8090560
table(SKN8090560$group4)
SKN8090560@meta.data$group5 <- "AD30"

SKN8090561<-seuratObject.list$SKN8090561
table(SKN8090561$group4)
SKN8090561@meta.data$group5 <- "AD31"

SKN8090563<-seuratObject.list$SKN8090563
table(SKN8090563$group4)
SKN8090563@meta.data$group5 <- "AD_NL15"

SKN8090564<-seuratObject.list$SKN8090564
table(SKN8090564$group4)
SKN8090564@meta.data$group5 <- "AD32"

SKN8090565<-seuratObject.list$SKN8090565
table(SKN8090565$group4)
SKN8090565@meta.data$group5 <- "AD33"

SKN8090567<-seuratObject.list$SKN8090567
table(SKN8090567$group4)
SKN8090567@meta.data$group5 <- "AD_NL16"

SKN8090576<-seuratObject.list$SKN8090576
table(SKN8090576$group4)
SKN8090576@meta.data$group5 <- "Pso9"

SKN8090577<-seuratObject.list$SKN8090577
table(SKN8090577$group4)
SKN8090577@meta.data$group5 <- "Pso10"

SKN8090578<-seuratObject.list$SKN8090578
table(SKN8090578$group4)
SKN8090578@meta.data$group5 <- "Pso11"

SKN8090579<-seuratObject.list$SKN8090579
table(SKN8090579$group4)
SKN8090579@meta.data$group5 <- "Pso12"

SKN8090580<-seuratObject.list$SKN8090580
table(SKN8090580$group4)
SKN8090580@meta.data$group5 <- "Pso_NL1"

SKN8090581<-seuratObject.list$SKN8090581
table(SKN8090581$group4)
SKN8090581@meta.data$group5 <- "Pso_NL2"

SKN8090582<-seuratObject.list$SKN8090582
table(SKN8090582$group4)
SKN8090582@meta.data$group5 <- "Pso_NL3"

SKN8090583<-seuratObject.list$SKN8090583
table(SKN8090583$group4)
SKN8090583@meta.data$group5 <- "Pso_NL4"

SKN8090588<-seuratObject.list$SKN8090588
table(SKN8090588$group4)
SKN8090588@meta.data$group5 <- "Pso13"

SKN8090589<-seuratObject.list$SKN8090589
table(SKN8090589$group4)
SKN8090589@meta.data$group5 <- "Pso14"

SKN8090590<-seuratObject.list$SKN8090590
table(SKN8090590$group4)
SKN8090590@meta.data$group5 <- "Pso15"

SKN8090591<-seuratObject.list$SKN8090591
table(SKN8090591$group4)
SKN8090591@meta.data$group5 <- "Pso16"

SKN8090592<-seuratObject.list$SKN8090592
table(SKN8090592$group4)
SKN8090592@meta.data$group5 <- "Pso_NL5"

SKN8090593<-seuratObject.list$SKN8090593
table(SKN8090593$group4)
SKN8090593@meta.data$group5 <- "Pso_NL6"

SKN8090594<-seuratObject.list$SKN8090594
table(SKN8090594$group4)
SKN8090594@meta.data$group5 <- "Pso_NL7"

SKN8090595<-seuratObject.list$SKN8090595
table(SKN8090595$group4)
SKN8090595@meta.data$group5 <- "Pso_NL8"

SKN8090600<-seuratObject.list$SKN8090600
table(SKN8090600$group4)
SKN8090600@meta.data$group5 <- "Pso_NL9"

SKN8090601<-seuratObject.list$SKN8090601
table(SKN8090601$group4)
SKN8090601@meta.data$group5 <- "Pso_NL10"

SKN8090602<-seuratObject.list$SKN8090602
table(SKN8090602$group4)
SKN8090602@meta.data$group5 <- "Pso_NL11"

SKN8090603<-seuratObject.list$SKN8090603
table(SKN8090603$group4)
SKN8090603@meta.data$group5 <- "Pso_NL12"

SKN8090604<-seuratObject.list$SKN8090604
table(SKN8090604$group4)
SKN8090604@meta.data$group5 <- "Pso17"

SKN8090605<-seuratObject.list$SKN8090605
table(SKN8090605$group4)
SKN8090605@meta.data$group5 <- "Pso18"

SKN8090606<-seuratObject.list$SKN8090606
table(SKN8090606$group4)
SKN8090606@meta.data$group5 <- "Pso19"

SKN8090607<-seuratObject.list$SKN8090607
table(SKN8090607$group4)
SKN8090607@meta.data$group5 <- "Pso20"

SKN8104894<-seuratObject.list$SKN8104894
table(SKN8104894$group4)
SKN8104894@meta.data$group5 <- "Normal26"

SKN8104895<-seuratObject.list$SKN8104895
table(SKN8104895$group4)
SKN8104895@meta.data$group5 <- "Normal27"

SKN8104896<-seuratObject.list$SKN8104896
table(SKN8104896$group4)
SKN8104896@meta.data$group5 <- "Normal28"

SKN8104897<-seuratObject.list$SKN8104897
table(SKN8104897$group4)
SKN8104897@meta.data$group5 <- "Normal29"

SKN8104899<-seuratObject.list$SKN8104899
table(SKN8104899$group4)
SKN8104899@meta.data$group5 <- "Normal30"

SKN8104900<-seuratObject.list$SKN8104900
table(SKN8104900$group4)
SKN8104900@meta.data$group5 <- "Normal31"

SKN8104901<-seuratObject.list$SKN8104901
table(SKN8104901$group4)
SKN8104901@meta.data$group5 <- "Normal32"

SKN8104902<-seuratObject.list$SKN8104902
table(SKN8104902$group4)
SKN8104902@meta.data$group5 <- "Normal33"

SKN8105192<-seuratObject.list$SKN8105192
table(SKN8105192$group4)
SKN8105192@meta.data$group5 <- "Normal34"

SKN8105193<-seuratObject.list$SKN8105193
table(SKN8105193$group4)
SKN8105193@meta.data$group5 <- "Normal35"

SKN8105194<-seuratObject.list$SKN8105194
table(SKN8105194$group4)
SKN8105194@meta.data$group5 <- "Normal36"

SKN8105195<-seuratObject.list$SKN8105195
table(SKN8105195$group4)
SKN8105195@meta.data$group5 <- "Normal37"

SKN8105197<-seuratObject.list$SKN8105197
table(SKN8105197$group4)
SKN8105197@meta.data$group5 <- "Normal38"

SKN8105198<-seuratObject.list$SKN8105198
table(SKN8105198$group4)
SKN8105198@meta.data$group5 <- "Normal39"

SKN8105199<-seuratObject.list$SKN8105199
table(SKN8105199$group4)
SKN8105199@meta.data$group5 <- "Normal40"

SKN8105200<-seuratObject.list$SKN8105200
table(SKN8105200$group4)
SKN8105200@meta.data$group5 <- "Normal41"

`4820STDY7388991`<-seuratObject.list$`4820STDY7388991`
table(`4820STDY7388991`$group4)
`4820STDY7388991`@meta.data$group5 <- "Normal42"

`4820STDY7388992`<-seuratObject.list$`4820STDY7388992`
table(`4820STDY7388992`$group4)
`4820STDY7388992`@meta.data$group5 <- "Normal43"

`4820STDY7388993`<-seuratObject.list$`4820STDY7388993`
table(`4820STDY7388993`$group4)
`4820STDY7388993`@meta.data$group5 <- "Normal44"

`4820STDY7388994`<-seuratObject.list$`4820STDY7388994`
table(`4820STDY7388994`$group4)
`4820STDY7388994`@meta.data$group5 <- "Normal45"

`4820STDY7388995`<-seuratObject.list$`4820STDY7388995`
table(`4820STDY7388995`$group4)
`4820STDY7388995`@meta.data$group5 <- "Normal46"

`4820STDY7388996`<-seuratObject.list$`4820STDY7388996`
table(`4820STDY7388996`$group4)
`4820STDY7388996`@meta.data$group5 <- "Normal47"

`4820STDY7388997`<-seuratObject.list$`4820STDY7388997`
table(`4820STDY7388997`$group4)
`4820STDY7388997`@meta.data$group5 <- "Normal48"

`4820STDY7388998`<-seuratObject.list$`4820STDY7388998`
table(`4820STDY7388998`$group4)
`4820STDY7388998`@meta.data$group5 <- "Normal49"

`4820STDY7388999`<-seuratObject.list$`4820STDY7388999`
table(`4820STDY7388999`$group4)
`4820STDY7388999`@meta.data$group5 <- "Normal50"

`4820STDY7389000`<-seuratObject.list$`4820STDY7389000`
table(`4820STDY7389000`$group4)
`4820STDY7389000`@meta.data$group5 <- "Normal51"

`4820STDY7389001`<-seuratObject.list$`4820STDY7389001`
table(`4820STDY7389001`$group4)
`4820STDY7389001`@meta.data$group5 <- "Normal52"

`4820STDY7389002`<-seuratObject.list$`4820STDY7389002`
table(`4820STDY7389002`$group4)
`4820STDY7389002`@meta.data$group5 <- "Normal53"

`4820STDY7389003`<-seuratObject.list$`4820STDY7389003`
table(`4820STDY7389003`$group4)
`4820STDY7389003`@meta.data$group5 <- "Normal54"

`4820STDY7389004`<-seuratObject.list$`4820STDY7389004`
table(`4820STDY7389004`$group4)
`4820STDY7389004`@meta.data$group5 <- "Normal55"

`4820STDY7389005`<-seuratObject.list$`4820STDY7389005`
table(`4820STDY7389005`$group4)
`4820STDY7389005`@meta.data$group5 <- "Normal56"

`4820STDY7389006`<-seuratObject.list$`4820STDY7389006`
table(`4820STDY7389006`$group4)
`4820STDY7389006`@meta.data$group5 <- "Normal57"

`4820STDY7389007`<-seuratObject.list$`4820STDY7389007`
table(`4820STDY7389007`$group4)
`4820STDY7389007`@meta.data$group5 <- "Normal58"

`4820STDY7389008`<-seuratObject.list$`4820STDY7389008`
table(`4820STDY7389008`$group4)
`4820STDY7389008`@meta.data$group5 <- "Normal59"

`4820STDY7389009`<-seuratObject.list$`4820STDY7389009`
table(`4820STDY7389009`$group4)
`4820STDY7389009`@meta.data$group5 <- "Normal60"

`4820STDY7389010`<-seuratObject.list$`4820STDY7389010`
table(`4820STDY7389010`$group4)
`4820STDY7389010`@meta.data$group5 <- "Normal61"

`4820STDY7389011`<-seuratObject.list$`4820STDY7389011`
table(`4820STDY7389011`$group4)
`4820STDY7389011`@meta.data$group5 <- "Normal62"

`4820STDY7389012`<-seuratObject.list$`4820STDY7389012`
table(`4820STDY7389012`$group4)
`4820STDY7389012`@meta.data$group5 <- "Normal63"

`4820STDY7389013`<-seuratObject.list$`4820STDY7389013`
table(`4820STDY7389013`$group4)
`4820STDY7389013`@meta.data$group5 <- "Normal64"

`4820STDY7389014`<-seuratObject.list$`4820STDY7389014`
table(`4820STDY7389014`$group4)
`4820STDY7389014`@meta.data$group5 <- "Normal65"

seurat<-merge(SKN8090524,y=c(SKN8090525,SKN8090527,SKN8090528,SKN8090529,SKN8090530,
                             SKN8090531,SKN8090537,SKN8090538,SKN8090539,SKN8090540,SKN8090542,
                             SKN8090543,SKN8090548,SKN8090549,SKN8090550,SKN8090551,SKN8090552,
                             SKN8090553,SKN8090554,SKN8090555,SKN8090560,SKN8090561,SKN8090563,
                             SKN8090564,SKN8090565,SKN8090567,SKN8090576,SKN8090577,SKN8090578,
                             SKN8090579,SKN8090580,SKN8090581,SKN8090582,SKN8090583,SKN8090588,
                             SKN8090589,SKN8090590,SKN8090591,SKN8090592,SKN8090593,SKN8090594,
                             SKN8090595,SKN8090600,SKN8090601,SKN8090602,SKN8090603,SKN8090604,
                             SKN8090605,SKN8090606,SKN8090607,SKN8104894,SKN8104895,SKN8104896,
                             SKN8104897,SKN8104899,SKN8104900,SKN8104901,SKN8104902,SKN8105192,
                             SKN8105193,SKN8105194,SKN8105195,SKN8105197,SKN8105198,SKN8105199,
                             SKN8105200,`4820STDY7388991`,`4820STDY7388992`,`4820STDY7388993`,
                             `4820STDY7388994`,`4820STDY7388995`,`4820STDY7388996`,`4820STDY7388997`,
                             `4820STDY7388998`,`4820STDY7388999`,`4820STDY7389000`,`4820STDY7389001`,
                             `4820STDY7389002`,`4820STDY7389003`,`4820STDY7389004`,`4820STDY7389005`,
                             `4820STDY7389006`,`4820STDY7389007`,`4820STDY7389008`,`4820STDY7389009`,
                             `4820STDY7389010`,`4820STDY7389011`,`4820STDY7389012`,`4820STDY7389013`,`4820STDY7389014`))

rm(SKN8090524,SKN8090525,SKN8090527,SKN8090528,SKN8090529,SKN8090530,
   SKN8090531,SKN8090537,SKN8090538,SKN8090539,SKN8090540,SKN8090542,
   SKN8090543,SKN8090548,SKN8090549,SKN8090550,SKN8090551,SKN8090552,
   SKN8090553,SKN8090554,SKN8090555,SKN8090560,SKN8090561,SKN8090563,
   SKN8090564,SKN8090565,SKN8090567,SKN8090576,SKN8090577,SKN8090578,
   SKN8090579,SKN8090580,SKN8090581,SKN8090582,SKN8090583,SKN8090588,
   SKN8090589,SKN8090590,SKN8090591,SKN8090592,SKN8090593,SKN8090594,
   SKN8090595,SKN8090600,SKN8090601,SKN8090602,SKN8090603,SKN8090604,
   SKN8090605,SKN8090606,SKN8090607,SKN8104894,SKN8104895,SKN8104896,
   SKN8104897,SKN8104899,SKN8104900,SKN8104901,SKN8104902,SKN8105192,
   SKN8105193,SKN8105194,SKN8105195,SKN8105197,SKN8105198,SKN8105199,
   SKN8105200,`4820STDY7388991`,`4820STDY7388992`,`4820STDY7388993`,
   `4820STDY7388994`,`4820STDY7388995`,`4820STDY7388996`,`4820STDY7388997`,
   `4820STDY7388998`,`4820STDY7388999`,`4820STDY7389000`,`4820STDY7389001`,
   `4820STDY7389002`,`4820STDY7389003`,`4820STDY7389004`,`4820STDY7389005`,
   `4820STDY7389006`,`4820STDY7389007`,`4820STDY7389008`,`4820STDY7389009`,
   `4820STDY7389010`,`4820STDY7389011`,`4820STDY7389012`,`4820STDY7389013`,`4820STDY7389014`)
rm(seuratObject.list)
table(seurat$sample_id)
counts<-seurat@assays$RNA@counts
sciencedata<-CreateSeuratObject(counts =counts,project = "science",min.cells = 3, min.features = 200)

sciencedata$group<-seurat$group
sciencedata$group2<-seurat$group2
sciencedata$group3 <- seurat$group3
sciencedata$group4<-seurat$group4
sciencedata$group5<-seurat$group5
sciencedata[["RNA"]]@meta.features<- data.frame(row.names = rownames(seurat[["RNA"]]))
rm(seurat)
##sciencedata
science_nFeature<-as.numeric(sciencedata@meta.data$nFeature_RNA)
class(science_nFeature)
science_nCount<-sciencedata@meta.data$nCount_RNA
class(science_nCount)
min(science_nFeature)
max(science_nFeature)
VlnPlot(object = sciencedata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)#绘制小提琴图
## 10%-90% 812.2-5421.8
sciencedata<-subset(sciencedata,subset = nFeature_RNA > 812.2 & nFeature_RNA < 5421.8)
table(sciencedata$group4)

saveRDS(sciencedata, "sciencedata_new(group5).rds")

##doublefinder 
library(DoubletFinder)
table(sciencedata$group3)
sciencedata.list <- SplitObject(sciencedata, split.by = "group3")
sciencedata.list <- lapply(X = sciencedata.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(sciencedata.list)){sciencedata.list[[i]]@meta.data$DF_classification = sciencedata.list[[i]]@meta.data[,12]}

SKN8090524_Eczema<-sciencedata.list$SKN8090524_Eczema
sciencedata.list$SKN8090524_Eczema<-NULL

sciencedata.list_seuratobject<-merge(SKN8090524_Eczema,sciencedata.list)
table(sciencedata.list_seuratobject$DF_classification)

sciencedata.list_seuratobject<-subset(sciencedata.list_seuratobject,DF_classification=="Singlet")
rm(sciencedata.list)
rm(sciencedata)


##GSE222840 PN 
GSE222840 <- readRDS("/data/syh/AA_new/GSE222840.rds")
table(GSE222840$group3)
GSE222840.list <- SplitObject(GSE222840, split.by = "group3")
GSE222840.list <- lapply(X = GSE222840.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE222840.list)){GSE222840.list[[i]]@meta.data$DF_classification = GSE222840.list[[i]]@meta.data[,12]}

GSE222840_PN1<-GSE222840.list$GSE222840_PN1
GSE222840.list$GSE222840_PN1<-NULL

GSE222840.list_seuratobject<-merge(GSE222840_PN1,GSE222840.list)
table(GSE222840.list_seuratobject$DF_classification)
GSE222840.list_seuratobject<-subset(GSE222840.list_seuratobject,DF_classification=="Singlet")
rm(GSE222840,GSE222840.list)
rm(GSE222840_PN1)

##AA GSE212450
GSE212450 <- readRDS("/data/syh/AA_new/GSE212450_AA.rds")

GSE212450.list <- SplitObject(GSE212450, split.by = "group3")
GSE212450.list <- lapply(X = GSE212450.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE212450.list)){GSE212450.list[[i]]@meta.data$DF_classification = GSE212450.list[[i]]@meta.data[,12]}

GSE212450_AA2<-GSE212450.list$GSE212450_AA2
GSE212450.list$GSE212450_AA2<-NULL

GSE212450.list_seuratobject<-merge(GSE212450_AA2,GSE212450.list)
table(GSE212450.list_seuratobject$DF_classification)
GSE212450.list_seuratobject<-subset(GSE212450.list_seuratobject,DF_classification=="Singlet")
rm(GSE212450,GSE212450.list)
rm(GSE212450_AA2)

##fibad GSE147424
GSE147424 <- readRDS("/data/syh/AA_new/GSE147424_AD.rds")

GSE147424.list <- SplitObject(GSE147424, split.by = "group3")
GSE147424.list <- lapply(X = GSE147424.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE147424.list)){GSE147424.list[[i]]@meta.data$DF_classification = GSE147424.list[[i]]@meta.data[,12]}

GSE147424_ad1<-GSE147424.list$GSE147424_ad1
GSE147424.list$GSE147424_ad1<-NULL

GSE147424.list_seuratobject<-merge(GSE147424_ad1,GSE147424.list)
table(GSE147424.list_seuratobject$DF_classification)
GSE147424.list_seuratobject<-subset(GSE147424.list_seuratobject,DF_classification=="Singlet")
table(GSE147424.list_seuratobject$DF_classification)

rm(GSE147424.list,GSE147424)
rm(GSE147424_ad1)

##GSE150672 PSO
GSE150672 <- readRDS("/data/syh/AA_new/GSE150672_PSO.rds")

GSE150672.list <- SplitObject(GSE150672, split.by = "group3")
GSE150672.list <- lapply(X = GSE150672.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE150672.list)){GSE150672.list[[i]]@meta.data$DF_classification = GSE150672.list[[i]]@meta.data[,12]}

GSE150672_Normal1<-GSE150672.list$GSE150672_Normal1
GSE150672.list$GSE150672_Normal1<-NULL

GSE150672.list_seuratobject<-merge(GSE150672_Normal1,GSE150672.list)
table(GSE150672.list_seuratobject$DF_classification)
GSE150672.list_seuratobject<-subset(GSE150672.list_seuratobject,DF_classification=="Singlet")
rm(GSE150672_Normal1,GSE150672.list,GSE150672)

##GSE153760
GSE153760 <- readRDS("/data/syh/AA_new/GSE153760_AD.rds")

GSE153760.list <- SplitObject(GSE153760, split.by = "group3")
GSE153760.list <- lapply(X = GSE153760.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,npcs = 40,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE153760.list)){GSE153760.list[[i]]@meta.data$DF_classification = GSE153760.list[[i]]@meta.data[,12]}

GSE153760_HC1<-GSE153760.list$GSE153760_HC1
GSE153760.list$GSE153760_HC1<-NULL

GSE153760.list_seuratobject<-merge(GSE153760_HC1,GSE153760.list)
table(GSE153760.list_seuratobject$DF_classification)
GSE153760.list_seuratobject<-subset(GSE153760.list_seuratobject,DF_classification=="Singlet")
rm(GSE153760.list,GSE153760_HC1)
rm(GSE153760)

##GSE162183
GSE162183 <- readRDS("/data/syh/AA_new/GSE162183_PSO.rds")

GSE162183.list <- SplitObject(GSE162183, split.by = "group3")
GSE162183.list <- lapply(X = GSE162183.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,npcs = 40,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE162183.list)){GSE162183.list[[i]]@meta.data$DF_classification = GSE162183.list[[i]]@meta.data[,12]}

GSE162183_Ctrl1<-GSE162183.list$GSE162183_Ctrl1
GSE162183.list$GSE162183_Ctrl1<-NULL

GSE162183.list_seuratobject<-merge(GSE162183_Ctrl1,GSE162183.list)
table(GSE162183.list_seuratobject$DF_classification)
GSE162183.list_seuratobject<-subset(GSE162183.list_seuratobject,DF_classification=="Singlet")
rm(GSE162183.list,GSE162183_Ctrl1)
rm(GSE162183)

##HS_1 GSE175990
GSE175990 <- readRDS("/data/syh/AA_new/GSE175990_HS.rds")

GSE175990.list <- SplitObject(GSE175990, split.by = "group3")
GSE175990.list <- lapply(X = GSE175990.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE175990.list)){GSE175990.list[[i]]@meta.data$DF_classification = GSE175990.list[[i]]@meta.data[,12]}

`project_GSE175990-HS_1_1`<-GSE175990.list$`project_GSE175990-HS_1_1`
GSE175990.list$`project_GSE175990-HS_1_1`<-NULL

GSE175990.list_seuratobject<-merge(`project_GSE175990-HS_1_1`,GSE175990.list)
table(GSE175990.list_seuratobject$DF_classification)
GSE175990.list_seuratobject<-subset(GSE175990.list_seuratobject,DF_classification=="Singlet")
rm(`project_GSE175990-HS_1_1`,GSE175990.list)
rm(GSE175990)


##HS_2 GSE154775
GSE154775 <- readRDS("/data/syh/AA_new/GSE154775_HS.rds")

GSE154775.list <- SplitObject(GSE154775, split.by = "group3")
GSE154775.list <- lapply(X = GSE154775.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 2000)
  x<-ScaleData(x) 
  x<-RunPCA(x,verbose=FALSE)
  x <- RunUMAP(x, reduction = "pca", dims = 1:20)
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <-FindClusters(x,resolution = 0.4)
  sweep.res.list<- paramSweep_v3(x, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- x@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <-doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) 
})

for(i in 1:length(GSE154775.list)){GSE154775.list[[i]]@meta.data$DF_classification = GSE154775.list[[i]]@meta.data[,12]}

GSE154775_HS_2_1<-GSE154775.list$GSE154775_HS_2_1
GSE154775.list$GSE154775_HS_2_1<-NULL

GSE154775.list_seuratobject<-merge(GSE154775_HS_2_1,GSE154775.list)
table(GSE154775.list_seuratobject$DF_classification)
GSE154775.list_seuratobject<-subset(GSE154775.list_seuratobject,DF_classification=="Singlet")
rm(GSE154775_HS_2_1,GSE154775.list)
rm(GSE154775)

saveRDS(sciencedata.list_seuratobject, "sciencedata.list_seuratob ject.single.rds")
saveRDS(GSE222840.list_seuratobject, "GSE222840.list_seuratobject.single.rds")

saveRDS(GSE212450.list_seuratobject, "GSE212450.list_seuratobject.single.rds")
saveRDS(GSE147424.list_seuratobject, "GSE147424.list_seuratobject.single.rds")

saveRDS(GSE150672.list_seuratobject, "GSE150672.list_seuratobject.single.rds")

saveRDS(GSE153760.list_seuratobject, "GSE153760.list_seuratobject.single.rds")

saveRDS(GSE162183.list_seuratobject, "GSE162183.list_seuratobject.single.rds")
saveRDS(GSE175990.list_seuratobject, "GSE175990.list_seuratobject.single.rds")

saveRDS(GSE154775.list_seuratobject, "GSE154775.list_seuratobject.single.rds")


##combine the data
Infla_scRNA_harmony<-merge(sciencedata.list_seuratobject.single,
                           y=c(GSE222840.list_seuratobject,
                               GSE212450.list_seuratobject,
                               GSE147424.list_seuratobject,
                               GSE150672.list_seuratobject,
                               GSE153760.list_seuratobject,
                               GSE162183.list_seuratobject,
                               GSE175990.list_seuratobject,
                               GSE154775.list_seuratobject
                           ))
rm(sciencedata.list_seuratobject.single,GSE222840.list_seuratobject,
   GSE212450.list_seuratobject,
   GSE147424.list_seuratobject,
   GSE150672.list_seuratobject,
   GSE153760.list_seuratobject,
   GSE162183.list_seuratobject,
   GSE175990.list_seuratobject,
   GSE154775.list_seuratobject)

table(Infla_scRNA_harmony$group)
table(Infla_scRNA_harmony$group2)
table(Infla_scRNA_harmony$group3)
table(Infla_scRNA_harmony$group4)
table(Infla_scRNA_harmony$group5)

group<-Infla_scRNA_harmony@meta.data$group
group4<-Infla_scRNA_harmony@meta.data$group4

group<-gsub("NL_ad", "AD", group)
group<-gsub("Eczema", "AD", group)
group<-gsub("Healthy", "Normal", group)
group<-gsub("psoriasis", "Psoriasis", group)
group<-gsub("normal", "Normal", group)
group<-gsub("ad", "AD", group)
table(group)
Infla_scRNA_harmony@meta.data$"type"=group

table(group4)
group4<-gsub("Eczema_lesion", "AD", group4)
group4<-gsub("Eczema_non_lesion", "AD_NL", group4)
group4<-gsub("Healthy_non_lesion", "Normal", group4)
group4<-gsub("Psoriasis_lesion", "Psoriasis", group4)
group4<-gsub("Psoriasis_non_lesion", "Psoriasis_NL", group4)
group4<-gsub("Psoriais", "Psoriasis", group4)
table(group4)
Infla_scRNA_harmony@meta.data$"type_NL"=group4

##filter
##percent.mito
Infla_scRNA_harmony[["percent.mito"]] <- PercentageFeatureSet(Infla_scRNA_harmony, pattern = "^MT",assay = 'RNA')
mito.genes <- grep(pattern = "^MT", x = rownames(x = Infla_scRNA_harmony@assays$RNA), 
                   value = TRUE)
percent.mito <- Matrix::colSums( Infla_scRNA_harmony@assays$RNA[mito.genes, ])/
  Matrix::colSums(Infla_scRNA_harmony@assays$RNA)

Infla_scRNA_harmony<- AddMetaData(object = Infla_scRNA_harmony, metadata = percent.mito,
                                  col.name = "percent.mito")

##percent.ribo
rb.genes <- rownames(Infla_scRNA_harmony)[grep("^RP[SL]",rownames(Infla_scRNA_harmony))]
percent.ribo <- Matrix::colSums( Infla_scRNA_harmony@assays$RNA[rb.genes,])/Matrix::colSums(Infla_scRNA_harmony@assays$RNA)
Infla_scRNA_harmony <- AddMetaData(Infla_scRNA_harmony, percent.ribo, col.name = "percent.ribo")


VlnPlot(object = Infla_scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 5)#绘制小提琴图

Infla_scRNA_harmony <-subset(Infla_scRNA_harmony, subset =  percent.mito>0&percent.mito < 0.1&percent.ribo>0.05&nCount_RNA<50000)
VlnPlot(object = Infla_scRNA_harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo","percent.HB"), ncol = 5)#绘制小提琴图
VlnPlot(object = Infla_scRNA_harmony_res, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo","percent.HB"), ncol = 5)#绘制小提琴图

##Harmony
Infla_scRNA_harmony <- NormalizeData(Infla_scRNA_harmony) 
Infla_scRNA_harmony <-FindVariableFeatures(Infla_scRNA_harmony) 
Infla_scRNA_harmony<-ScaleData(Infla_scRNA_harmony ) 
Infla_scRNA_harmony<-RunPCA(Infla_scRNA_harmony,verbose=FALSE)
Infla_scRNA_harmony <- RunHarmony(Infla_scRNA_harmony, group.by.vars = "group3",plot_convergence = TRUE)
saveRDS(Infla_scRNA_harmony, "Infla_scRNA_harmony.rds")
rm(Infla_scRNA_harmony)
##dim 1:20
Infla_scRNA_harmony_res0.3 <- RunUMAP(Infla_scRNA_harmony, reduction = "harmony", dims = 1:20)
Infla_scRNA_harmony_res0.3 <- FindNeighbors(Infla_scRNA_harmony_res0.3, reduction = "harmony", dims = 1:20)
Infla_scRNA_harmony_res0.3 <- FindClusters(Infla_scRNA_harmony_res0.3, resolution = 0.3)
table(Infla_scRNA_harmony_res0.3$RNA_snn_res.0.3)  ##25群
Infla_scRNA_harmony_res0.3 <- FindClusters(Infla_scRNA_harmony_res0.3, resolution = 0.1)
table(Infla_scRNA_harmony_res$RNA_snn_res.0.1)  ##16群 选用的是这种

markers_res0.1=FindAllMarkers(Infla_scRNA_harmony_res0.3,only.pos = TRUE)
write.xlsx(markers_res0.1, file = "markers_res0.1.xlsx",rowNames=TRUE)

marker1<-c("KRT5","KRT15","KRT14","KRT1","KRT10","KRT6A","KRT16","UBE2C","TOP2A","GJB2","DCN","COL1A1","LUM","PECAM1",'VWF',"ACTA2","TAGLN","PMEL","MLANA","CD3D",'CD3G','CD3E','CD8A','LYZ','HLA-DRA','CD68','ITGAX',"IGHG1","IGKC","TPSB2","TPSAB1")
DotPlot(Infla_scRNA_harmony_res0.3,features = marker1,cols = c('#f6eff7','#1c9099'),assay="RNA")+RotatedAxis()
DimPlot(Infla_scRNA_harmony_res0.3,reduction = "umap", label=TRUE,cols = palette1)

Infla_scRNA_harmony_res0.3 <- RunUMAP(Infla_scRNA_harmony_res0.3, reduction = "harmony", dims = 1:10)
Infla_scRNA_harmony_res0.3 <- FindNeighbors(Infla_scRNA_harmony_res0.3, reduction = "harmony", dims = 1:10)
Infla_scRNA_harmony_res0.3 <- FindClusters(Infla_scRNA_harmony_res0.3, resolution = 0.2)
Infla_scRNA_harmony_res0.3<- RenameIdents(Infla_scRNA_harmony_res0.3, "0"="KC.1","1"="T_C","2"="FIB","3"="MAC_DC","4"="EC","5"="KC.2","6"="PC_vSMC","7"="KC.3","8"="EC","9"="MEL","10"="B/Mast_C",
                                          "11"='KC.4',"12"='KC.4',"13"='KC.4',"15"='KC.4',"14"='MAC_DC',"16"="T_C")

levels(Infla_scRNA_harmony_res0.3)<-c("KC.1","KC.2","KC.3","KC.4",
                                      "FIB","EC","PC_vSMC","MEL",
                                      "T_C","MAC_DC","B/Mast_C")
Infla_scRNA_harmony_res0.3$cell.type<-Idents(Infla_scRNA_harmony_res0.3)

DimPlot(Infla_scRNA_harmony_res0.3,reduction = "umap", label=TRUE,cols = palette3)
DimPlot(Infla_scRNA_harmony_res0.3,reduction = "umap", label=FALSE,cols = palette3)
table(Infla_scRNA_harmony_res0.3$type)
Idents(Infla_scRNA_harmony_res0.3)<-Infla_scRNA_harmony_res0.3$type
levels(Infla_scRNA_harmony_res0.3)<-c("Normal","AD","Psoriasis","PN","AA","HS")
Idents(Infla_scRNA_harmony_res0.3)<-Infla_scRNA_harmony_res0.3$cell.type

Infla_scRNA_harmony_res0.3@meta.data$type <- factor(Infla_scRNA_harmony_res0.3@meta.data$type, 
                                                    levels=c("Normal","AD","Psoriasis","PN","AA","HS"))

DotPlot(Infla_scRNA_harmony_res0.3,features = marker1,cols = c('#f6eff7','#1c9099'),assay="RNA")+RotatedAxis()
DimPlot(Infla_scRNA_harmony_res0.3,reduction = "umap", label=TRUE,cols = palette3,split.by = "type")
P<-DoHeatmap(Infla_scRNA_harmony_res0.3, features = marker1,group.colors=palette3) + NoLegend()
saveRDS(Infla_scRNA_harmony_res0.3, "Infla_scRNA_harmony_res0.2.rds")

##inflammatory pathways
th1<-c("IFNG","CXCL9","CXCL10","STAT1","IL-12P35","IL-12RB1","IL-27P28","CXCL11","CCL2","CCL3","CCL4") 
names(th1)<-c("th1")
infla<-getGmt("/data/scRNA_mix/zrh_mix/gene_infla.gmt")
infla_gene<-read.gmt("/data/scRNA_mix/zrh_mix/gene_infla.gmt")
infla_gene2<-read.gmt("/data/scRNA_mix/zrh_mix/gene_infla.gmt")
infla_gene<-infla_gene[,2]
exp_infla<-subset(Infla_scRNA_harmony_res0.2,features=infla_gene)
exp_infla<-as.matrix(exp_infla@assays$RNA@data)
cells_rankings_infla<- AUCell_buildRankings(exp_infla, nCores=7, plotStats=TRUE) 
names(infla)
cells_AUC_infla <- AUCell_calcAUC(infla, cells_rankings_infla,nCores = 5,normAUC = TRUE,aucMaxRank = ceiling(1 * nrow(cells_rankings_infla)))
auc_infla<-getAUC(cells_AUC_infla)
auc_infla_assay<-CreateAssayObject(auc_infla)
Infla_scRNA_harmony_res0.2@assays$infla<-auc_infla_assay
Infla_scRNA_harmony_res0.2@assays$infla@key<-"infla_"
Infla_scRNA_harmony_res0.2<-ScaleData(Infla_scRNA_harmony_res0.2,assay ="infla")
options(scipen = 200)
table(Infla_scRNA_harmony_res$orig.ident,Infla_scRNA_harmony_res$type)

##metabolism
kegg<-getGmt("/data/scRNA_mix/zrh_mix/KEGG_metabolism.gmt")
kegg_gene<-read.gmt("/data/scRNA_mix/zrh_mix/KEGG_metabolism.gmt")
kegg_gene<-kegg_gene[,2]
exp<-subset(Infla_scRNA_harmony_res0.2,features=kegg_gene)
exp<-as.matrix(exp@assays$RNA@data)
cells_rankings_scale<- AUCell_buildRankings(exp, nCores=7, plotStats=TRUE) 
names(kegg)
cells_AUC <- AUCell_calcAUC(kegg, cells_rankings_scale,nCores = 5,normAUC = TRUE,aucMaxRank = ceiling(1 * nrow(cells_rankings_scale)))
auc<-getAUC(cells_AUC)
auc_assay<-CreateAssayObject(auc)
Infla_scRNA_harmony_res0.2@assays$metabolic<-auc_assay
Infla_scRNA_harmony_res0.2@assays$metabolic@key<-"metabolic_"
Infla_scRNA_harmony_res0.2<-ScaleData(Infla_scRNA_harmony_res0.2,assay ="metabolic")

Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$group5
auc_ave_Infla_group5<-AverageExpression(Infla_scRNA_harmony_res0.2,assay ="infla",slot="scale.data")
write.xlsx(auc_ave_Infla_group5, file = "auc_ave_Infla_group5.xlsx",rowNames=TRUE)

Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$type_NL
auc_ave_Infla_type<-AverageExpression(Infla_scRNA_harmony_res0.2,assay ="infla",slot="data")
Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$group5
auc_ave_Infla_group5<-AverageExpression(Infla_scRNA_harmony_res0.2,assay ="infla",slot="data")
rm(exp_infla)
write.xlsx(auc_ave_Infla_group5, file = "auc_ave_Infla_group5.data.average.xlsx",rowNames=TRUE)
rm(Infla_scRNA_harmony_res0.2)

##inflammatory phenotypes
auc_infla_classification <- read_excel("auc_infla_classification.xlsx")
auc_infla_classification<-auc_infla_classification%>% 
  dplyr::select(...1,classification)
names(auc_infla_classification)<-c("col","classification")

auc_infla<-auc_ave_Infla_group5$infla
auc_infla<-data.frame(auc_infla)
auc_infla<-auc_infla%>% 
  select(-Normal1,-Normal2,-Normal3,-Normal4,-Normal5,-Normal6,-Normal7,-Normal8,-Normal9,-Normal10,
         -Normal11,-Normal12,-Normal13,-Normal14,-Normal15,-Normal16,-Normal17,-Normal18,-Normal19,-Normal20,
         -Normal21,-Normal22,-Normal24,-Normal25,-Normal27,-Normal28,-Normal29,-Normal30,
         -Normal31,-Normal32,-Normal33,-Normal34,-Normal35,-Normal36,-Normal37,-Normal38,-Normal39,-Normal40,
         -Normal41,-Normal42,-Normal43,-Normal44,-Normal45,-Normal46,-Normal47,-Normal48,-Normal49,-Normal50,
         -Normal51,-Normal52,-Normal53,-Normal54,-Normal55,-Normal56,-Normal57,-Normal58,-Normal59,-Normal60,
         -Normal61,-Normal62,-Normal63,-Normal64,-Normal65)

##Inflammatory PCA
col<-colnames(auc_infla)
col<-data.frame(col)

col<-merge(col,auc_infla_classification,by="col",sort=F)
col<-col%>% 
  dplyr::filter(classification!=c("ns"))

auc_infla<-auc_infla%>% 
  dplyr::select(-c("AD32","Pso_NL12","AD7"))

annotation_col = data.frame(col)
rownames(annotation_col)=colnames((auc_infla))
colnames(annotation_col)=c("CellType","group")

annotation_col<-cbind(rownames(annotation_col),annotation_col)
group_infor<-annotation_col
names(group_infor)<-c("colnames","group1","group2")
class(group_infor)
options(ggrepel.max.overlaps = Inf)
data_pca=t(auc_infla)
data_pca=cbind(data_pca,group_infor)
res.pca <- PCA(data_pca[,1:3], graph = TRUE)
rownames(group_infor)=group_infor[,1]
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
res.pca.p<- fviz_pca_ind(res.pca, geom = "point", col.ind =data_pca$group2 )
ggpubr::ggpar(res.pca.p,
              title = "Principal Component Analysis",
              xlab = "PC1", ylab = "PC2", shape ="data_pca$group2",
              legend.title = "Species", legend.position = "top", addEllipses = TRUE,
              palette = "jco",xlim=c(-5,5),ylim=c(-3,3))+stat_ellipse( linetype = 2,level = 0.65,aes(group = data_pca$group2, colour = data_pca$group2))

fviz_pca_ind(res.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800", col.ind =data_pca$group2, addEllipses = TRUE,
             repel = TRUE,ylim=c(-10,10) # Avoid text overlapping (slow if many points)
)

var <- get_pca_var(res.pca)
coord=var$coord
corrplot(var$cos2, is.corr=FALSE)

##cell proportion
cell.prop<-as.data.frame(prop.table(table(Infla_scRNA_harmony_res0.2$cell.type,Infla_scRNA_harmony_res0.2@meta.data$type)))
names(cell.prop)<-c("cluster","patient","proportion")
ggplot(cell.prop,aes(patient,proportion,fill=cluster))+
  
  geom_bar(stat="identity",position="fill")+
  
  ggtitle("")+
  
  theme_bw()+
  
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = palette3)

##Add inflammatory phenotypes
group6<-Infla_scRNA_harmony_res0.2$group5
table(group6)
group6<-gsub("Pso_NL12", "type3", group6) ##gsub会把带有字符的都算进去，有问题，不完全匹配

for(i in length(col$col):1){
  A<-col$col[i]   
  B<-col$classification[i]
  C<- which(group6==A,arr.ind = T)
  group6[C]<-B
}
table(group6)
group6<-gsub("Pso_NL12", "ns", group6) 
group6<-gsub("AD7", "ns", group6) 
group6<-gsub("AD32", "ns", group6) 
Infla_scRNA_harmony_res0.2@meta.data$group6=group6
Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$cell.type
DimPlot(Infla_scRNA_harmony_res0.2,reduction = "umap", label=TRUE,cols = palette3,split.by = "group6")

group7<-Infla_scRNA_harmony_res0.2@meta.data$group6
table(group7)
group7<-gsub("Normal1", "Normal", group7) 
group7<-gsub("Normal2", "Normal", group7) 
group7<-gsub("Normal3", "Normal", group7) 
group7<-gsub("Normal4", "Normal", group7) 
group7<-gsub("Normal5", "Normal", group7) 
group7<-gsub("Normal6", "Normal", group7) 
group7<-gsub("Normal7", "Normal", group7) 
group7<-gsub("Normal8", "Normal", group7) 
group7<-gsub("Normal9", "Normal", group7) 
group7<-gsub("Normal10", "Normal", group7) 
group7<-gsub("Normal11", "Normal", group7) 
group7<-gsub("Normal12", "Normal", group7) 
group7<-gsub("Normal13", "Normal", group7) 
group7<-gsub("Normal14", "Normal", group7) 
group7<-gsub("Normal15", "Normal", group7) 
group7<-gsub("Normal16", "Normal", group7) 
group7<-gsub("Normal17", "Normal", group7) 
group7<-gsub("Normal18", "Normal", group7) 
group7<-gsub("Normal19", "Normal", group7) 
group7<-gsub("Normal20", "Normal", group7) 
group7<-gsub("Normal21", "Normal", group7) 
group7<-gsub("Normal22", "Normal", group7) 
group7<-gsub("Normal23", "Normal", group7) 
group7<-gsub("Normal24", "Normal", group7) 
group7<-gsub("Normal25", "Normal", group7) 
group7<-gsub("Normal26", "Normal", group7) 
group7<-gsub("Normal27", "Normal", group7) 
group7<-gsub("Normal28", "Normal", group7) 
group7<-gsub("Normal29", "Normal", group7) 
group7<-gsub("Normal30", "Normal", group7) 
group7<-gsub("Normal31", "Normal", group7) 
group7<-gsub("Normal32", "Normal", group7) 
group7<-gsub("Normal33", "Normal", group7) 
group7<-gsub("Normal34", "Normal", group7) 
group7<-gsub("Normal35", "Normal", group7) 
group7<-gsub("Normal36", "Normal", group7) 
group7<-gsub("Normal37", "Normal", group7) 
group7<-gsub("Normal38", "Normal", group7) 
group7<-gsub("Normal39", "Normal", group7) 
group7<-gsub("Normal40", "Normal", group7) 
group7<-gsub("Normal41", "Normal", group7) 
group7<-gsub("Normal42", "Normal", group7) 
group7<-gsub("Normal43", "Normal", group7) 
group7<-gsub("Normal44", "Normal", group7) 
group7<-gsub("Normal45", "Normal", group7) 
group7<-gsub("Normal46", "Normal", group7) 
group7<-gsub("Normal47", "Normal", group7) 
group7<-gsub("Normal48", "Normal", group7) 
group7<-gsub("Normal49", "Normal", group7) 
group7<-gsub("Normal50", "Normal", group7) 
group7<-gsub("Normal51", "Normal", group7) 
group7<-gsub("Normal52", "Normal", group7) 
group7<-gsub("Normal53", "Normal", group7) 
group7<-gsub("Normal54", "Normal", group7) 
group7<-gsub("Normal55", "Normal", group7) 
group7<-gsub("Normal56", "Normal", group7) 
group7<-gsub("Normal57", "Normal", group7) 
group7<-gsub("Normal58", "Normal", group7) 
group7<-gsub("Normal59", "Normal", group7) 
group7<-gsub("Normal60", "Normal", group7) 
group7<-gsub("Normal61", "Normal", group7) 
group7<-gsub("Normal62", "Normal", group7) 
group7<-gsub("Normal63", "Normal", group7) 
group7<-gsub("Normal64", "Normal", group7) 
group7<-gsub("Normal65", "Normal", group7) 
group7<-gsub("Normal0", "Normal", group7) 
group7<-gsub("Normal1", "Normal", group7) 
group7<-gsub("Normal2", "Normal", group7) 
group7<-gsub("Normal3", "Normal", group7) 
group7<-gsub("Normal4", "Normal", group7) 
group7<-gsub("Normal5", "Normal", group7) 
table(group7)

Infla_scRNA_harmony_res0.2@meta.data$group7=group7
Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$cell.type
DimPlot(Infla_scRNA_harmony_res0.2,reduction = "umap", label=TRUE,cols = palette3,split.by = "group7")

saveRDS(Infla_scRNA_harmony_res0.2, "Infla_scRNA_harmony_res0.2.rds")

Infla_scRNA_harmony_res0.2<-subset(Infla_scRNA_harmony_res0.2,
                                   group7==c("Normal","type1",
                                             "type1/type2","type1/type3","type2",
                                             "type2/type3","type3"))

DimPlot(Infla_scRNA_harmony_res0.2,reduction = "umap", label=TRUE,cols = palette3,split.by = "group7")

cell.prop<-as.data.frame(prop.table(table(Infla_scRNA_harmony_res0.2$cell.type,Infla_scRNA_harmony_res0.2@meta.data$group7)))
names(cell.prop)<-c("cluster","patient","proportion")
ggplot(cell.prop,aes(patient,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = palette3)

##Normal\NL\Infla tissues
group8<-Infla_scRNA_harmony_res0.2@meta.data$type_NL
table(group8)
group8<-gsub("AD_NL", "NL", group8) 
group8<-gsub("AA_NL", "NL", group8) 
group8<-gsub("Psoriasis_NL", "NL", group8) 
group8<-gsub("AA", "Infla", group8) 
group8<-gsub("AD", "Infla", group8) 
group8<-gsub("HS", "Infla", group8) 
group8<-gsub("PN", "Infla", group8) 
group8<-gsub("Psoriasis", "Infla", group8) 
Infla_scRNA_harmony_res0.2$group8<-group8
Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$group8
Idents(Infla_scRNA_harmony_res0.2)<-Infla_scRNA_harmony_res0.2$cell.type
auc_ave_Infla<-AverageExpression(Infla_scRNA_harmony_res0.2,assay ="infla",add.ident = "group8",slot="scale.data")
auc_ave_Infla<-auc_ave_Infla$infla
auc_ave_Infla<-data.frame(t(auc_ave_Infla))
auc_ave_Infla<-auc_ave_Infla[,-4]

pheatmap(auc_ave_Infla,cluster_rows = T,cluster_cols = F,color=
           paletteer_c("scico::broc",n=100))
pheatmap(auc_ave_Infla,cluster_rows = F,cluster_cols = F,color=
           paletteer_c("scico::broc",n=100))

auc_ave_Infla<-AverageExpression(Infla_scRNA_harmony_res0.2,assay ="infla",add.ident = "type",slot="scale.data")
auc_ave_Infla<-auc_ave_Infla$infla
auc_ave_Infla<-data.frame(t(auc_ave_Infla))
auc_ave_Infla<-auc_ave_Infla[,-4]
pheatmap(auc_ave_Infla,cluster_rows = T,cluster_cols = F,color=
           paletteer_c("scico::broc",n=100))