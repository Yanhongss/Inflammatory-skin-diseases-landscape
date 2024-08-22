##Mapping 
##KC
Infla_scRNA_harmony_res0.2<- readRDS("/data/syh/AA_new/Infla_scRNA_harmony_res0.2_metabolic.rds")
Idents(Infla_scRNA_harmony_res0.2)
KC<-subset(Infla_scRNA_harmony_res0.2,idents=c("KC.1","KC.2","KC.3","KC.4"))
meta.data<-KC@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

KC<-CreateSeuratObject(counts = KC@assays$RNA@counts, meta.data = meta.data)
KC <- NormalizeData(KC) 
KC <-FindVariableFeatures(KC) 
KC<-ScaleData(KC ) 
KC<-RunPCA(KC,verbose=FALSE)
KC <- RunHarmony(KC, group.by.vars = "group3",plot_convergence = TRUE)
ElbowPlot(KC)

##dim 1:15
KC <- RunUMAP(KC, reduction = "harmony", dims = 1:15)
KC <- FindNeighbors(KC, reduction = "harmony", dims = 1:15)
KC <- FindClusters(KC, resolution = 0.2)
DimPlot(KC,reduction = "umap", label=TRUE)
DotPlot(KC,features = marker_KC,cols = c('#f6eff7','#1c9099'),assay="RNA")+RotatedAxis()

KC_marker<-FindAllMarkers(KC,only.pos = TRUE)
write.xlsx(KC_marker, file = "KC_marker_dim15_res0.2.xlsx",rowNames=TRUE)

KC <- FindClusters(object = KC,resolution = c(seq(.1,1,.1)))
clustree(KC@meta.data, prefix = "RNA_snn_res.")

##update cellname
Idents(KC)<-KC$seurat_clusters
KC<- RenameIdents(KC, "0"="KC2","1"="KC1","2"="KC3","3"="KC5","4"="KC7","5"="KC6","6"="others","7"="others","8"="KC2","9"="KC4","10"="others", "11"='KC2',"12"='others',"13"='others')
levels(KC)<-c("KC1","KC2","KC3","KC4","KC5","KC6","KC7","others")

##T
Tcells<-subset(Infla_scRNA_harmony_res0.2,idents=c("T_C"))
meta.data<-Tcells@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

Tcells<-CreateSeuratObject(counts = Tcells@assays$RNA@counts, meta.data = meta.data)
Tcells <- NormalizeData(Tcells) 
Tcells <-FindVariableFeatures(Tcells) 
Tcells<-ScaleData(Tcells ) 
Tcells<-RunPCA(Tcells,verbose=FALSE)
Tcells <- RunHarmony(Tcells, group.by.vars = "group3",plot_convergence = TRUE)
ElbowPlot(Tcells)

##dim 1:15
Tcells <- RunUMAP(Tcells, reduction = "harmony", dims = 1:15)
Tcells <- FindNeighbors(Tcells, reduction = "harmony", dims = 1:15)
Tcells <- FindClusters(Tcells, resolution = 1.2)

Tcells <- FindClusters(object = Tcells, resolution = c(seq(.1,1,.1)))
clustree(Tcells@meta.data, prefix = "RNA_snn_res.")

Tcells <- FindClusters(Tcells, resolution = 1.2)
KC_marker<-FindAllMarkers(KC,only.pos = TRUE)
Tcells_marker<-FindAllMarkers(Tcells,only.pos = TRUE)
write.xlsx(Tcells_marker, file = "Tcells_marker_res1.2.xlsx",rowNames=TRUE)

##update cellname
Idents(Tcells)<-Tcells$seurat_clusters
Tcells<- RenameIdents(Tcells, "0"="Trm1","1"="Trm3","2"="Tstr","3"="Trm1","4"="CTLex","5"="Treg", "6"="Trm1","7"="Treg","8"="CTLax","9"="NKT","10"="Trm1", "11"='NKT',"12"='others',"13"='Trm5',"14"="Trm6","15"="Trm4", "16"="NKT","17"="NKT","18"="others","19"="others","20"="Trm2","21"="Tmm","22"="Tcm","23"="others")
levels(Tcells)<-c("CTLax","CTLex","Treg","NKT", "Trm1","Trm2","Trm3","Trm4","Trm5","Trm6","Tmm","Tcm","Tstr","others")

##FIB
FIB<-subset(Infla_scRNA_harmony_res0.2,idents=c("FIB"))
meta.data<-FIB@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

FIB<-CreateSeuratObject(counts = FIB@assays$RNA@counts, meta.data = meta.data)
FIB <- NormalizeData(FIB) 
FIB <-FindVariableFeatures(FIB) 
FIB<-ScaleData(FIB ) 
FIB<-RunPCA(FIB,verbose=FALSE)
FIB <- RunHarmony(FIB, group.by.vars = "group3",plot_convergence = TRUE)
ElbowPlot(FIB)

##dim 1:15
FIB <- RunUMAP(FIB, reduction = "harmony", dims = 1:15)
FIB <- FindNeighbors(FIB, reduction = "harmony", dims = 1:15)
library(clustree)
FIB <- FindClusters(object = FIB,resolution = c(seq(.1,1,.1)))
clustree(FIB@meta.data, prefix = "RNA_snn_res.")

FIB <- FindClusters(FIB, resolution = 0.7)
FIB_marker<-FindAllMarkers(FIB,only.pos = TRUE)
write.xlsx(FIB_marker, file = "FIB_marker_dim15_res0.7.xlsx",rowNames=TRUE)

##update cellname
Idents(FIB)
FIB<- RenameIdents(FIB, "0"="FIB1","1"="FIB3","2"="FIB3","3"="FIB2","4"="FIB1", "5"="FIB3", "6"="FIB4","7"="FIB5","8"="FIB2","9"="FIB3","10"="FIB4","11"='FIB6',"12"='FIB7',"13"='FIB8',"14"="FIB9","15"="FIB10","16"="others","17"="others")
levels(FIB)<-c("FIB1","FIB2","FIB3","FIB4","FIB5","FIB6", "FIB7",'FIB8',"FIB9","FIB10", "others")

##EC
EC<-subset(Infla_scRNA_harmony_res0.2,idents=c("EC"))
meta.data<-EC@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

EC<-CreateSeuratObject(counts = EC@assays$RNA@counts, meta.data = meta.data)
EC <- NormalizeData(EC) 
EC <-FindVariableFeatures(EC) 
EC<-ScaleData(EC ) 
EC<-RunPCA(EC,verbose=FALSE)
EC <- RunHarmony(EC, group.by.vars = "group3",plot_convergence = TRUE)
##dim 1:15
EC <- RunUMAP(EC, reduction = "harmony", dims = 1:15)
EC <- FindNeighbors(EC, reduction = "harmony", dims = 1:15)
library(clustree)
EC <- FindClusters( object = EC,resolution = c(seq(.1,1,.1)))
clustree(EC@meta.data, prefix = "RNA_snn_res.")
EC <- FindClusters(EC, resolution = 0.2)
EC_marker<-FindAllMarkers(EC,only.pos = TRUE)
write.xlsx(EC_marker, file = "EC_marker_dim15_res0.2.xlsx",rowNames=TRUE)

##update cellname
EC<- RenameIdents(EC, "0"="VE1","1"="VE2","2"="LE","3"="VE3","4"="VE4","5"="others", "6"="VE5","7"="others","8"="VE1","9"="LE")
levels(EC)<-c("VE1","VE2","VE3","VE4","VE5","LE","others")

##MAC_DC
MAC_DC<-subset(Infla_scRNA_harmony_res0.2,idents=c("MAC_DC"))
meta.data<-MAC_DC@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

MAC_DC<-CreateSeuratObject(counts = MAC_DC@assays$RNA@counts, meta.data = meta.data)
MAC_DC <- NormalizeData(MAC_DC) 
MAC_DC <-FindVariableFeatures(MAC_DC) 
MAC_DC<-ScaleData(MAC_DC ) 
MAC_DC<-RunPCA(MAC_DC,verbose=FALSE)
MAC_DC <- RunHarmony(MAC_DC, group.by.vars = "group3",plot_convergence = TRUE)
##dim 1:15
MAC_DC <- RunUMAP(MAC_DC, reduction = "harmony", dims = 1:15)
MAC_DC <- FindNeighbors(MAC_DC, reduction = "harmony", dims = 1:15)
library(clustree)
MAC_DC <- FindClusters( object = MAC_DC, resolution = c(seq(.1,1,.1)))
clustree(MAC_DC@meta.data, prefix = "RNA_snn_res.")
MAC_DC <- FindClusters(MAC_DC, resolution = 0.4)
MAC_DC_marker<-FindAllMarkers(MAC_DC,only.pos = TRUE)
write.xlsx(MAC_DC_marker, file = "MAC_DC_marker_dim15_res0.4.xlsx",rowNames=TRUE)

##update cellnames
Idents(MAC_DC)<-MAC_DC$seurat_clusters
MAC_DC<- RenameIdents(MAC_DC,"0"="cDC1", "1"="migDC1","2"="Mac2", "3"="migDC2","4"="migDC3","5"="Inf.MAC1","6"="migDC5","7"="cDC4","8"="Inf.MAC2","9"="LC","10"="cDC1", "11"="cDC2",'12'="others",'13'="others", "14"="Mac1","15"="cDC5","16"="others","17"="cDC3","18"="pDC","19"="others","20"="Inf.MAC3","21"="migDC4","22"="others","23"="cDC1")
levels(MAC_DC)<-c("Mac1","Mac2","Inf.MAC1","Inf.MAC2","Inf.MAC3","cDC1","cDC2","cDC3","cDC4","cDC5","migDC1","migDC2","migDC3","migDC4","migDC5","pDC","LC","others")

##B/Mast_C
B_Mast<-subset(Infla_scRNA_harmony_res0.2,idents=c("B/Mast_C"))
meta.data<-B_Mast@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 


B_Mast<-CreateSeuratObject(counts = B_Mast@assays$RNA@counts, meta.data = meta.data)
B_Mast <- NormalizeData(B_Mast) 
B_Mast <-FindVariableFeatures(B_Mast) 
B_Mast<-ScaleData(B_Mast ) 
B_Mast<-RunPCA(B_Mast,verbose=FALSE)

B_Mast <- RunHarmony(B_Mast, group.by.vars = "group3",plot_convergence = TRUE)
##dim 1:15
B_Mast <- RunUMAP(B_Mast, reduction = "harmony", dims = 1:15)
B_Mast <- FindNeighbors(B_Mast, reduction = "harmony", dims = 1:15)
library(clustree)
B_Mast <- FindClusters(object = B_Mast, resolution = c(seq(.1,1,.1)))
clustree(B_Mast@meta.data, prefix = "RNA_snn_res.")
B_Mast <- FindClusters(B_Mast, resolution = 0.2)
B_Mast_marker<-FindAllMarkers(B_Mast,only.pos = TRUE)
write.xlsx(B_Mast_marker, file = "B_Mast_marker_dim15_res0.2.xlsx",rowNames=TRUE)

Idents(B_Mast)<-B_Mast$seurat_clusters
B_Mast<- RenameIdents(B_Mast,"0"="Plasma", "1"="others","2"="Mast", "3"="others","4"="others","5"="others","6"="others","7"="others","8"="B","9"="others","10"="others", "11"="others",'12'="others")
levels(B_Mast)<-c("Mast","B","Plasma","others")

##PC_vSMC
PC_vSMC<-subset(Infla_scRNA_harmony_res0.2,idents=c("PC_vSMC"))
meta.data<-PC_vSMC@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 

PC_vSMC<-CreateSeuratObject(counts = PC_vSMC@assays$RNA@counts,meta.data = meta.data)
PC_vSMC <- NormalizeData(PC_vSMC) 
PC_vSMC <-FindVariableFeatures(PC_vSMC) 
PC_vSMC<-ScaleData(PC_vSMC ) 
PC_vSMC<-RunPCA(PC_vSMC,verbose=FALSE)
PC_vSMC <- RunHarmony(PC_vSMC, group.by.vars = "group3",plot_convergence = TRUE)
##dim 1:15
PC_vSMC <- RunUMAP(PC_vSMC, reduction = "harmony", dims = 1:15)
PC_vSMC <- FindNeighbors(PC_vSMC, reduction = "harmony", dims = 1:15)
library(clustree)
PC_vSMC <- FindClusters( object = PC_vSMC,resolution = c(seq(.1,1,.1)))
clustree(PC_vSMC@meta.data, prefix = "RNA_snn_res.")
PC_vSMC <- FindClusters(PC_vSMC, resolution = 0.2)
PC_vSMC_marker<-FindAllMarkers(PC_vSMC,only.pos = TRUE)
library(openxlsx)
write.xlsx(PC_vSMC_marker, file = "PC_vSMC_marker_dim15_res0.2.xlsx",rowNames=TRUE)

PC_vSMC<- RenameIdents(PC_vSMC,"0"="PC1", "1"="vSMC1","2"="PC2","3"="PC3","4"="vSMC2","5"="others","6"="others", "7"="others","8"="vSMC3","9"="others","10"="others")
levels(PC_vSMC)<-c("PC1","PC2","PC3","vSMC1","vSMC2","vSMC3","others")

##MEL
MEL<-subset(Infla_scRNA_harmony_res0.2,idents=c("MEL"))
meta.data<-MEL@meta.data
meta.data<-meta.data%>% 
  dplyr::select(orig.ident,nCount_RNA,nFeature_RNA,group,group2,group3,group4,group5,group6,group7,type,type_NL) 


MEL<-CreateSeuratObject(counts = MEL@assays$RNA@counts,meta.data = meta.data)
MEL <- NormalizeData(MEL) 
MEL <-FindVariableFeatures(MEL) 
MEL<-ScaleData(MEL ) 
MEL<-RunPCA(MEL,verbose=FALSE)
MEL <- RunHarmony(MEL, group.by.vars = "group3",plot_convergence = TRUE)
##dim 1:15
MEL <- RunUMAP(MEL, reduction = "harmony", dims = 1:15)
MEL <- FindNeighbors(MEL, reduction = "harmony", dims = 1:15)
library(clustree)
MEL <- FindClusters( object = MEL, resolution = c(seq(.1,1,.1)))
clustree(MEL@meta.data, prefix = "RNA_snn_res.")
MEL <- FindClusters(MEL, resolution = 0.1)
MEL_marker<-FindAllMarkers(MEL,only.pos = TRUE)
library(openxlsx)
write.xlsx(MEL_marker, file = "MEL_marker_dim15_res0.1.xlsx",rowNames=TRUE)

MEL<- RenameIdents(MEL,"0"="MEL1", "1"="MEL2","2"="MEL3","3"="others","4"="others","5"="MEL4","6"="MEL5","7"="others","8"="Schwann","9"="Neuron")
levels(MEL)<-c("MEL1","MEL2","MEL3","MEL4","MEL5","Neuron","Schwann","others")
