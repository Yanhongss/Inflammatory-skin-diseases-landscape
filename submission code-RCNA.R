##rcna  Covarying Neighborhood Analysis
library(devtools)
library(rcna)
library(glue)
library(ggthemes)
library(patchwork)
library(ggplot2)
library(ggbeeswarm)

##Tcells for example
library(Seurat)
Idents(Tcells)<-Tcells$group7
Tcells_subset<-subset(Tcells,idents=c("type1","type1/type2","type1/type3","type2", "type2/type3","type3"))
Tcells_subset <- FindNeighbors(Tcells_subset, reduction = "harmony", dims = 1:15)

meta.data<-Tcells_subset@meta.data
meta.data$type1 <- ifelse(meta.data$group7=="type1", '1','0')
meta.data$type1_type2 <- ifelse(meta.data$group7=="type1/type2", '1','0')
meta.data$type1_type3 <- ifelse(meta.data$group7=="type1/type3", '1','0')
meta.data$type2 <- ifelse(meta.data$group7=="type2", '1','0')
meta.data$type2_type3 <- ifelse(meta.data$group7=="type2/type3", '1','0')
meta.data$type3 <- ifelse(meta.data$group7=="type3", '1','0')
meta.data$type1<-as.numeric(meta.data$type1)
meta.data$type1_type2<-as.numeric(meta.data$type1_type2)
meta.data$type1_type3<-as.numeric(meta.data$type1_type3)
meta.data$type2<-as.numeric(meta.data$type2)
meta.data$type2_type3<-as.numeric(meta.data$type2_type3)
meta.data$type3<-as.numeric(meta.data$type3)

##Tcells 
##type3
Tcells_subset@meta.data<-meta.data
table(Tcells_subset$type1)
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type3', samplem_key = 'group5', graph_use = 'RNA_nn',  verbose = TRUE,batches = NULL, covs = NULL )
options(repr.plot.width=14, repr.plot.height=4)
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type3')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')


##type1
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type1', samplem_key = 'group5', graph_use = 'RNA_nn', verbose = TRUE,batches = NULL,  covs = NULL)
options(repr.plot.width=14, repr.plot.height=4)
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type1')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')

##type2
table(Tcells_subset$group7)
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type2', samplem_key = 'group5', graph_use = 'RNA_nn', verbose = TRUE,batches = NULL, covs = NULL )
options(repr.plot.width=14, repr.plot.height=4)
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type2')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')

##type1_type2
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type1_type2', samplem_key = 'group5', graph_use = 'RNA_nn', verbose = TRUE,batches = NULL )
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type1_type2')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')

##type1_type3
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type1_type3', samplem_key = 'group5', graph_use = 'RNA_nn', verbose = TRUE,batches = NULL,covs = NULL)
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type1_type3')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')

##type2_type3
obj <- association.Seurat(seurat_object = Tcells_subset, test_var = 'type2_type3', samplem_key = 'group5', graph_use = 'RNA_nn', verbose = TRUE,batches = NULL,covs = NULL )
FeaturePlot(obj, features = c('cna_ncorrs'))[[1]] + 
  scale_color_gradient2_tableau() + 
  labs(
    title = 'CNA disease association', color = 'Correlation'
  ) + 
  FeaturePlot(obj, features = c('cna_ncorrs_fdr10'))[[1]] + 
  scale_color_gradient2_tableau() + 
  scale_color_viridis_c()+
  labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
  DimPlot(obj, group.by = 'type2_type3')[[1]] + 
  scale_color_tableau() + labs(title = 'Disease status') + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a')