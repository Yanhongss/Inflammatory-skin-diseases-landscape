library(reshape2) 
library(vegan)
library(ggplot2)
library(stringr)
library(MASS)
library(tidyr)
library(dplyr)


Cellratio <- prop.table(table(Idents(Infla_scRNA_harmony_res0.2), Infla_scRNA_harmony_res0.2$group5), margin = 2)
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]

table(Infla_scRNA_harmony_res0.2$type,Infla_scRNA_harmony_res0.2$orig.ident)
sample <- rownames(cellper) 
group <- c(rep("NL",3),rep("infla",4),rep("NL",15),rep("infla",33),rep("infla",6),rep("Normal",63),rep("infla",7),rep("NL",12),rep("infla",20) )
samples <- data.frame(sample, group)
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']

pplist = list()
sce_groups = c("NL","infla","Normal")

##shannon
cellper_2 = cellper %>%
  dplyr::filter(group=='Normal')
cellper_2= cellper_2%>%
  dplyr::select(-group)

cellper_2 <- melt(cellper_2, id = 'sample')
erb.mat <- acast(cellper_2,  formula =   variable~sample ,  value.var = "value",  fill = 0)
Shannon.Wiener <- diversity(erb.mat, index = "shannon")
Simpson <- diversity(erb.mat, index = "simpson")
Inverse.Simpson <- diversity(erb.mat, index = "inv")
S <- specnumber(erb.mat)

J <- Shannon.Wiener/log(S)
J 