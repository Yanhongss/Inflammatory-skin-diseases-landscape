setwd("/data/ext/syh")
palette1=c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
           "#E15759","#FF9D9A","#79706E","#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6",
           '#8dd3c7','#fee391','#bebada','#fb8072','#d9d9d9','#fdb462','#b3de69','#80b1d3','#bc80bd',
           '#ccebc5','#ffed6f','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f')

Infla_scRNA_harmony_res0.2_metabolic <- readRDS("/data/syh/AA_new/Infla_scRNA_harmony_res0.2_metabolic.rds")

##AA
AA<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("AA"))
expr<-AA@assays$RNA
meta<-AA@meta.data
table(meta$group2)
rogue.res.AA<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.AA)
rogue_df_sample <- data.frame(Sample = names(rogue.res.AA),
                              ROGUE = unlist(rogue.res.AA))
ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#AA_NL
AA_NL<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("AA_NL"))
expr<-AA_NL@assays$RNA
meta<-AA_NL@meta.data
table(meta$group2)
rogue.res.AA_NL<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.AA_NL)
rogue_df_sample <- data.frame(Sample = names(rogue.res.AA_NL),
                              ROGUE = unlist(rogue.res.AA_NL))
ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#AD
AD<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("AD"))
expr<-AD@assays$RNA
meta<-AD@meta.data
table(meta$group2)
rogue.res.AD<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.AD)
rogue_df_sample <- data.frame(Sample = names(rogue.res.AD),
                              ROGUE = unlist(rogue.res.AD))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8)

#AD_NL
AD_NL<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("AD_NL"))
expr<-AD_NL@assays$RNA
meta<-AD_NL@meta.data
table(meta$group5)
rogue.res.AD_NL<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.AD_NL)
rogue_df_sample <- data.frame(Sample = names(rogue.res.AD_NL),
                              ROGUE = unlist(rogue.res.AD_NL))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#HS
HS<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("HS"))
expr<-HS@assays$RNA
meta<-HS@meta.data
table(meta$group2)
rogue.res.HS<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.HS)
rogue_df_sample <- data.frame(Sample = names(rogue.res.HS),
                              ROGUE = unlist(rogue.res.HS))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#PN
PN<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("PN"))
expr<-PN@assays$RNA
meta<-PN@meta.data
table(meta$group2)
rogue.res.PN<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.PN)
rogue_df_sample <- data.frame(Sample = names(rogue.res.PN),
                              ROGUE = unlist(rogue.res.PN))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#Psoriasis
Psoriasis<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("Psoriasis"))
expr<-Psoriasis@assays$RNA
meta<-Psoriasis@meta.data
table(meta$group2)
rogue.res.Psoriasis<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.Psoriasis)
rogue_df_sample <- data.frame(Sample = names(rogue.res.Psoriasis),
                              ROGUE = unlist(rogue.res.Psoriasis))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette2)

#Psoriasis_NL
table(Infla_scRNA_harmony_res0.2_metabolic$type_NL)
Psoriasis_NL<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("Psoriasis_NL"))
expr<-Psoriasis_NL@assays$RNA
meta<-Psoriasis_NL@meta.data
table(meta$group2)
rogue.res.Psoriasis_NL<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.Psoriasis_NL)
rogue_df_sample <- data.frame(Sample = names(rogue.res.Psoriasis_NL),
                              ROGUE = unlist(rogue.res.Psoriasis_NL))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#Normal
Normal<-subset(Infla_scRNA_harmony_res0.2_metabolic,idents=c("Normal"))
expr<-Normal@assays$RNA
meta<-Normal@meta.data
table(meta$group5)
table(meta$group3)

table(meta$cell.type,meta$group5)

rogue.res.Normal<- rogue(expr, labels = meta$cell.type, samples =meta$group5,platform = "UMI")

rogue.boxplot(rogue.res.Normal)
rogue_df_sample <- data.frame(Sample = names(rogue.res.Normal),
                              ROGUE = unlist(rogue.res.Normal))

ggplot(rogue_df_sample, aes(x = Sample, y = ROGUE, fill = Sample)) +
  theme_minimal() +
  labs(title = "Sample Consistency (ROGUE)", 
       x = "Sample", 
       y = "ROGUE Value")+
  geom_col(position = "dodge")+
  stat_compare_means(method = "anova")+
  geom_point(aes(color = Sample),  # Add points for individual data
             position = position_jitter(width = 0.15, height = 0), 
             size = 2, alpha = 0.8) +
  scale_fill_manual(values = palette1)

#coefficient of variation calculation
cv <- function(x) {
  return(sd(x,na.rm=TRUE) / mean(x,na.rm=TRUE) * 100)
}

##here
x <- matrix(NA, nrow = 3, ncol = length(rogue.res.HS))

# 循环计算
for (i in 1:11) {
  
  # 提取数据
  data <- rogue.res.HS[[i]]
  
  # 计算标准差、均值和变异系数
  x[1, i] <- sd(data, na.rm = TRUE)
  x[2, i] <- mean(data, na.rm = TRUE)
  x[3, i] <- cv(data)
}

# results
x
colnames(x)<-colnames(rogue.res.HS)
rownames(x)<-c("sd",'mean','cv')
x
x<-data.frame(x)

