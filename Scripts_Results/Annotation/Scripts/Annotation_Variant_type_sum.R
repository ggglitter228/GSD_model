setwd("E:/GSD_VCF_model/Scripts_Results/Annotation/Scripts")

library(dplyr)
library(stringr)

VCF_GSD_anno_1 <- read.delim("../data/VCF_GSD_anno_out.hg19_multianno.txt")
VCF_GSD_anno_1 <- as.data.frame(VCF_GSD_anno_1)
VCF_GSD_anno_1 <- VCF_GSD_anno_1[,c(1:11,19,27,28,34,37,40,44,47,50,53,56,58,
                                    61,65,68,71,73,79,81,84,93,112,117,120,
                                    122,126,129,132,135,137,152,154,156,160,
                                    162,167,187)]

VCF_GSD_anno_1 <- VCF_GSD_anno_1 %>%
  mutate(
    CLNSIG = coalesce(str_match(Otherinfo11,"CLNSIG=(.+?);")[,2],"")
  )

VCF_GSD_anno_1 <- VCF_GSD_anno_1 %>%
  mutate(Pathogenic = case_when(
    CLNSIG %in% c("Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic") ~ 1,
    CLNSIG %in% c("Benign","Likely_benign","Benign/Likely_benign") ~ 0,
    TRUE ~ NA_real_      # 对于其他值设为NA
  ))
VCF_GSD_anno_1$Pathogenic <- as.character(VCF_GSD_anno_1$Pathogenic)


library(tidyr)
library(ggplot2)

variant_type_1 <- VCF_GSD_anno_1 %>%
  group_by(Func.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(Func.refGene,Pathogenic,fill = list(amount = 0)) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))


ggplot(variant_type_1,aes(x = Func.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Variant_Type_1",
       x = "variant_type",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 8,
                                    angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(variant_type_1$amount) * 1.15))




variant_type_2 <- VCF_GSD_anno_1 %>%
  group_by(ExonicFunc.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(ExonicFunc.refGene,Pathogenic,fill = list(amount = 0))


variant_type_2 <- variant_type_2 %>%
  add_row(ExonicFunc.refGene = "splicing", Pathogenic = "0", amount = 3,.after = 2) %>%
  add_row(ExonicFunc.refGene = "splicing", Pathogenic = "1", amount = 1085,.after = 3) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))

variant_type_2[1,"amount"] <- 10899
variant_type_2[2,"amount"] <- 207
variant_type_2[1:2,"ExonicFunc.refGene"] <- "Other"

ggplot(variant_type_2,aes(x = ExonicFunc.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Variant_Type_2",
       x = "variant_type",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 9,
                                   angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(variant_type_2$amount) * 1.15))



top_genes <- VCF_GSD_anno_1 %>%
  count(Gene.refGene,name = "total_count") %>%
  slice_max(total_count,n = 10) %>%
  pull(Gene.refGene)

gene_sum_1 <- VCF_GSD_anno_1 %>%
  filter(Gene.refGene %in% top_genes) %>%
  group_by(Gene.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(Gene.refGene,Pathogenic,fill = list(amount = 0)) %>%
  mutate(Gene.refGene = factor(Gene.refGene,levels = top_genes)) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))

ggplot(gene_sum_1,aes(x = Gene.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Most_Genes_Summary",
       x = "Genes",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 9,
                                   angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(gene_sum_1$amount) * 1.15))



VCF_anno_clean <- VCF_GSD_anno_1[,-c(48,49)]
VCF_anno_clean <- subset(VCF_anno_clean,AF < 0.001 | is.na(AF))

VCF_anno_clean[VCF_anno_clean == "."] <- NA
VCF_anno_clean[,11:47] <- lapply(VCF_anno_clean[,11:47],as.numeric)
VCF_anno_clean$Pathogenic <- as.factor(VCF_anno_clean$Pathogenic)

set.seed(123)

ind <- sample(2,nrow(VCF_anno_clean),replace = TRUE,prob = c(0.7,0.3))
# smple(2) 从1：2中抽样，nrow()抽样次数等于数据行数，replace = TRUE有放回抽样
#prob=c(0.7,0.3)抽样概率，1的概率70%，2的概率30%
Train_data <- VCF_anno_clean[ind == 1,]
Test_data <- VCF_anno_clean[ind == 2,]

save(VCF_GSD_anno_1,file = "../data/VCF_anno.Rdata")
save(VCF_anno_clean,file = "../data/VCF_anno_clean.Rdata")
save(Train_data,file = "../data/Train_data.Rdata")
save(Test_data,file = "../data/Test_data.Rdata")

variant_type_3 <- VCF_anno_clean %>%
  group_by(Func.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(Func.refGene,Pathogenic,fill = list(amount = 0)) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))


ggplot(variant_type_3,aes(x = Func.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Variant_Type_1",
       x = "variant_type",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 8,
                                   angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(variant_type_3$amount) * 1.15))




variant_type_4 <- VCF_anno_clean %>%
  group_by(ExonicFunc.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(ExonicFunc.refGene,Pathogenic,fill = list(amount = 0))

variant_type_4 <- variant_type_4 %>%
  add_row(ExonicFunc.refGene = "splicing", Pathogenic = "0", amount = 1,.after = 2) %>%
  add_row(ExonicFunc.refGene = "splicing", Pathogenic = "1", amount = 1007,.after = 3) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))

variant_type_4[1,"amount"] <- 1742
variant_type_4[2,"amount"] <- 182
variant_type_4[1:2,"ExonicFunc.refGene"] <- "Other"


ggplot(variant_type_4,aes(x = ExonicFunc.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Variant_Type_2",
       x = "variant_type",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 9,
                                   angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(variant_type_4$amount) * 1.15))



top_genes_1 <- VCF_anno_clean %>%
  count(Gene.refGene,name = "total_count") %>%
  slice_max(total_count,n = 10) %>%
  pull(Gene.refGene)

gene_sum_2 <- VCF_anno_clean %>%
  filter(Gene.refGene %in% top_genes) %>%
  group_by(Gene.refGene,Pathogenic) %>%
  summarise(amount = n(),.groups = 'drop') %>%
  complete(Gene.refGene,Pathogenic,fill = list(amount = 0)) %>%
  mutate(Gene.refGene = factor(Gene.refGene,levels = top_genes)) %>%
  mutate(label = ifelse(Pathogenic == "1","P/LP", "B/LB"))

ggplot(gene_sum_2,aes(x = Gene.refGene,y = amount,fill = label)) +
  geom_bar(stat = "identity",position = position_dodge(0.8),width = 0.7) +
  geom_text(aes(label = amount),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 2.5) +
  scale_fill_manual(values = c("P/LP" = "#FF6B6B","B/LB" = "#4ECDC4")) +
  labs(title = "Most_Genes_Summary",
       x = "Genes",
       y = "Freq",
       fill = "Pathogenic") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold",margin = margin(b = 15)),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(face = "italic",size = 9,
                                   angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(face = "bold",size = 10),
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)),
                     limits = c(0,max(gene_sum_2$amount) * 1.15))




