setwd("E:/GSD_VCF_model/Scripts_Results/Raw-data-collection/Scripts")

# install.packages("vcfR")
library(vcfR)
library(stringr)
library(dplyr)
library(tidyr)


#

#  -----------------------------读取分割信息，第一次筛选------------------------------------------



#读取ClinVar下载的总数据库数据
clinvar_data <- read.vcfR("../data/clinvar_20250907_GRCh37.vcf/clinvar_20250907_GRCh37.vcf")
total_raw_data_37 <- clinvar_data@fix
total_raw_data_37 <- as.data.frame(total_raw_data_37)
# 将total_raw_data_37从matrix转换成dataframe类型
INFO <- total_raw_data_37[,8]
#将INFO列有关的信息分出来，包括AlleleID，变异种类，变异审核通过标准，致病性与良性的区分等，由于数据量庞大，geneID和名称后续再进行拆分
total_raw_data_37 <- total_raw_data_37 %>%   #%>%管道操作符，将左侧对象的值传递给右侧函数的第一个参数位置
  mutate(
    #mutate是dplyr包中的函数，用于添加新列或修改现有列
    ALLELEID = coalesce(str_match(INFO, "ALLELEID=(\\d+);")[, 2], ""),
    CLNDISDB = coalesce(str_match(INFO,"CLNDISDB=(.+?);")[,2],""),
    CLNDN = coalesce(str_match(INFO,"CLNDN=(.+?);")[,2],""),
    CLNHGVS = coalesce(str_match(INFO, "CLNHGVS=(.+?);")[, 2], ""),
    CLNREVSTAT = coalesce(str_match(INFO,"CLNREVSTAT=(.+?);")[,2], ""),
    CLNSIG = coalesce(str_match(INFO, "CLNSIG=(.+?);")[, 2], ""),
    CLNVC = coalesce(str_match(INFO, "CLNVC=(\\w+);")[, 2], "")
    # GEN_name = coalesce(str_match(INFO,"GENEINFO=(\\w+):(\\d+);")[,2],""),
    # GEN_id = coalesce(str_match(INFO,"GENEINFO=(\\w+):(\\d+);")[,3],"")
  )
# coalesce()是dplyr包中的函数，用于从多个向量中选取第一个非缺失值。它逐个元素地检查提供的向量，返回每个位置上的第一个非NA值
# str_match是stringr包中的一个函数，用于从字符串中提取正则表达式匹配的捕获组。它返回一个字符矩阵，其中包含匹配的完整字符串和各个捕获组的内容
# \\d - 匹配任何数字[0-9];\\w - 匹配任何字母、数字或下划线；. - 匹配任何字符（除了换行符）
# + 匹配前面元素一次或多次；? 匹配前面的元素零次或一次（是匹配变为非贪婪匹配。）



df <- head(total_raw_data_37,200)    #检查数据情况使用


#筛选出单核苷酸变异，Pathogenic/Likely Pathogenic，Benign/Likely Benign，两星以上评级的变异条目
clinvar_vcf_37 <- subset(total_raw_data_37,CLNVC == "single_nucleotide_variant")
clinvar_vcf_37 <- subset(clinvar_vcf_37, CLNSIG %in% c("Pathogenic", "Likely_pathogenic",
                                                                             "Pathogenic/Likely_pathogenic", "Benign",
                                                                             "Likely_benign", "Benign/Likely_benign"))
clinvar_vcf_37 <- subset(clinvar_vcf_37, CLNREVSTAT %in% c("criteria_provided,_multiple_submitters,_no_conflicts",
                                                                     "reviewed_by_expert_panel","practice_guideline"))




#  ---------------------------分割gene信息，筛选与GSD-gene相关的--------------------------------------------


#将INFO列中geneID以及geneSymol拆分出来

INFO_filter <- clinvar_vcf_37[,8]
geneinfo_obj = str_match(INFO_filter,"GENEINFO=(.+?);")
geneinfo_str = geneinfo_obj[,2]


geneinfo_name <- c()
geneinfo_id <- c()
for (gene_str in geneinfo_str) {
  # print(gene_str)
  if (grepl("\\|",gene_str)) {  #要用双反斜杠来表示转义字符
    geneinfo_split <- unlist(strsplit(gene_str,"\\|"))
    geneinfo_name_list <- c()
    geneinfo_id_list <- c()
    for (geneinfo in geneinfo_split) {
      # print(geneinfo)
      parts <- unlist(strsplit(geneinfo,":"))
      geneinfo_name_list <- c(geneinfo_name_list,parts[1])
      geneinfo_id_list <- c(geneinfo_id_list,parts[2])
    }
    # print(geneinfo_name_list)
    geneinfo_name <- append(geneinfo_name,paste(geneinfo_name_list,collapse = ";"))
    geneinfo_id <- append(geneinfo_id,paste(geneinfo_id_list,collapse = ";"))
  } else {
    # print(gene_str)
    parts <- unlist(strsplit(gene_str,":"))
    # print(parts)
    geneinfo_name <- append(geneinfo_name,parts[1])
    geneinfo_id <- append(geneinfo_id,parts[2])
  }
}

#加入到列表中。
clinvar_vcf_37 <- clinvar_vcf_37 %>%
  mutate(
    Symbol = coalesce(geneinfo_name,""),
    EntrezID = coalesce(geneinfo_id,"")
  )


#读取遗传性骨病（Genetic Skeletal Dieases，GSD）的基因-表型对照表
library(readxl)
library(tidyverse)
library(writexl)
library(tibble)

GSD_gene_PHO <- read_excel("../data/GSD_Gene_PHO_simple.xlsx")
clinvar_vcf_filter <- clinvar_vcf_37 %>% separate_rows(Symbol,sep = ";")
GSD_gene_PHO <- GSD_gene_PHO %>% separate_rows(`Gene or locus`,sep = ",")
GSD_gene_PHO <- GSD_gene_PHO %>% separate_rows(`MIM No.`, sep = ",")



VCF_GSD <- clinvar_vcf_filter %>% mutate(
  if_OMIM = str_detect(CLNDISDB,"OMIM:\\d+"),
  OMIM_codes = ifelse(if_OMIM,
                      map_chr(str_extract_all(CLNDISDB,"OMIM:\\d+"),
                              ~ paste(str_remove(.x,"OMIM:"), collapse = ";")),
                      NA_character_),
  if_OMIM = NULL
)
# mutate()在数据框中创建或者修改新列；str_detect()检测字符串中是否含有匹配模式，是返回TRUE，否返回FALSE；
# if_else(condition,yes_value,no_value)；str_extract_all()从字符串中提取所有匹配正则表达式的部分；
# map_chr()对列表的每个元素应用函数，确保返回字符向量。
# ~paste(str_remove(.x,"OMIM:"), collapse = ";") 等价于：function(.x) {paste(str_remove(.x, "OMIM:"), collapse = ";")}
# ~-purrr公式语法（匿名函数简写）~expression(.x) 等价于 function(.x){expression(.x)}
# str_remove()移除匹配部分；paste(...,sep = "", collapse = NULL)，...是一个或多个要连接的对象，sep连接多个向量时使用的分隔符，collapse连接向量内多个元素时使用的分隔符。

VCF_GSD_PHO <- VCF_GSD %>%
  separate_rows(OMIM_codes,sep = ";") %>%
  mutate(OMIM = trimws(OMIM_codes))    # trimws()去除字符串两端的空格
VCF_GSD_PHO_na <- VCF_GSD_PHO %>% filter(!is.na(OMIM_codes))

GSD_gene_PHO <- rename(GSD_gene_PHO,Symbol = `Gene or locus`,OMIM = `MIM No.`)

VCF_GSD_PHO_match <- inner_join(VCF_GSD_PHO_na, GSD_gene_PHO,by = c("Symbol","OMIM"), relationship = "many-to-many") %>% distinct()
# inner_join()内连接，只保留两个数据集中都存在的行；distinct()去除重复行
VCF_GSD_PHO_match <- VCF_GSD_PHO_match %>% group_by(across(-OMIM)) %>%
  summarise(OMIM = paste(unique(OMIM),collapse = ";"), .groups = "drop")
# group_by()：按指定列分组；across(-OMIM)按除了OMIM列之外的所有列分组；summarise()对每组数据进行汇总
# paste()用分号链接多个OMIM编号；.group = "drop" 取消分组。



#  -----------------------添加没有OMIM表型的B/LB的变异------------------------------------------------


VCF_B <- VCF_GSD_PHO %>% filter(CLNSIG %in% c("Benign/Likely_benign","Benign","Likely_benign"))
VCF_B <- VCF_B %>% filter(is.na(OMIM))
symbol <- c(unique(GSD_gene_PHO$Symbol))

VCF_B <- VCF_B %>% filter(Symbol %in% symbol)

final_VCF_GSD <- rbind(VCF_GSD_PHO_match[,c(1:17,23)],VCF_B[,c(1:17,19)])

final_VCF_GSD <- final_VCF_GSD[!duplicated(final_VCF_GSD$ID),]


final_VCF_GSD$ALT[is.na(final_VCF_GSD$ALT)] <- "N"
if (!is.numeric(final_VCF_GSD$QUAL)) {
  final_VCF_GSD$QUAL <- as.numeric(final_VCF_GSD$QUAL)
}
final_VCF_GSD$QUAL[is.na(final_VCF_GSD$QUAL)] <- 0
final_VCF_GSD$FILTER[is.na(final_VCF_GSD$FILTER)] <- "unknown"



save(total_raw_data_37,clinvar_vcf_37,clinvar_vcf_filter,GSD_gene_PHO,VCF_GSD,VCF_GSD_PHO,file = "../data/Clinvar_filter_1.Rdata")
save(final_VCF_GSD, file = "../data/VCF_GSD.Rdata")


# - --------------------统计input的数据情况------------------------------------


#统计其中不同类型的变异分别有多少
MC <- str_match(final_VCF_GSD$INFO, "MC=(.+?);")[,2]
VCF_GSD_MC <- final_VCF_GSD %>% mutate(
  variants_split  = str_split(MC, ","),
  variant1 = map_chr(variants_split, ~ifelse(length(.x) > 0,
                                             str_extract(.x[1],"(?<=\\|)[^|]+$"),
                                             NA)),
  variant2 = map_chr(variants_split, ~ifelse(length(.x) > 1,
                                             str_extract(.x[2],"(?<=\\|)[^|]+$"),
                                             NA))
) %>%
  select( -variants_split)

save(VCF_GSD_MC, file = "../data/VCF_GSD_MC.Rdata")

VCF_GSD_missense  <- subset(VCF_GSD_MC, variant1 == "missense_variant" | variant2 == "missense_variant")
VCF_GSD_missense_noncoding <- subset(VCF_GSD_missense, !is.na(VCF_GSD_missense$variant1) & !is.na(VCF_GSD_missense$variant2))

save(VCF_GSD_missense,file = "../data/VCF_GSD_missense.Rdata")


variant_sum <- as.data.frame.matrix(table(VCF_GSD_MC$variant1,VCF_GSD_MC$variant2))
variant_sum <- rownames_to_column(variant_sum,var = "variant_type")
write_xlsx(variant_sum,"../data/variant_sum.xlsx", col_names = TRUE)



# - ----------------------将最终文件转换为VCF文件，方便AnnoVar注释-------------------------------------------------



#将文件转换成vcf格式输出方便进行注释。
options("repo" = c(CRAN ="HTTPS://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
library("VariantAnnotation")

VCF_GSD_anno_row <- VRanges(
  seqnames = final_VCF_GSD$CHROM,
  ranges = IRanges(start = as.numeric(final_VCF_GSD$POS),width = 1),
  ref = final_VCF_GSD$REF,
  alt = final_VCF_GSD$ALT
)
# VRanges是GenomicRanges的扩展，专门用于基因变异数据，ranges必须使用IRanges定义变异位置。width是编译影响的长度。



VCF_GSD_anno <- VCF(
  rowRanges = VCF_GSD_anno_row,
  fixed = DataFrame(
    REF = DNAStringSet(final_VCF_GSD$REF),
    ALT = DNAStringSetList(strsplit(final_VCF_GSD$ALT,",")),
    QUAL = final_VCF_GSD$QUAL,
    FILTER = final_VCF_GSD$FILTER
  ),
  info = DataFrame(INFO = final_VCF_GSD$INFO)
)
# VCF类是VariantAnnotation包中表示VCF文件数据的核心容器。
# rowRanges：包含变异位置信息的VRanges对象；fixed固定字段DataFrame，包含VCF文件每行都有的
# 固定字段。DNAStringSet包装，确保为有效的DNA序列。

header(VCF_GSD_anno) <- VCFHeader(
  samples = c("GSD gene-related VCF to annotate"),
  reference = "hg19"
)

writeVcf(VCF_GSD_anno,"../data/VCF_GSD_anno.vcf")
