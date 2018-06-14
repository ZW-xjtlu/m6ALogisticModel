#Package data future

#Group_list default

Feature_lst_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  YTHDC1_TREW = YTHDC1_TREW_gr,
  YTHDF1_TREW = YTHDF1_TREW_gr,
  YTHDF2_TREW = YTHDF2_TREW_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  miRanda = miRanda_hg19_gr,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr
)



group_list_default = list(
  UTR5 = c("UTR5", "Pos_UTR5", "length_UTR5"),
  CDS = c("CDS", "Pos_CDS", "length_CDS"),
  UTR3 = c("UTR3", "Pos_UTR3", "length_UTR3"),
  Exon = c("exons", "Pos_exons", "long_exon","Last_exons_50bp"),
  Gene = c("Pos_Tx","length_gene_ex","length_gene_full","Isoform_num","sncRNA","lncRNA","HK_genes"),
  LandMarks = c("m6Am","Start_codons","Stop_codons"),
  Motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"),
  Structure = c("struct_hybridize","struct_loop"),
  Evolution = c("PC_1nt","PC_101nt","FC_1nt","FC_101nt"),
  miRNA_RBP = c("HNRNPC_eCLIP", "YTHDC1_TREW", "YTHDF1_TREW", "YTHDF2_TREW", "miR_targeted_genes","TargetScan","Verified_miRtargets"),
  Batch = c("GC_cont_genes","GC_cont_101bp","Intercept")
)

names(group_list_default) = c("5'UTR","CDS","3'UTR","Exon","Gene","Marks","Motif","Struc","Evolution","miRNA & RBP","Batch")

readRDS("")
readRDS("")
readRDS("")

setwd("/Users/zhenwei/Documents/GitHub/m6ALogisticModel")
eCLIP_HNRNPC_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/eCLIP_HNRNPC.rds")
devtools::use_data(eCLIP_HNRNPC_gr,overwrite = TRUE)

YTHDC1_TREW_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/YTHDC1_TREW_gr.rds")
devtools::use_data(YTHDC1_TREW_gr)

YTHDF1_TREW_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/YTHDF1_TREW_gr.rds")
devtools::use_data(YTHDF1_TREW_gr)

YTHDF2_TREW_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/YTHDF2_TREW_gr.rds")
devtools::use_data(YTHDF2_TREW_gr)

miR_targeted_genes_grl <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/miR_targeted_genes_grl.rds")
devtools::use_data(miR_targeted_genes_grl)

TargetScan_hg19_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/TargetScan_hg19_gr.rds")
devtools::use_data(TargetScan_hg19_gr)

verified_targets_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/verified_targets.rds")
devtools::use_data(verified_targets_gr)


verified_targets_gr <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/I_GLM_feature_prep_2018_1_4/verified_targets.rds")

SE_Learning <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/J_GLM_Learning_2018_1_15/SE_Learning.rds")

mcols( SE_Learning ) = NULL

idx_1 <- assay(SE_Learning) == 1
idx_NA <- assay(SE_Learning) == 0
idx_n1 <- assay(SE_Learning) == -1

assay(SE_Learning)[idx_NA] = NA
assay(SE_Learning)[idx_n1] = 0
assay(SE_Learning)[idx_1] = 1

file.list <- grep("\\.",list.files("/Users/zhenwei/Documents/GitHub/TREW-cons/Targets_collection"),value = T, invert = T)

SE_Learning$ID = file.list

SE_example = SE_Learning

setwd("/Users/zhenwei/Documents/GitHub/m6ALogisticModel")
devtools::use_data(SE_example)

METTL16_CLIP <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/METTL16_CLIP.rds")
devtools::use_data(METTL16_CLIP)

ALKBH5_PARCLIP <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/PARCLIP_ALKBH5.rds")
devtools::use_data(ALKBH5_PARCLIP)

FTO_CLIP <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/CLIP_FTO.rds")
devtools::use_data(FTO_CLIP)

FTO_eCLIP <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/eCLIP_FTO.rds")
devtools::use_data(FTO_eCLIP)

METTL3_TREW <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/METTL3_TREW_gr.rds")
devtools::use_data(METTL3_TREW,overwrite = TRUE)

METTL14_TREW <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/METTL14_TREW_gr.rds")
devtools::use_data(METTL14_TREW,overwrite = TRUE)

WTAP_TREW <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/More\ Modeling\ triles/WTAP_TREW_gr.rds")
devtools::use_data(WTAP_TREW,overwrite = TRUE)

hg19_miCLIP_gr <- readRDS("/Users/zhenwei/Documents/GitHub/Project_X/hg19_miCLIP_ml.rds")
devtools::use_data(hg19_miCLIP_gr,overwrite = TRUE,pkg = "/Users/zhenwei/Documents/GitHub/m6ALogisticModel")

group_list_expanded = list(
  UTR5 = c("UTR5", "Pos_UTR5", "length_UTR5"),
  CDS = c("CDS", "Pos_CDS", "length_CDS"),
  UTR3 = c("UTR3", "Pos_UTR3", "length_UTR3"),
  Exon = c("exons", "Pos_exons", "long_exon","Last_exons_50bp"),
  Gene = c("Pos_Tx","length_gene_ex","length_gene_full","Isoform_num","sncRNA","lncRNA","HK_genes"),
  LandMarks = c("m6Am","Start_codons","Stop_codons"),
  Motif = c("AAACA","GAACA","AGACA","GGACA","AAACT","GAACT","AGACT","GGACT","AAACC","GAACC","AGACC","GGACC"),
  Structure = c("struct_hybridize","struct_loop"),
  Evolution = c("PC_1nt","PC_101nt","FC_1nt","FC_101nt"),
  miRNA_RBP = c("HNRNPC_eCLIP", "YTHDC1_TREW", "YTHDF1_TREW", "YTHDF2_TREW", "miR_targeted_genes","TargetScan","Verified_miRtargets"),
  Writers_Erasers = c("METTL3_TREW","METTL14_TREW","WTAP_TREW","METTL16_CLIP","ALKBH5_PARCLIP","FTO_CLIP","FTO_eCLIP"),
  Batch = c("GC_cont_genes","GC_cont_101bp","Intercept")
)

devtools::use_data(group_list_expanded)
