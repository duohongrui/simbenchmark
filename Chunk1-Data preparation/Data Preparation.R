################################################################################
#    &&&....&&&    % Project: A Systematic Evaluation of Single-Cell RNA-Seq   #
#  &&&&&&..&&&&&&                       Data Simulation Methods                #
#  &&&&&&..&&&&&&  % Chunk1: Data Preparation                                  #
#  &&&&&&&&&&&&&&  % Author: Duo Hongrui                                       #
#   &&&&&&&&&&&&   % Date: Oct. 26th, 2021                                     #
#     &&&&&&&&     % Environment: R version 4.1.1 (2021-08-10);                #
#       &&&&       % Platform: x86_64-w64-mingw32/x64 (64-bit)                 #
#        &         %           Intel(R)Core(TM)i5-8250U CPU@1.60GHz 1.80 GHz   #                                                                                                       #
################################################################################
library(Seurat)
library(dplyr)
library(stringr)
library(R.utils)
library(biomaRt)
library(GEOquery)

merged_fun <- function(dir, list){
  data_comb <- read.table(file = paste0(dir, '/', list[1]),
                          header = TRUE,
                          row.names = 1,
                          sep = "\t")
  for(file in list[2:length(list)]){
    data_tmp <- read.table(file = paste0(dir, '/', file),
                           header = TRUE,
                           row.names = 1,
                           sep = "\t")
    data_comb <- cbind(data_comb, data_tmp)
    message('Complete!')
  }
  print(dim(data_comb))
  return(data_comb)
}



################################################################################

#                         Download datasets from GEO

################################################################################
dir.create(path = './datasets')
setwd('./datasets/')
GEO_list <- c("GSE102827", "GSE106202", "GSE108097", "GSE95432",
              "GSE109205", "GSE112004", "GSE114727", "GSE117450", "GSE117617",
              "GSE133539", "GSE133542", "GSE167297", "GSE183590", "GSE184950",
              "GSE54006", "GSE54695", "GSE60361", "GSE62270", "GSE65525",
              "GSE67835", "GSE72857", "GSE78779", "GSE85241", "GSE85755",
              "GSE86469", "GSE87038", "GSE90047", "GSE95434", "GSE95436",
              "GSE95445", "GSE95448", "GSE95701", "GSE98816")
for(accession_number in GEO_list){
  getGEOSuppFiles(GEO = accession_number)
}



################################################################################

#                       Download datasets from FigShare

################################################################################
dir.create(path = './FigShare_10X')
dir.create(path = './FigShare_Smartseq2')
a <- reticulate::import("wget")
a$download(url = "https://figshare.com/ndownloader/files/10700167",
           out = "./FigShare_10X/droplet.zip")

a$download(url = "https://figshare.com/ndownloader/files/10700143",
           out = "./FigShare_Smartseq2/Smartseq2.zip")
## Downlaod cell annotation file
a$download(url = "https://figshare.com/ndownloader/files/10881908",
           out = "./FigShare_Smartseq2/annotations_FACS.csv")
a$download(url = "https://figshare.com/ndownloader/files/10881902",
           out = "./FigShare_10X/annotations_droplet.csv")



################################################################################

#                       Download datasets from zenodo

################################################################################
dir.create(path = './zenodo_datasets')
traj_datasets <- c("aging-hsc-old_kowalczyk.rds",               
                   "aging-hsc-young_kowalczyk.rds",             
                   "cell-cycle_buettner.rds",                   
                   "cellbench-SC1_luyitian.rds",                
                   "cellbench-SC2_luyitian.rds",                
                   "cellbench-SC3_luyitian.rds",                
                   "cellbench-SC4_luyitian.rds",                
                   "developing-dendritic-cells_schlitzer.rds",  
                   "germline-human-female-weeks_li.rds",        
                   "germline-human-male-weeks_li.rds",          
                   "hematopoiesis-gates_olsson.rds",            
                   "human-embryos_petropoulos.rds",             
                   "macrophage-salmonella_saliba.rds",          
                   "mesoderm-development_loh.rds",              
                   "myoblast-differentiation_trapnell.rds",     
                   "NKT-differentiation_engel.rds",             
                   "pancreatic-alpha-cell-maturation_zhang.rds",
                   "pancreatic-beta-cell-maturation_zhang.rds", 
                   "psc-astrocyte-maturation-glia_sloan.rds",   
                   "psc-astrocyte-maturation-neuron_sloan.rds", 
                   "stimulated-dendritic-cells-LPS_shalek.rds", 
                   "stimulated-dendritic-cells-PAM_shalek.rds", 
                   "stimulated-dendritic-cells-PIC_shalek.rds")
traj_datasets_url <- paste0("https://zenodo.org/record/1443566/files/real/gold/",
                            traj_datasets,
                            "?download=1")
for(i in 1:2){
  download.file(url = traj_datasets_url[i],
                destfile = paste0("./zenodo_datasets/","data", i+78, "_", traj_datasets[i]))
}


################################################################################

#            Preprocess the count matrix and convert to list object

################################################################################
dir.create(path = './preprocessed_data')
setwd("./datasets/")


########## data1 GSE54006
data <- read.table(file = './GSE54006/GSE54006_umitab.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
## cell information
cell_condition <- stringr::str_split(colnames(data), pattern = "_", simplify = TRUE)[, 1]
group1 <- paste0("X", 0:27)
group2 <- paste0("X", c(28, 28, 32, 33, 36, 37, 40, 41))
### save two groups
save_index <- ifelse(cell_condition %in% group1, FALSE, TRUE)
data <- data[, save_index]
cell_condition <- cell_condition[save_index]
### group
group_condition <- ifelse(cell_condition %in% group2, 1, 2)
treatment <- ifelse(group_condition == 1, "PBS", "LPS")
## Filter
index <- colSums(data) > 0
data <- data[, index]
group_condition <- group_condition[index]
treatment <- treatment[index]
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
colnames(data) <- paste0(treatment, 1:ncol(data))
## ERCC count
ERCC_count <- data[grep(rownames(data), pattern = "^ERCC-"), ]
## Gene transformation
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "transcript_biotype",
                                            "external_gene_name"),
                             mart = ensembl) %>%
  dplyr::filter("external_gene_name" != "")
gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
data <- data[gene_filter, ]
rownames(data) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
## combine
data <- rbind(data, ERCC_count)
## data info
data_info <- simutils::meta_info(id = "data1_GSE54006",
                                 repository = "GEO",
                                 accession_number = "GSE54006",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Spleen",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 ERCC = TRUE,
                                 dilution_factor = 50000,
                                 volume = 0.01,
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data1_GSE54006 <- list(data = as.matrix(data),
                       data_info = data_info)
saveRDS(data1_GSE54006, file = '../preprocessed_data/data1_GSE54006.rds')



########## data2 GSE72857
data <- read.table(file = './GSE72857/GSE72857_umitab.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
## rownames
gene_name <- str_split(rownames(data), pattern = ";", simplify = TRUE)[, 1]
rownames(data) <- gene_name
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
data_info <- simutils::meta_info(id = "data2_GSE72857",
                                 repository = "GEO",
                                 accession_number = "GSE72857",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
## Save
data2_GSE72857 <- list(data = as.matrix(data),
                       data_info = data_info)
saveRDS(data2_GSE72857, file = '../preprocessed_data/data2_GSE72857.rds')



########## data3-data4 GSE133542
data_list <- list.files(path = './GSE133542/', pattern = 'mat.tsv.gz')
Organism <- c("Homo sapiens", "Mus musculus")
organs <- c("Human PBMC", "Mouse Colon Cell Lines")
for(i in 1:length(data_list)){
  data <- read.table(file = paste0('./GSE133542/', data_list[i]),
                     header = TRUE,
                     row.names = 1,
                     sep = '\t')
  ## Filter
  data <- data[, colSums(data) > 0]
  print(table(colSums(data) > 0))
  ## rownames
  gene_name <- str_split(rownames(data), pattern = "[.]", simplify = TRUE)[, 1]
  data$gene_name <- gene_name
  data <- aggregate(.~data$gene_name, data[, -ncol(data)], mean)
  rownames(data) <- data$`data$gene_name`
  data <- round(data[, -1])
  rownames(data) <- str_replace_all(rownames(data),
                                    pattern = "_",
                                    replacement = "-")
  print(dim(data))
  data_info <- simutils::meta_info(id = paste0("data", i+2, '_GSE133542_subset', i),
                                   repository = "GEO",
                                   accession_number = "GSE133542",
                                   platform = "MARS-Seq",
                                   species = Organism[i],
                                   organ = organs[i],
                                   cell_num = ncol(data),
                                   gene_num = nrow(data))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', 'data', i+2, '_GSE133542_subset', i, '.rds'))
}



########## data5 SCP1729
data <- read.table("./SCP1729/rawData.txt",
                   header = TRUE,
                   row.names = 1)
cell_info <- read.table("./SCP1729/metadataInfo.txt",
                        header = TRUE,
                        row.names = 1, sep = "\t")
cell_info <- cell_info[-1, ]
data <- data[, rownames(cell_info)]
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## Treatment
treatment <- cell_info$age
## Group
group_condition <- ifelse(treatment == "Young", 1, 2)
## Filter
index <- colSums(data) > 0
table(index)
## Data meta
data_info <- simutils::meta_info(id = "data5_SCP1729",
                                 repository = "Single Cell Portal",
                                 accession_number = "SCP1729",
                                 platform = "10X Genomics",
                                 species = "Mus musculus",
                                 organ = "Mammary tissue",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data5_SCP1729.rds")



########## data6 SCP1821
data <- read.csv("./SCP1821/scp_raw_counts.csv",
                 header = TRUE,
                 row.names = 1)
cell_info <- read.csv("./SCP1821/scp_LEC_metadata.csv",
                      header = TRUE,
                      row.names = 1)
cell_info <- cell_info[-1, ]
data <- data[, rownames(cell_info)]
## Treatment
treatment <- cell_info$genotype
## Group
group_condition <- ifelse(treatment == "WT", 1, 2)
## cluster
cluster <- cell_info$LECsubtype_annotation
## Filter
index <- colSums(data) > 0
table(index)
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## Data meta
data_info <- simutils::meta_info(id = "data6_SCP1821",
                                 repository = "Single Cell Portal",
                                 accession_number = "SCP1821",
                                 platform = "10X Genomics",
                                 species = "Mus musculus",
                                 organ = "Mesentery",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment,
                                 cluster_info = cluster)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data6_SCP1821.rds")



########## data7 SCP1675
data <- read.csv("./SCP1675/rawls5532_countmatrix.csv",
                 header = TRUE,
                 row.names = 1)
cell_info <- read.table("./SCP1675/metaData_SCP.txt",
                        header = TRUE,
                        row.names = 1,
                        sep = "\t")
cell_info <- cell_info[-1, ]
data <- data[, rownames(cell_info)]
## Treatment
treatment <- cell_info$biosample_id
## Group
group_condition <- ifelse(treatment == "WT", 1, 2)
## Filter
index <- colSums(data) > 0
table(index)
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## Data meta
data_info <- simutils::meta_info(id = "data7_SCP1675",
                                 repository = "Single Cell Portal",
                                 accession_number = "SCP1675",
                                 platform = "10X Genomics",
                                 species = "Danio rerio",
                                 organ = "Larvae",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data7_SCP1675.rds")



########## data8 GSE95701
untar('./GSE95701/GSE95701_RAW.tar', exdir = './GSE95701/')
data_list <- list.files(path = './GSE95701/', pattern = 'txt.gz')[2:6]
## Merge data
data <- merged_fun(dir = './GSE95701', list = data_list)
treatment <- c(rep("case", 384*2),
               rep("control", 384*3))
group_condition <- ifelse(treatment == "case", 2, 1)
## Filter
index <- colSums(data) > 0
table(index)
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## Data meta
data_info <- simutils::meta_info(id = "data8_GSE95701",
                                 repository = "GEO",
                                 accession_number = "GSE95701",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Monocytes",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = paste0('../preprocessed_data/', 'data8_GSE95701','.rds'))



########## data9 SCP1618
data <- as.data.frame(Read10X("./SCP1618/seurat_file/"))
cell_info <- read.table("./SCP1618/metaData6.tsv",
                        header = TRUE,
                        row.names = 1,
                        sep = "\t")
cell_info <- cell_info[-1, ]
data <- data[, rownames(cell_info)]
## Treatment
treatment <- cell_info$condition
## Group
group_condition <- ifelse(treatment == "Control", 1, 2)
## Filter
index <- colSums(data) > 0
table(index)
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## Data meta
data_info <- simutils::meta_info(id = "data9_SCP1618",
                                 repository = "Single Cell Portal",
                                 accession_number = "SCP1618",
                                 platform = "10X Genomics",
                                 species = "Mus musculus",
                                 organ = "Retinae",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data9_SCP1618.rds")



########## data10 GSE98816
data <- read.table(file = './GSE98816/GSE98816_Brain_samples_raw_read_counts_matrix.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
data_info <- simutils::meta_info(id = "data10_GSE98816",
                                 repository = "GEO",
                                 accession_number = "GSE98816",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "Brain",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data10_GSE98816.rds')



########## data11 FigShare Smartseq2 and 10× batch (Aorta and Bladder)
untar('./FigShare_Smartseq2/Smartseq2.zip', exdir = './FigShare_Smartseq2')
untar('./FigShare_10X/droplet.zip', exdir = './FigShare_10X')
## smart-seq2 data list
cell_anno_smartseq2 <- read.csv("./FigShare_Smartseq2/annotations_FACS.csv")
data_list_smartseq2 <- list.files(path = './FigShare_Smartseq2/FACS/', pattern = '.csv')
## 10× data list
cell_anno_10X <- read.csv("./FigShare_10X/annotations_droplet.csv")
cell_anno_10X$cell <- stringr::str_split(cell_anno_10X$cell,
                                         pattern = "_",
                                         simplify = TRUE)[, 4]
data_list_10X <- list.files(path = './FigShare_10X/droplet/')
## Gene transformation data
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "transcript_biotype",
                                            "external_gene_name"),
                             mart = ensembl) %>%
  dplyr::filter("external_gene_name" != "")
## Aorta
Figshare_function <- function(organs,
                              cell_type_index,
                              droplet_file_index,
                              droplet_organs,
                              data_id){
  data1_smartseq2 <- read.csv(file = file.path("./FigShare_Smartseq2/FACS",
                                               data_list_smartseq2[str_which(data_list_smartseq2, pattern = organs)]),
                              header = TRUE,
                              row.names = 1)
  ## Read 10× data
  data1_file <- data_list_10X[droplet_file_index]
  data1_10X <- as.data.frame(Seurat::Read10X(file.path("./FigShare_10X/droplet", data1_file[1])))
  if(length(data1_file) != 1){
    for(i in data1_file[2:length(data1_file)]){
      data_tmp <- as.data.frame(Seurat::Read10X(file.path("./FigShare_10X/droplet", i)))
      data1_10X <- cbind(data1_10X, data_tmp)
    }
  }
  colnames(data1_10X) <- stringr::str_split(colnames(data1_10X),
                                            pattern = "-",
                                            simplify = TRUE)[, 1]
  anno_data1_smartseq2 <- cell_anno_smartseq2 %>% 
    filter(tissue == organs) %>% 
    filter(cell %in% intersect(colnames(data1_smartseq2), cell_anno_smartseq2$cell))
  anno_data1_10X <- cell_anno_10X %>% 
    filter(tissue == droplet_organs) %>% 
    filter(cell %in% intersect(colnames(data1_10X), cell))
  cell_type_smartseq2 <- unique(anno_data1_smartseq2$cell_ontology_class)
  cell_type_10X <- unique(anno_data1_10X$cell_ontology_class)
  cell_type <- intersect(cell_type_smartseq2, cell_type_10X)[seq_len(cell_type_index)]
  anno_data1_smartseq2 <- anno_data1_smartseq2 %>% 
    filter(cell_ontology_class %in% cell_type)
  data1_smartseq2 <- data1_smartseq2[, anno_data1_smartseq2$cell]
  #### Filter
  index <- colSums(data1_smartseq2) > 0
  data1_smartseq2 <- data1_smartseq2[, index]
  anno_data1_smartseq2 <- anno_data1_smartseq2[index, ]
  anno_data1_smartseq2 <- anno_data1_smartseq2$cell_ontology_class
  #### ERCC count
  ERCC_count <- data1_smartseq2[grep(rownames(data1_smartseq2), pattern = "^ERCC-"), ]
  #### gene transformation
  gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data1_smartseq2),
                                                                    id_convert$external_gene_name))]
  data1_smartseq2 <- data1_smartseq2[gene_filter, ]
  rownames(data1_smartseq2) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data1_smartseq2),
                                                                               id_convert$external_gene_name))]
  data1_smartseq2 <- rbind(data1_smartseq2, ERCC_count)
  message("Smart-seq2 data done!")
  
  ## Filter 10X
  anno_data1_10X <- anno_data1_10X %>% 
    filter(cell_ontology_class %in% cell_type)
  data1_10X <- data1_10X[, anno_data1_10X$cell]
  print(table(anno_data1_10X$cell_ontology_class))
  #### Filter
  index <- colSums(data1_10X) > 0
  data1_10X <- data1_10X[, index]
  anno_data1_10X <- anno_data1_10X[index, ]
  anno_data1_10X <- anno_data1_10X$cell_ontology_class
  #### ERCC count
  ERCC_count <- data1_10X[grep(rownames(data1_10X), pattern = "^ERCC-"), ]
  #### gene transformation
  gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data1_10X),
                                                                    id_convert$external_gene_name))]
  data1_10X <- data1_10X[gene_filter, ]
  rownames(data1_10X) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data1_10X),
                                                                         id_convert$external_gene_name))]
  data1_10X <- rbind(data1_10X, ERCC_count)
  message("10× Genomics data done!")
  
  ## Combination
  data <- cbind(data1_smartseq2,
                data1_10X)
  cluster <- c(anno_data1_smartseq2,
               anno_data1_10X)
  batch_info <- c(rep("Smart-seq2", ncol(data1_smartseq2)),
                  rep("10×", ncol(data1_10X)))
  data_id <- paste0(data_id, "_Figshare_", organs[1])
  data_info <- simutils::meta_info(id = data_id,
                                   repository = "FigShare",
                                   URL = c("https://figshare.com/ndownloader/files/10700143",
                                           "https://figshare.com/ndownloader/files/10700167"),
                                   platform = c("Smart-seq2", "10× Genomics"),
                                   species = "Mus musculus",
                                   organ = organs,
                                   cell_num = ncol(data),
                                   gene_num = nrow(data),
                                   ERCC = TRUE,
                                   dilution_factor = 600000,
                                   volume = 0.4,
                                   batch_info = batch_info,
                                   cluster_info = cluster)
  print(dim(data))
  print(length(table(data_info$cluster_info)))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', data_id,'.rds'))
}

Figshare_function(organs = "Aorta",
                  cell_type_index = 2,
                  droplet_file_index = 4,
                  droplet_organs = c("Heart_and_Aorta"),
                  data_id = "data11")

########## data12 FigShare Smartseq2 and 10× batch (Bladder)
Figshare_function(organs = "Bladder",
                  cell_type_index = 2,
                  droplet_file_index = 1:3,
                  droplet_organs = "Bladder",
                  data_id = "data12")


########## data13 FigShare Smartseq2 and 10× batch (Lung)
Figshare_function(organs = "Lung",
                  cell_type_index = 2,
                  droplet_file_index = 13,
                  droplet_organs = "Lung",
                  data_id = "data13")

########## data14 FigShare Smartseq2 and 10× batch (Marrow)
Figshare_function(organs = "Marrow",
                  cell_type_index = 2,
                  droplet_file_index = c(19, 20),
                  droplet_organs = "Marrow",
                  data_id = "data14")

########## data15 FigShare Smartseq2 and 10× batch (Tongue)
Figshare_function(organs = "Tongue",
                  cell_type_index = 2,
                  droplet_file_index = 24,
                  droplet_organs = "Tongue",
                  data_id = "data15")

########## data16 FigShare Smartseq2 and 10× batch (Thymus)
Figshare_function(organs = "Thymus",
                  cell_type_index = 2,
                  droplet_file_index = 23,
                  droplet_organs = "Thymus",
                  data_id = "data16")

########## data17 FigShare Smartseq2 and 10× batch (Kidney)
Figshare_function(organs = "Kidney",
                  cell_type_index = 2,
                  droplet_file_index = 5:7,
                  droplet_organs = "Kidney",
                  data_id = "data17")

########## data18 FigShare Smartseq2 and 10× batch (Liver)
Figshare_function(organs = "Liver",
                  cell_type_index = 2,
                  droplet_file_index = 10:12,
                  droplet_organs = "Liver",
                  data_id = "data18")

########## data19 FigShare Smartseq2 and 10× batch (Limb_Muscle)
Figshare_function(organs = "Limb_Muscle",
                  cell_type_index = 2,
                  droplet_file_index = 8:9,
                  droplet_organs = "Limb_Muscle",
                  data_id = "data19")

########## data20 FigShare Smartseq2 and 10× batch (Spleen)
Figshare_function(organs = "Spleen",
                  cell_type_index = 2,
                  droplet_file_index = 21:22,
                  droplet_organs = "Spleen",
                  data_id = "data20")



########## data21 GSE54695
data <- read.table(file = './GSE54695/GSE54695_data_transcript_counts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, c(1:80, 161:240)]
## cell group
treatment <- c(rep("2i", 80), rep("serum", 80))
group_condition <- ifelse(treatment == "2i", 2, 1)
## ERCC count
ERCC_count <- data[grep(rownames(data), pattern = "^ERCC-"), ]
## Gene transformation
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "transcript_biotype",
                                            "external_gene_name"),
                             mart = ensembl) %>%
  dplyr::filter("external_gene_name" != "")
gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
data <- data[gene_filter, ]
rownames(data) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
## rownames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## combine
data <- rbind(data, ERCC_count)
## data information
data_info <- simutils::meta_info(id = "data21_GSE54695",
                                 repository = "GEO",
                                 accession_number = "GSE54695",
                                 platform = "CEL-seq",
                                 species = "Mus musculus",
                                 organ = "Embryonic Stem Cells",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 ERCC = TRUE,
                                 dilution_factor = 50000,
                                 volume = 1,
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data21_GSE54695.rds')



########## data22 GSE86469
data <- read.csv(file = './GSE86469/GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz',
                 header = TRUE,
                 row.names = 1) %>% round()
colnames(data) <- str_replace_all(colnames(data), pattern = "[.]", replacement = "_")
colnames(data) <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 1]
control <- paste0("X", c(1:3, 5:9), "th")
## treatment
treatment <- ifelse(colnames(data) %in% control, "control", "case")
## group
group_condition <- ifelse(treatment == "control", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data22_GSE86469",
                                 repository = "GEO",
                                 accession_number = "GSE86469",
                                 platform = "Fluidigm C1",
                                 species = "Homo sapiens",
                                 organ = "Pancreas",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 group_condition = group_condition,
                                 treatment = treatment)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data22_GSE86469.rds')



########## data23 GSE62270
data_list <- list.files(path = './GSE62270/', pattern = 'Whole_Organoid')
data <- read.table(file = paste0('./GSE62270/', data_list[1]),
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
for(i in 2:length(data_list)){
  data_tmp <- read.table(file = paste0('./GSE62270/', data_list[i]),
                         header = TRUE,
                         row.names = 1,
                         sep = '\t')
  data <- cbind(data, data_tmp)
  print(dim(data))
}
data <- round(data)
## Gene replicates removal
data$gene <- str_split(rownames(data), pattern = "__", simplify = TRUE)[, 1]
data <- aggregate(.~data$gene, data[, -673], mean)
rownames(data) <- data$`data$gene`
data <- round(data[, -1])
## ERCC count
ERCC_count <- data[grep(rownames(data), pattern = "^ERCC-"), ]
## Gene transformation
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "transcript_biotype",
                                            "external_gene_name"),
                             mart = ensembl) %>%
  dplyr::filter("external_gene_name" != "")
gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
data <- data[gene_filter, ]
rownames(data) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
## combine
data <- rbind(data, ERCC_count)
## cell batch
batch <- c(rep("batch1", 288), rep("batch2", 384))
## data information
data_info <- simutils::meta_info(id = "data23_GSE62270",
                                 repository = "GEO",
                                 accession_number = "GSE62270",
                                 platform = "CEL-seq",
                                 species = "Mus musculus",
                                 organ = "Embryonic Stem Cells",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 ERCC = TRUE,
                                 dilution_factor = 50000,
                                 volume = 0.03,
                                 batch_info = batch)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data23_GSE62270.rds')



########## data24 GSE85755
data <- read.csv(file = './GSE85755/GSE85755_Merged_CEL-Seq_AllMiceAndLibrariesMerged.csv.gz',
                 header = TRUE,
                 row.names = 1,
                 sep = '\t')
## rownames
gene_name <- str_split(rownames(data), pattern = "__", simplify = TRUE)[, 1]
data$gene_name <- gene_name
data <- aggregate(.~ data$gene_name, data[, -ncol(data)], mean)
rownames(data) <- data$`data$gene_name`
data <- round(data[, -1])
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data24_GSE85755",
                                 repository = "GEO",
                                 accession_number = "GSE85755",
                                 platform = "CEL-seq",
                                 species = "Mus musculus",
                                 organ = "Muscle Stem Cells",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data24_GSE85755.rds')



########## data25 GSE85241
data <- read.csv(file = './GSE85241/GSE85241_cellsystems_dataset_4donors_updated.csv.gz',
                 header = TRUE,
                 row.names = 1,
                 sep = '\t') %>% 
  round()
gene_name <- str_split(rownames(data), pattern = "__", simplify = TRUE)[, 1]
data$gene_name <- gene_name
data <- aggregate(.~ data$gene_name, data[, -ncol(data)], mean)
rownames(data) <- data$`data$gene_name`
data <- round(data[, -1])
dim(data)
## data information
data_info <- simutils::meta_info(id = "data25_GSE85241",
                                 repository = "GEO",
                                 accession_number = "GSE85241",
                                 platform = "CEL-seq2",
                                 species = "Homo sapiens",
                                 organ = "Pancreas",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data25_GSE85241.rds')



########## data26 GSE78779
data1 <- read.table("./GSE78779/GSE78779_CS1_manual.txt.gz",
                    header = TRUE,
                    row.names = 1)
data2 <- read.table("./GSE78779/GSE78779_CS2_manual.txt.gz",
                    header = TRUE,
                    row.names = 1)
data <- cbind(data1, data2)
## ERCC count
ERCC_index <- grep(rownames(data), pattern = "^ERCC-")
ERCC_count <- data[ERCC_index, ]
data <- rbind(data[-ERCC_index, ], ERCC_count)
## batch
batch <- c(rep("CEL-seq", 24),
           rep("CEL-seq2", 20))
## data information
data_info <- simutils::meta_info(id = "data26_GSE78779",
                                 repository = "GEO",
                                 accession_number = "GSE78779",
                                 platform = c("CEL-seq", "CEL-seq2"),
                                 species = "Mus musculus",
                                 organ = "Fibroblast cells",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 ERCC = TRUE,
                                 dilution_factor = 1000000,
                                 volume = 1,
                                 batch_info = batch)
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data26_GSE78779.rds')



########## data27 GSE117617
untar(tarfile = './GSE117617/GSE117617_RAW.tar', exdir = './GSE117617')
data <- read.csv(file = './GSE117617/GSM3305230_CelSeq2_Mixture_Sample_counts.csv.gz',
                 header = TRUE,
                 row.names = 1)
## data information
data_info <- simutils::meta_info(id = "data27_GSE117617",
                                 repository = "GEO",
                                 accession_number = "GSE117617",
                                 platform = "CEL-seq2",
                                 species = "Homo sapiens",
                                 organ = "Lung adenocarcinoma",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data27_GSE117617.rds')



########## data28-data32 GSE117450
untar(tarfile = './GSE117450/GSE117450_RAW.tar', exdir = './GSE117450')
data_list <- list.files(path = './GSE117450/', pattern = 'csv.gz')
for(i in 1:length(data_list)){
  data <- read.table(file = paste0('./GSE117450/', data_list[i]),
                     header = TRUE,
                     row.names = 1,
                     sep = ',')
  ## data information
  data_info <- simutils::meta_info(id = paste0("data", i+27, "_GSE117450_subset", i),
                                   repository = "GEO",
                                   accession_number = "GSE117617",
                                   platform = "CEL-seq2",
                                   species = "Homo sapiens",
                                   organ = "Lung adenocarcinoma",
                                   cell_num = ncol(data),
                                   gene_num = nrow(data))
  print(dim(data))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', "data", i+27, "_GSE117450_subset", i, '.rds'))
}



########## data33-data34 GSE133539
data_list <- list.files(path = './GSE133539/', pattern = 'mat.tsv.gz')
Organism <- c("Homo sapiens", "Mus musculus")
sources <- c("Human PBMC", 
             "Mouse Colon Cell Lines")
for(i in 1:length(data_list)){
  data <- read.table(file = paste0('./GSE133539/',
                                   data_list[i]),
                     header = TRUE,
                     row.names = 1,
                     sep = '\t')
  gene_name <- str_split(rownames(data), pattern = "[.]", simplify = TRUE)[, 1]
  data$gene_name <- gene_name
  data <- aggregate(.~ data$gene_name, data[, -ncol(data)], mean)
  rownames(data) <- data$`data$gene_name`
  data <- round(data[, -1])
  dim(data)
  ## data information
  data_info <- simutils::meta_info(id = paste0("data", i+32, "_GSE133539_subset", i),
                                   repository = "GEO",
                                   accession_number = "GSE133539",
                                   platform = "CEL-seq2",
                                   species = Organism[i],
                                   organ = sources[i],
                                   cell_num = ncol(data),
                                   gene_num = nrow(data))
  print(dim(data))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', "data", i+32, "_GSE133539_subset", i, '.rds'))
}



########## data35-data39 GSE109205
untar(tarfile = './GSE109205/GSE109205_RAW.tar', exdir = './GSE109205')
data_list <- list.files(path = './GSE109205/', pattern = 'txt.gz')
for(i in 1:length(data_list)){
  data <- read.table(file = paste0('./GSE109205/', data_list[i]),
                     header = TRUE,
                     row.names = 1,
                     sep = '\t')
  ## rownames
  rownames(data) <- str_replace_all(rownames(data),
                                    pattern = "_",
                                    replacement = "-")
  ## Filter
  data <- data[, colSums(data) > 0]
  ## data information
  data_info <- simutils::meta_info(id = paste0("data", i+34, "_GSE109205_subset", i),
                                   repository = "GEO",
                                   accession_number = "GSE109205",
                                   platform = "Drop-seq",
                                   species = "Homo sapiens",
                                   organ = "Kidney",
                                   cell_num = ncol(data),
                                   gene_num = nrow(data))
  print(dim(data))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', "data", i+34, "_GSE109205_subset", i, '.rds'))
}



########## data40 GSE65525 subset1
untar(tarfile = './GSE65525/GSE65525_RAW.tar', exdir = './GSE65525')
data1 <- read.csv(file = "./GSE65525/GSM1599495_ES_d0_biorep_techrep1.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data2 <- read.csv(file = "./GSE65525/GSM1599496_ES_d0_biorep_techrep2.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data <- cbind(data1, data2)
## batch
batch <- c(rep("rep1", ncol(data1)),
           rep("rep2", ncol(data2)))
## data information
data_info <- simutils::meta_info(id = "data40_GSE65525_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE65525",
                                 platform = "inDrop",
                                 species = "Mus musculus",
                                 organ = "Embryonic stem cell",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 batch_info = batch)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data40_GSE65525_subset1.rds')



########## data41 GSE65525 subset2
untar(tarfile = './GSE65525/GSE65525_RAW.tar', exdir = './GSE65525')
data1 <- read.csv(file = "./GSE65525/GSM1599497_ES_d2_LIFminus.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data2 <- read.csv(file = "./GSE65525/GSM1599499_ES_d7_LIFminus.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data <- cbind(data1, data2)
## treatment
treatment <- c(rep("d2", ncol(data1)),
               rep("d7", ncol(data2)))
group_condition <- ifelse(treatment == "d2", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data41_GSE65525_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE65525",
                                 platform = "inDrop",
                                 species = "Mus musculus",
                                 organ = "Embryonic stem cell",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data41_GSE65525_subset2.rds')



########## data42 GSE65525 subset3
untar(tarfile = './GSE65525/GSE65525_RAW.tar', exdir = './GSE65525')
data1 <- read.csv(file = "./GSE65525/GSM1599500_K562_cells.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data2 <- read.csv(file = "./GSE65525/GSM1599501_K562_pure_RNA.csv.bz2",
                  header = FALSE,
                  row.names = 1)
data <- cbind(data1, data2)
## treatment
treatment <- c(rep("case", ncol(data1)),
               rep("control", ncol(data2)))
group_condition <- ifelse(treatment == "case", 2, 1)
## ERCC count
ERCC_count <- data[grep(rownames(data), pattern = "^ERCC-"), ]
## Gene transformation
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
id_convert <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                            "transcript_biotype",
                                            "external_gene_name"),
                             mart = ensembl) %>%
  dplyr::filter("external_gene_name" != "")
gene_filter <- id_convert$external_gene_name[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
data <- data[gene_filter, ]
rownames(data) <- id_convert$ensembl_gene_id[stats::na.omit(match(rownames(data),
                                                                  id_convert$external_gene_name))]
data <- rbind(data, ERCC_count)
## data information
data_info <- simutils::meta_info(id = "data42_GSE65525_subset3",
                                 repository = "GEO",
                                 accession_number = "GSE65525",
                                 platform = "inDrop",
                                 species = "Homo sapiens",
                                 organ = "Human lymphoblastoma",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 ERCC = TRUE,
                                 dilution_factor = 5000,
                                 volume = 0.001,
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data42_GSE65525_subset3.rds')



########## data43 GSE183590
data <- read.table(file = './GSE183590/GSE183590_Read_count_matrix_05092021.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
## data information
data_info <- simutils::meta_info(id = "data43_GSE183590",
                                 repository = "GEO",
                                 accession_number = "GSE183590",
                                 platform = "Fluidigm C1",
                                 species = "Homo sapiens",
                                 organ = "Lung",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data43_GSE183590.rds')



########## data44 GSE67835
untar(tarfile = './GSE67835/GSE67835_RAW.tar', exdir = './GSE67835')
data_list <- list.files(path = './GSE67835', pattern = "csv.gz")
data <- merged_fun(dir = './GSE67835', list = data_list)
colnames(data) <- paste0("Cell", 1:ncol(data))
## rownames
rownames(data) <- gsub(pattern = " ", replacement = "", rownames(data))
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data44_GSE67835",
                                 repository = "GEO",
                                 accession_number = "GSE67835",
                                 platform = "Fluidigm C1",
                                 species = "Homo sapiens",
                                 organ = "Brain",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data44_GSE67835.rds')



########## data45 GSE60361
data <- read.table(file = './GSE60361/GSE60361_C1-3005-Expression.txt.gz',
                   header = TRUE,
                   sep = '\t')
data$cell_id[which(data$cell_id %in% c("Mar-01", "Mar-02"))] <- c("March2", "March1", "Marchf2", "Marchf1")
rownames(data) <- data$cell_id
data <- data[, -1]
## data information
data_info <- simutils::meta_info(id = "data45_GSE60361",
                                 repository = "GEO",
                                 accession_number = "GSE60361",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Cerebral cortex/Hippocampus",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data))
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data45_GSE60361.rds')



########## data46 GSE95432
data <- read.table(file = './GSE95432/GSE95432_GeneCounts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, -1]
## cluster
cluster <- c(rep("basal cell", 96), rep("luminal cell", 90))
## data information
data_info <- simutils::meta_info(id = "data46_GSE95432",
                                 repository = "GEO",
                                 accession_number = "GSE95432",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Mammary Gland",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 cluster_info = cluster)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data46_GSE95432.rds')



########## data47 GSE95434
data <- read.table(file = './GSE95434/GSE95434_GeneCounts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, -1]
## treatment
treatment <- c(rep("week2", 144),
               rep("week4", 136),
               rep("week10", 66))
group_condition <- ifelse(treatment == "week2", 1, ifelse(treatment == "week4", 2, 3))
## data information
data_info <- simutils::meta_info(id = "data47_GSE95434",
                                 repository = "GEO",
                                 accession_number = "GSE95434",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Mammary Gland",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data47_GSE95434.rds')



########## data48 GSE95436
data <- read.table(file = './GSE95436/GSE95436_GeneCounts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, -1]
## cluster
cluster <- c(rep("luminal cell", 43),
             rep("basal cell", 112),
             rep("epithelial cell", 123))
## data information
data_info <- simutils::meta_info(id = "data48_GSE95436",
                                 repository = "GEO",
                                 accession_number = "GSE95436",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Mammary Gland",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 cluster_info = cluster)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data48_GSE95436.rds')



########## data49 GSE95445
data <- read.table(file = './GSE95445/GSE95445_GeneCounts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, -1]
## cluster
cluster <- c(rep("puberty epithelial cell", 221),
             rep("adult epithelial cell", 223))
## data information
data_info <- simutils::meta_info(id = "data49_GSE95445",
                                 repository = "GEO",
                                 accession_number = "GSE95445",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Mammary Gland",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 cluster_info = cluster)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data49_GSE95445.rds')



########## data50 GSE95448
data <- read.table(file = './GSE95448/GSE95448_GeneCounts.txt.gz',
                   header = TRUE,
                   row.names = 1,
                   sep = '\t')
data <- data[, -1]
## cluster
cluster <- c(rep("Adult basal cell", 237),
             rep("Pregnancy basal cell", 75))
## data information
data_info <- simutils::meta_info(id = "data50_GSE95448",
                                 repository = "GEO",
                                 accession_number = "GSE95448",
                                 platform = "Fluidigm C1",
                                 species = "Mus musculus",
                                 organ = "Mammary Gland",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 cluster_info = cluster)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data50_GSE95448.rds')



########## data51-data60 GSE108097
untar('./GSE108097/GSE108097_RAW.tar', exdir = './GSE108097')
data_list <- list.files(path = './GSE108097/', pattern = "txt.gz")[c(1,2,4,6,9)]
organs <- c("Bladder",
            "Prostate",
            "Mammary Gland",
            "Bone Marrow",
            "Bone Marrow")
for(i in 1:5){
  data <- read.table(file = paste0('./GSE108097/', data_list[i]),
                     header = TRUE,
                     sep = ' ',
                     row.names = 1)
  ## rownames and colnames
  rownames(data) <- str_replace_all(rownames(data),
                                    pattern = "_",
                                    replacement = "-")
  ## data information
  data_info <- simutils::meta_info(id = paste0("data", i+50, "_GSE108097_subset", i),
                                   repository = "GEO",
                                   accession_number = "GSE108097",
                                   platform = "Microwell-seq",
                                   species = "Mus musculus",
                                   organ = organs[i],
                                   cell_num = ncol(data),
                                   gene_num = nrow(data))
  print(dim(data))
  ## Save
  data <- list(data = as.matrix(data),
               data_info = data_info)
  saveRDS(data, file = paste0('../preprocessed_data/', "data", i+50, "_GSE108097_subset", i, '.rds'))
}



########## data56 GSE167297 subset1
untar("./GSE167297/GSE167297_RAW.tar", exdir = "./GSE167297")
data1 <- read.table("./GSE167297/GSM5101013_Pt1_Normal_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data2 <- read.table("./GSE167297/GSM5101015_Pt1_Deep_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data <- cbind(data1, data2)
## treatment
treatment <- c(rep("control", ncol(data1)),
               rep("case", ncol(data2)))
## group
group_condition <- ifelse(treatment == "control", 1, 2)
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data56_GSE167297_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE167297",
                                 platform = "10X Genomics",
                                 species = "Homo sapiens",
                                 organ = "Gastric Cancer",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data56_GSE167297_subset1.rds")



########## data57 GSE167297 subset2
data1 <- read.table("./GSE167297/GSM5101018_Pt3_Normal_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data2 <- read.table("./GSE167297/GSM5101020_Pt3_Deep_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data <- cbind(data1, data2)
## treatment
treatment <- c(rep("control", ncol(data1)),
               rep("case", ncol(data2)))
## group
group_condition <- ifelse(treatment == "control", 1, 2)
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data57_GSE167297_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE167297",
                                 platform = "10X Genomics",
                                 species = "Homo sapiens",
                                 organ = "Gastric Cancer",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data57_GSE167297_subset2.rds")



########## data58 GSE167297 subset3
data1 <- read.table("./GSE167297/GSM5101021_Pt4_Normal_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data2 <- read.table("./GSE167297/GSM5101023_Pt4_Deep_CountMatrix.txt.gz",
                    header = TRUE,
                    row.names = 1)
data <- cbind(data1, data2)
## treatment
treatment <- c(rep("control", ncol(data1)),
               rep("case", ncol(data2)))
## group
group_condition <- ifelse(treatment == "control", 1, 2)
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data58_GSE167297_subset3",
                                 repository = "GEO",
                                 accession_number = "GSE167297",
                                 platform = "10X Genomics",
                                 species = "Homo sapiens",
                                 organ = "Gastric Cancer",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data58_GSE167297_subset3.rds")



########## data59 GSE184950
dir.create("./GSE184950")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5602315&format=file&file=GSM5602315%5FA10%2Etar%2Egz",
              destfile = "./GSE184950/GSM5602315_A10.tar.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5602322&format=file&file=GSM5602322%5FB13%2Etar%2Egz",
              destfile = "./GSE184950/GSM5602322_B13.tar.gz")

gunzip(filename = "./GSE184950/GSM5602315_A10.tar.gz")
untar("./GSE184950/GSM5602315_A10.tar", exdir = "./GSE184950/control")
gunzip(filename = "./GSE184950/GSM5602322_B13.tar.gz")
untar("./GSE184950/GSM5602322_B13.tar", exdir = "./GSE184950/case")

data1 <- as.data.frame(Read10X("./GSE184950/control/filtered_feature_bc_matrix/"))
data2 <- as.data.frame(Read10X("./GSE184950/case//filtered_feature_bc_matrix/"))
data <- cbind(data1, data2)

## treatment
treatment <- c(rep("control", ncol(data1)),
               rep("case", ncol(data2)))
## group
group_condition <- ifelse(treatment == "control", 1, 2)
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data59_GSE184950",
                                 repository = "GEO",
                                 accession_number = "GSE184950",
                                 platform = "10X Genomics",
                                 species = "Homo sapiens",
                                 organ = "Brain",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data59_GSE184950.rds")



########## data60 SCP1261
data <- read.csv("/Users/duohongrui/Downloads/zuizhong_raw_data.csv",
                 header = TRUE,
                 row.names = 1)
cell_info <- read.csv("/Users/duohongrui/Downloads/zuizhong_metada.csv",
                      header = TRUE,
                      row.names = 1)
cell_info <- cell_info[-1, ]
## treatment
treatment <- cell_info$donor_id
## group
group_condition <- ifelse(treatment == "Old_mice", 1, 2)
## cluster
cluster <- cell_info$cell_type_label
## rownames and colnames
rownames(data) <- str_replace_all(rownames(data),
                                  pattern = "_",
                                  replacement = "-")
## data information
data_info <- simutils::meta_info(id = "data60_SCP1261",
                                 repository = "GEO",
                                 accession_number = "GSE108097",
                                 platform = "10X Genomics",
                                 species = "Mus musculus",
                                 organ = "Tumor",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition,
                                 cluster_info = cluster)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = "../preprocessed_data/data60_SCP1261.rds")



## data61 GSE114727 subset1
untar("./GSE114727/GSE114727_RAW.tar", exdir = "./GSE114727/")
data <- read.csv(file = "./GSE114727/GSM3148587_BC01_NORMAL1_counts.csv.gz",
                 header = TRUE)
data[is.na(data)] <- 0
data <- data[, -1] %>% t() %>% round()
## Filter
index <- colSums(data) > 0
data <- data[, index]
colnames(data) <- paste0("control_cell", 1:ncol(data))
### case
data_case <- read.csv(file = "./GSE114727/GSM3148591_BC01_TUMOR1_counts.csv.gz",
                      header = TRUE)
data_case[is.na(data_case)] <- 0
data_case <- data_case[, -1] %>% t() %>% round()
## Filter
index <- colSums(data_case) > 0
data_case <- data_case[, index]
colnames(data_case) <- paste0("case_cell", 1:ncol(data_case))
### combine
gene <- intersect(rownames(data), rownames(data_case))
data <- cbind(data[gene, ], data_case[gene, ])
### treatment
treatment <- c(rep("normal", 2993),
               rep("tumor", 1589))
### group
group_condition <- ifelse(treatment == "normal", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data61_GSE114727_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE114727",
                                 platform = "inDrop",
                                 species = "Homo sapiens",
                                 organ = "Breast Tumor",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data61_GSE114727_subset1.rds')



## data62 GSE114727 subset2
data <- read.csv(file = "./GSE114727/GSM3148608_BC03_NORMAL1_counts.csv.gz",
                 header = TRUE)
data[is.na(data)] <- 0
data <- data[, -1] %>% t() %>% round()
## Filter
index <- colSums(data) > 0
data <- data[, index]
colnames(data) <- paste0("control_cell", 1:ncol(data))
### case
data_case <- read.csv(file = "./GSE114727/GSM3148611_BC03_TUMOR1_counts.csv.gz",
                      header = TRUE)
data_case[is.na(data_case)] <- 0
data_case <- data_case[, -1] %>% t() %>% round()
## Filter
index <- colSums(data_case) > 0
data_case <- data_case[, index]
colnames(data_case) <- paste0("case_cell", 1:ncol(data_case))
### combine
gene <- intersect(rownames(data), rownames(data_case))
data <- cbind(data[gene, ], data_case[gene, ])
### treatment
treatment <- c(rep("normal", 2480),
               rep("tumor", 402))
### group
group_condition <- ifelse(treatment == "normal", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data62_GSE114727_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE114727",
                                 platform = "inDrop",
                                 species = "Homo sapiens",
                                 organ = "Breast Tumor",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data62_GSE114727_subset2.rds')



## data63 GSE106202
data <- read.table("./GSE106202/GSE106202_FRUCTOSE_GLUCOSE_merge.DGE.txt.gz",
                   header = TRUE,
                   row.names = 1)
### treatment
treatment <- c(rep("Fructose", 799),
               rep("Glucose", 800))
### group
group_condition <- ifelse(treatment == "Fructose", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data63_GSE106202",
                                 repository = "GEO",
                                 accession_number = "GSE106202",
                                 platform = "Drop-seq",
                                 species = "Homo sapiens",
                                 organ = "MDA-MB-231 cells",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data63_GSE106202.rds')



## data64 GSE102827
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102827/matrix/GSE102827_series_matrix.txt.gz",
              destfile = "./GSE102827/GSE102827_series_matrix.txt.gz")
sample_info <- getGEO(filename = "./GSE102827/GSE102827_series_matrix.txt.gz", getGPL = FALSE)
cell_info <- pData(sample_info)
### data
untar("./GSE102827/GSE102827_RAW.tar", exdir = "./GSE102827/untar_file")
data1 <- read.table("./GSE102827/untar_file/GSM2746895_B1_1_0h_counts.txt.gz",
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")
data1 <- data1[, -c(1:3)]
data2 <- read.table("./GSE102827/untar_file/GSM2746896_B1_2_1h_counts.txt.gz",
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")
data2 <- data2[, -c(1:3)]
### combine
data <- cbind(data1, data2)
### treatment
treatment <- c(rep("control", ncol(data1)),
               rep("case", ncol(data2)))
### group
group_condition <- ifelse(treatment == "control", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data64_GSE102827",
                                 repository = "GEO",
                                 accession_number = "GSE102827",
                                 platform = "inDrop",
                                 species = "Mus musculus",
                                 organ = "Visual Cortex",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data64_GSE102827.rds')



## data65 GSE87038 subset1 (forebrain)
GSE87038 <- openxlsx::read.xlsx("./GSE87038/GSE87038_Mouse_Organogenesis_UMI_counts_matrix.xlsx", rowNames = TRUE)
forebrain_index <- str_starts(colnames(GSE87038), pattern = "^forebrain")
data <- GSE87038[, forebrain_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 96),
               rep("E11.5", 64))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data65_GSE87038_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "forebrain",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data65_GSE87038_subset1.rds')



## data66 GSE87038 subset2 (heart)
heart_index <- str_starts(colnames(GSE87038), pattern = "^heart")
data <- GSE87038[, heart_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 96),
               rep("E11.5", 96))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data66_GSE87038_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "heart",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data66_GSE87038_subset2.rds')


## data67 GSE87038 subset3 (somite)
somite_index <- str_starts(colnames(GSE87038), pattern = "^somite")
data <- GSE87038[, somite_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 92),
               rep("E11.5", 64))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data67_GSE87038_subset3",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "somite",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data67_GSE87038_subset3.rds')



## data68 GSE87038 subset4 (intestine)
intestine_index <- str_starts(colnames(GSE87038), pattern = "^intestine")
data <- GSE87038[, intestine_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 89),
               rep("E11.5", 64))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data68_GSE87038_subset4",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "intestine",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data68_GSE87038_subset4.rds')



## data69 GSE87038 subset5 (hindbrain)
hindbrain_index <- str_starts(colnames(GSE87038), pattern = "^hindbrain")
data <- GSE87038[, hindbrain_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 87),
               rep("E11.5", 64))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data69_GSE87038_subset5",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "hindbrain",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data69_GSE87038_subset5.rds')



## data70 GSE87038 subset6 (skin)
skin_index <- str_starts(colnames(GSE87038), pattern = "^skin")
data <- GSE87038[, skin_index]
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 2] %in% c("E9.5", "E11.5")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E9.5", 68),
               rep("E11.5", 64))
### group
group_condition <- ifelse(treatment == "E9.5", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data70_GSE87038_subset6",
                                 repository = "GEO",
                                 accession_number = "GSE87038",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "skin",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data70_GSE87038_subset6.rds')



## data71 GSE90047 subset1
data <- read.table("./GSE90047/GSE90047_Single-cell_RNA-seq_Read_Count.txt.gz",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t")
data <- data[, -c(1,2)]
colnames(data) <- str_replace_all(colnames(data), pattern = "[.]", replacement = "_")
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 1] %in% c("E10", "E17")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E10", 54),
               rep("E17", 70))
### group
group_condition <- ifelse(treatment == "E10", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data71_GSE90047_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE90047",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "liver",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data71_GSE90047_subset1.rds')



## data72 GSE90047 subset2
data <- read.table("./GSE90047/GSE90047_Single-cell_RNA-seq_Read_Count.txt.gz",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t")
data <- data[, -c(1,2)]
colnames(data) <- str_replace_all(colnames(data), pattern = "[.]", replacement = "_")
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 1] %in% c("E11", "E15")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E11", 70),
               rep("E15", 77))
### group
group_condition <- ifelse(treatment == "E11", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data72_GSE90047_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE90047",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "liver",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data72_GSE90047_subset2.rds')



## data73 GSE90047 subset3
data <- read.table("./GSE90047/GSE90047_Single-cell_RNA-seq_Read_Count.txt.gz",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t")
data <- data[, -c(1,2)]
colnames(data) <- str_replace_all(colnames(data), pattern = "[.]", replacement = "_")
filter_index <- str_split(colnames(data), pattern = "_", simplify = TRUE)[, 1] %in% c("E12", "E13")
data <- data[, filter_index]
### treatment
treatment <- c(rep("E12", 41),
               rep("E13", 65))
### group
group_condition <- ifelse(treatment == "E12", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data73_GSE90047_subset3",
                                 repository = "GEO",
                                 accession_number = "GSE90047",
                                 platform = "Smart-seq2",
                                 species = "Mus musculus",
                                 organ = "liver",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data73_GSE90047_subset3.rds')



## data74 GSE112004 subset1
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112004/matrix/GSE112004_series_matrix.txt.gz",
              destfile = "./GSE112004/GSE112004_series_matrix.txt.gz")
sample_info <- getGEO(filename = "./GSE112004/GSE112004_series_matrix.txt.gz", getGPL = FALSE)
cell_info <- pData(sample_info)
data_list <- list.files("./GSE112004/", pattern = "tsv.gz")
GSE112004 <- read.table(file.path("./GSE112004", data_list[1]),
                        header = TRUE,
                        row.names = 1)

gene_list <- purrr::map(data_list, function(x){
  data_tmp <- read.table(file.path("./GSE112004", x),
                         header = TRUE,
                         row.names = 1)
  rownames(data_tmp)
})
gene <- intersect(gene_list[[1]], gene_list[[2]])
for(i in 3:20){
  gene <- intersect(gene, gene_list[[i]])
  print(i)
}
## combine
for(i in data_list[2:20]){
  data_tmp <- read.table(file.path("./GSE112004", i),
                         header = TRUE,
                         row.names = 1)
  print(dim(data_tmp))
  print(i)
  GSE112004 <- cbind(GSE112004[gene, ], data_tmp[gene, ])
}
colnames(GSE112004) <- str_remove_all(colnames(GSE112004), pattern = "[X]")
colnames(GSE112004) <- str_replace_all(colnames(GSE112004), pattern = "[.]", replacement = "_")
## rownames
gene_name <- str_split(rownames(GSE112004), pattern = "[.]", simplify = TRUE)[, 1]
GSE112004$gene_name <- gene_name
GSE112004 <- aggregate(.~GSE112004$gene_name, GSE112004[, -ncol(GSE112004)], mean)
rownames(GSE112004) <- GSE112004$`GSE112004$gene_name`
GSE112004 <- round(GSE112004[, -1])
## time 0h vs 6h
index <- cell_info$characteristics_ch1.4 %in% c("time point: 0h", "time point: 6h")
data <- GSE112004[, index]
### treatment
treatment <- c(rep("0h", 384),
               rep("6h", 384))
### group
group_condition <- ifelse(treatment == "0h", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data74_GSE112004_subset1",
                                 repository = "GEO",
                                 accession_number = "GSE112004",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data74_GSE112004_subset1.rds')



## data75 GSE112004 subset2
## time 18h vs 42h
index <- cell_info$characteristics_ch1.4 %in% c("time point: 18h", "time point: 42h")
data <- GSE112004[, index]
### treatment
treatment <- c(rep("18h", 384),
               rep("42h", 384))
### group
group_condition <- ifelse(treatment == "18h", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data75_GSE112004_subset2",
                                 repository = "GEO",
                                 accession_number = "GSE112004",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data75_GSE112004_subset2.rds')



## data76 GSE112004 subset3
## time 66h vs 114h
index <- cell_info$characteristics_ch1.4 %in% c("time point: 66h", "time point: 114h")
data <- GSE112004[, index]
### treatment
treatment <- c(rep("66h", 384),
               rep("114h", 384))
### group
group_condition <- ifelse(treatment == "66h", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data76_GSE112004_subset3",
                                 repository = "GEO",
                                 accession_number = "GSE112004",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data76_GSE112004_subset3.rds')



## data77 GSE112004 subset4
## time day 2 vs day 4
index <- cell_info$characteristics_ch1.4 %in% c("time point: day 2", "time point: day 4")
data <- GSE112004[, index]
### treatment
treatment <- c(rep("2d", 384),
               rep("4d", 384))
### group
group_condition <- ifelse(treatment == "2d", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data77_GSE112004_subset4",
                                 repository = "GEO",
                                 accession_number = "GSE112004",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data77_GSE112004_subset4.rds')




## data78 GSE112004 subset4
## time day 6 vs day 8
index <- cell_info$characteristics_ch1.4 %in% c("time point: day 6", "time point: day 8")
data <- GSE112004[, index]
### treatment
treatment <- c(rep("6d", 384),
               rep("8d", 384))
### group
group_condition <- ifelse(treatment == "6d", 1, 2)
## data information
data_info <- simutils::meta_info(id = "data78_GSE112004_subset5",
                                 repository = "GEO",
                                 accession_number = "GSE112004",
                                 platform = "MARS-Seq",
                                 species = "Mus musculus",
                                 organ = "Bone Marrow",
                                 cell_num = ncol(data),
                                 gene_num = nrow(data),
                                 treatment = treatment,
                                 group_condition = group_condition)
print(dim(data))
## Save
data <- list(data = as.matrix(data),
             data_info = data_info)
saveRDS(data, file = '../preprocessed_data/data78_GSE112004_subset5.rds')



########## data79-data101
data_list <- list.files(path = './zenodo_datasets/')
plaf <- c("Smart-seq",
          "Smart-seq",
          "Smart-seq",
          "Smart-seq",
          "Fluidigm C1",
          "CEL-seq2",
          "CEL-seq2",
          "CEL-seq2",
          "CEL-seq2",
          "Fluidigm C1",
          "Smart-seq2",
          "Smart-seq2",
          "Fluidigm C1",
          "Fluidigm C1",
          "Smart-seq2",
          "Fluidigm C1",
          "Fluidigm C1",
          "Fluidigm C1",
          "Smart-seq2",
          "Smart-seq2",
          "Smart-seq2",
          "Smart-seq2",
          "Smart-seq")
org <- c("Mus musculus",
         "Mus musculus",
         "Mus musculus",
         "Mus musculus",
         "Mus musculus",
         "Homo sapiens",
         "Homo sapiens",
         "Homo sapiens",
         "Homo sapiens",
         "Mus musculus",
         "Homo sapiens",
         "Homo sapiens",
         "Mus musculus",
         "Homo sapiens",
         "Mus musculus",
         "Homo sapiens",
         "Homo sapiens",
         "Mus musculus",
         "Mus musculus",
         "Mus musculus",
         "Homo sapiens",
         "Homo sapiens",
         "Mus musculus")
sources <- c("Cell lines",
             "Cell lines",
             "Hematopoietic stem and progenitor cells",
             "Hematopoietic stem and progenitor cells",
             "mESC",
             "Adenocarcinoma cell lines",
             "Adenocarcinoma cell lines",
             "Adenocarcinoma cell lines",
             "Adenocarcinoma cell lines",
             "Bone Marrow",
             "Fetal germ cells",
             "Fetal germ cells",
             "Bone Marrow",
             "Embryos",
             "Bone Marrow",
             "Embryonic stem cells",
             "Human Skeletal Muscle Myoblasts",
             "Cell lines",
             "Pancreas",
             "Pancreas",
             "Human cortical spheroid-derived cells",
             "Human cortical spheroid-derived cells",
             "Cell lines"
)


for(i in 1:23){
  data <- readRDS(paste0('./zenodo_datasets/', data_list[i]))
  if(i %in% c(19, 20)){
    ERCC <- TRUE
    dilution_factor <- 500000
    volume <- 0.05
  }else{
    ERCC <- FALSE
    dilution_factor <- NULL
    volume <- NULL
  }
  data_info <- simutils::meta_info(id = data_list[i],
                                   repository = "zenodo",
                                   URL = "https://zenodo.org/record/1443566",
                                   platform = plaf[i],
                                   species = org[i],
                                   organ = sources[i],
                                   cell_num = dim(data[["counts"]])[1],
                                   gene_num = dim(data[["counts"]])[2],
                                   data_type = "count",
                                   ERCC = ERCC,
                                   dilution_factor = dilution_factor,
                                   volume = volume,
                                   group_condition = as.numeric(as.factor(data[["grouping"]])),
                                   cluster_info = data[["grouping"]])
  data <- list(data = data,
               data_info = data_info)
  print(data_info$id)
  print(data_info$platform)
  print(data_info$species)
  print(data_info$organ)
  print(data_info$cell_num)
  print(data_info$gene_num)
  print(data_info$dilution_factor)
  print(data_info$volume)
  print(length(unique(data_info$group_condition)))
  print(table(rowSums(data[["data"]][["counts"]]) > 0))
  saveRDS(data, file = paste0('../preprocessed_data/', data_list[i]))
}


################################################################################

#                            Spatial Datasets

################################################################################

########## data102





