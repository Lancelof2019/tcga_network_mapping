options(stringsAsFactors = F)
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19")
cancer_type<-"SKCM"
#exp_intgr <- readRDS("./data/tcga_data_processed/SKCM_exp_intgr.RData")
#mty_intgr <- readRDS("./data/tcga_data_processed/SKCM_mty_intgr.RData")

exp_intgr<- readRDS(paste0("./data/",cancer_type,"_exp_assy_all_test.RData"))
mty_intgr<- readRDS(paste0("./data/",cancer_type,"_mty_mat_all_test.RData"))
#"./data/mty_mat_all.RData")
dim(exp_intgr)
dim(mty_intgr)
#saveRDS(exp_assy_trial,"./data/exp_assy_all.RData")


hfeat <- colnames(exp_intgr) # The features of gene expression data
mfeat <- colnames(mty_intgr) # The features of DNA methylation data
#View(hfeat)
#print(hfeat)
selected_features <- matrix(ncol = 3, nrow=30000) # Matrix of 3 columns; column1: community, column2: genomic type, column3: mapped component

#print(getwd())
melanet_cmt<-readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

#length(melanet_cmt)
#print(melanet_cmt)

len <- 0
j <- 1
indx <- NULL
for(i in 1:length(melanet_cmt)){
  cmt = melanet_cmt[[i]]
  # Start mapping
  # Gene expression
  j = j + length(indx)
  indx = NULL
  indx = which(hfeat %in% cmt)
  if (length(indx) != 0){
    len = len +length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 1
    selected_features[j:len,3] = indx
  }
  
  # DNA methylation
  j = j+length(indx)
  indx = NULL
  indx = which(mfeat %in% cmt)
  if (length(indx) != 0){
    len = len + length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 2
    selected_features[j:len,3] = indx
  }
}

selected_features <- matrix(ncol = 3, nrow=30000) # Matrix of 3 columns; column1: community, column2: genomic type, column3: mapped component
len <- 0
j <- 1
indx <- NULL
for(i in 1:length(melanet_cmt)){
  cmt = melanet_cmt[[i]]
  # Start mapping
  # Gene expression
  j = j + length(indx)
  indx = NULL
  indx = which(hfeat %in% cmt)
  if (length(indx) != 0){
    len = len +length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 1
    selected_features[j:len,3] = indx
  }
  
  # DNA methylation
  j = j+length(indx)
  indx = NULL
  indx = which(mfeat %in% cmt)
  if (length(indx) != 0){
    len = len + length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 2
    selected_features[j:len,3] = indx
  }
}
selected_features <- na.omit(selected_features)
#write.csv(selected_features, "./data/python_related/data/selected_features03.csv", row.names = F)





#selected_features <- na.omit(selected_features)
write.csv(selected_features, paste0("./data/python_related/data/",cancer_type,"_selected_features_test03.csv"), row.names = F)
write.csv(colnames(exp_intgr), paste0("./data/python_related/data/",cancer_type,"_exp_feature_names_test03.csv"), row.names = F)
write.csv(colnames(mty_intgr), paste0("./data/python_related/data/",cancer_type,"_mty_feature_names_test03.csv"), row.names = F) 
write.csv(exp_intgr, paste0("./data/python_related/data/",cancer_type,"_exp_intgr_test03.csv"))
write.csv(mty_intgr, paste0("./data/python_related/data/",cancer_type,"_mty_intgr_test03.csv"))
