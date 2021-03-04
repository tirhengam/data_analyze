#step1: load libraries
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("limma")
install.packages("devtools")
devtools::install_github("bartongroup/proteusLabelFree")

data(proteusLabelFree)
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
BiocManager::install("DEP")
library(proteusLabelFree)
#library("proteus")
library("SummarizedExperiment")
library("DEP")
library("dplyr")



#step2: load data
data = read.delim("proteinGroups_2.txt", stringsAsFactors = FALSE, as.is=TRUE)

#Explore data to find information such as sample size, number of proteins..etc
str(data)
dim(data)
glimpse(data)
colnames(data)
View(data)
#Find the number of samples based on the names of LFQ intensities
grep("^LFQ.intensity", names(data), value = TRUE)
#Step3: Quality control checks on data
#filter for contaminant proteins and decoy database hits, which are indicated by "+"
#in the columns and one hit wonders

data <- data %>% filter(Potential.contaminant != "+") %>% filter(Reverse!= "+") %>% filter(Only.identified.by.site != "+") <br>

#Check false dicsovery rate (here Q.value) cut offs used to control type 1 error  
  
summary(as.numeric(data$Q.value)) 

#Remove rows with duplicate names
data$Protein.IDs %>% duplicated() %>% any()

data %>% group_by(Protein.IDs) %>% summarize(frequency = n()) %>%  arrange(desc(frequency)) %>% filter(frequency > 1)


#Make unique names using the annotation in the "Protein Id" column as primary names and the
#annotation in "Peptides" as name for those that do not have an Protein name.
data_unique <- make_unique(data, "Protein.IDs", "Peptides", delim = ";")

#Are there any duplicated names?
data$name %>% duplicated() %>% any()

#Generate a SummarizedExperiment object by parsing condition information from the
#column names
LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_se<- make_se_parse(data_unique, LFQ_columns)
data_se
#Check replicate quality of data
pairs(assay(data_se[,1:24]))
#Check standard deviation of values in datset with respect to mean
meanSdPlot(data_norm)

#Filter proteins containing missing values
plot_frequency(data_se)
#Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)
plot_frequency(data_filt)
#Plot a barplot showing the number of identified proteins per samples
plot_numbers(data_filt)


#Step4: Standardisation
#Check if data have skewed distribution
hist(2^assay(data_se)[,"AMP__1"],n=1000)
hist(assay(data_se)[,"AMP__1"],n=1000)

hist(2^assay(data_se)[,"AMP__2"],n=1000)
hist(assay(data_se)[,"AMP__2"],n=1000)

hist(2^assay(data_se)[,"AMP__3"],n=1000)
hist(assay(data_se)[,"AMP__3"],n=1000)

plot(density(assay(data_se)[,"AMP__1"],na.rm=TRUE))
plot(density(assay(data_se)[,"AMP__2"],na.rm=TRUE))
plot(density(assay(data_se)[,"AMP__3"],na.rm=TRUE))

#Step5 Normalization
#data_norm <- normalize_vsn(data_filt)
#data_norm_med <- data_filt
#assay(data_norm_med)<-proteus::normalizeMedian(assay(data_filt))
data_norm_quant <- data_filt
assay(data_norm_quant)<-limma::normalizeQuantiles(assay(data_filt))
plot_normalization(data_filt, data_norm)

#Step6 Imputation of missing values

#check missing values in dataset
plot_missval(data_filt)
plot_detect(data_filt)
#impute(data_norm, fun = "QRILC")
#impute(data_norm, fun = "QRILC")
impute(data_norm, fun = "zero")
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)
plot_imputation(data_norm, data_imp)


#STEP7: Differential expression analysis
#for AMP
data_diff_AMP <- test_diff(data_imp, type = "control", control = "K_") 
#Denote significant proteins based on user defined cutoffs
dep_AMP <- add_rejections(data_diff_AMP, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_AMP, contrast = "AMP__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_AMP, proteins = c("Q88ES8", "Q88QK8" , "Q88L36"))

plot_single(dep_AMP, proteins = c("Q88ES8", "Q88QK8" , "Q88L36"), type = "centered")

plot_cond(dep_AMP)

#STEP8:Visualisations of data
plot_pca(dep_AMP, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_AMP, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_AMP, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_AMP, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_AMP, dep_AMP, file = "data_AMP.RData")

#############################################################
#for CHL
data_diff_CHL <- test_diff(data_imp, type = "control", control = "CHL_") 

#Denote significant proteins based on user defined cutoffs
dep_CHL <- add_rejections(data_diff_CHL, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_CHL, contrast = "CHL__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_CHL, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_CHL, proteins = "Q88ES8", type = "centered")

plot_cond(dep_CHL)

#STEP7:Visualisations of data
plot_pca(dep_CHL, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_CHL, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_CHL, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_CHL, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_CHL, dep_CHL, file = "data_CHL.RData")

#############################################################
#for CIP
data_diff_CIP <- test_diff(data_imp, type = "control", control = "CIP_") 

#Denote significant proteins based on user defined cutoffs
dep_CIP <- add_rejections(data_diff_CIP, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_CIP, contrast = "CIP__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_CIP, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_CIP, proteins = "Q88ES8", type = "centered")

plot_cond(dep_CIP)

#STEP7:Visualisations of data
plot_pca(dep_CIP, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_CIP, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_CIP, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_CIP, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_CIP, dep_CIP, file = "data_CIP.RData")

#############################################################
#for ERY
data_diff_ERY <- test_diff(data_imp, type = "control", control = "ERY_") 

#Denote significant proteins based on user defined cutoffs
dep_ERY <- add_rejections(data_diff_ERY, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_ERY, contrast = "ERY__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_ERY, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_ERY, proteins = "Q88ES8", type = "centered")

plot_cond(dep_ERY)

#STEP7:Visualisations of data
plot_pca(dep_ERY, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_ERY, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_ERY, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_ERY, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_ERY, dep_ERY, file = "data_ERY.RData")


#############################################################
#for KAN
data_diff_KAN <- test_diff(data_imp, type = "control", control = "KAN_") 

#Denote significant proteins based on user defined cutoffs
dep_KAN <- add_rejections(data_diff_KAN, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_KAN, contrast = "KAN__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_KAN, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_KAN, proteins = "Q88ES8", type = "centered")

plot_cond(dep_KAN)

#STEP7:Visualisations of data
plot_pca(dep_KAN, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_KAN, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_KAN, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_KAN, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_KAN, dep_KAN, file = "data_KAN.RData")

#############################################################
#for SPEC
data_diff_SPEC <- test_diff(data_imp, type = "control", control = "SPEC_") 

#Denote significant proteins based on user defined cutoffs
dep_SPEC <- add_rejections(data_diff_SPEC, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_SPEC, contrast = "SPEC__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_SPEC, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_SPEC, proteins = "Q88ES8", type = "centered")

plot_cond(dep_SPEC)

#STEP7:Visualisations of data
plot_pca(dep_SPEC, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_SPEC, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_SPEC, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_SPEC, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_SPEC, dep_SPEC, file = "data_SPEC.RData")

#############################################################
#for TET
data_diff_TET <- test_diff(data_imp, type = "control", control = "TET_") 

#Denote significant proteins based on user defined cutoffs
dep_TET <- add_rejections(data_diff_TET, alpha = 0.05, lfc = log2(1.5))

#Visualization of the results
plot_volcano(dep_TET, contrast = "SPEC__vs_K_", label_size = 2, add_names = TRUE)

plot_single(dep_TET, proteins = c("Q88ES8", "Q88QK8"))

plot_single(dep_TET, proteins = "Q88ES8", type = "centered")

plot_cond(dep_TET)

#STEP7:Visualisations of data
plot_pca(dep_TET, x = 1, y = 2, n = 500, point_size = 4)
plot_cor(dep_TET, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_TET, type = "centered", kmeans = TRUE,k = 6, col_limit = 4, show_row_names = FALSE,indicate = c("condition", "replicate"))
plot_heatmap(dep_TET, type = "contrast", kmeans = TRUE,k = 6, col_limit = 10, show_row_names = FALSE)
save(data_se, data_norm, data_imp, data_diff_SPEC, dep_TET, file = "data_TET.RData")




load("data.RData")
