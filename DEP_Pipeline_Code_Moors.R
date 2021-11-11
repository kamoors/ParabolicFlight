library(limma)
library(gplots)
library(ggfortify)
library(ggplot2)
library(readxl)
library(writexl)
library("DEP")
library(rlang)
library(SNFtool)
library(factoextra)
library(ggpubr)
library(dplyr)
library(VennDiagram)
library(ellipse)
library(RColorBrewer)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(qvalue)


#SET WD TO DATA LOCATION!
#########################################################
#Analysis description:

# 1) Load data
# 2) Filter data to have max. 1 NA in each treatment
# 3) Find if all treatment samples appear of reasonable quality
# 4) Remove GC4 (PCA shows large distance from main cluster)
# 5) Run DEP pipeline
# 6) Manually adjust qvalue (DEP does something weird when adjusting - qvalue package gives very different results)

#Load data----------------------------------------------------------

data_full <- read_xlsx("200217_Dra parabolic flight proteomics.xlsx")

zeroes <- data_full == 0
data_full[zeroes] <- NA





label <- colnames(data_full[,7:18])

condition <- gsub("_\\d", "", colnames(data_full[,7:18]))

replicate <- rep(1:4, 3)




Experimental_design <- data.frame(label, stringsAsFactors = F)
Experimental_design$condition <- condition
Experimental_design$replicate <- replicate
# Experimental_design$time <- Time
# Experimental_design$location <- Location

# Filter data to have max. 1 NA in each treatment------------------------------------


nas <- is.na(data_full[,7:18])
NoNA <-  rowSums(nas) < 2
data_noNA <- data_full[NoNA,]


data_filtered_full <- data_full

data_filtered_full$KEEP <- 0

for (i in 1:nrow(nas)) {
  
  if (sum(nas[i,1:4]) < 2 & sum(nas[i,5:8]) < 2 & sum(nas[i,9:12]) < 2) {
    
    
    data_filtered_full$KEEP[i] = 1
    
  }
  
}


data_filtered_full <- data_filtered_full[data_filtered_full$KEEP == 1,]



#DEP with filtered data-----------------------------------------------------


# Remove GC4:

Without_GC4 <- data_filtered_full[,-18]



#Re-assign EX. Design without GC4

label_2 <- colnames(data_full[,7:17])

condition_2 <- gsub("_\\d", "", colnames(data_full[,7:17]))

replicate_2 <- c(1,2,3,4,1,2,3,4,1,2,3)


Experimental_design2 <- data.frame(label_2, stringsAsFactors = F)
Experimental_design2$condition <- condition_2
Experimental_design2$replicate <- replicate_2

colnames(Experimental_design2) <- c("label", "condition", "replicate")


data_filtered_full$Protein %>% duplicated() %>% any() #Find if there are duplicate protein IDs





data_unique <- make_unique(Without_GC4, "Protein", "Gene ID", delim = ";")


data_se <- make_se(data_unique,7:18, Experimental_design2)


# Plot a barplot of the protein identification overlap between samples

plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 3)


# Less stringent filtering:
# Filter for proteins that are identified in 4/5 replicates of at least one condition
# data_filt2 <- filter_missval(data_se, thr = 1)


# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)




# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)



# Normalize the data
data_norm <- normalize_vsn(data_filt)


# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)



# Plot a heatmap of proteins with missing values
missval_heatmap <- plot_missval(data_filt)
missval_heatmap <- plot_missval(data_norm)


# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)


#proteins with missing values have on average low intensities. This data (MNAR and close to the detection limit) should be imputed 
#by a left-censored imputation method, such as the quantile regression-based left-censored function ("QRILC") or random draws from 
#a left-shifted distribution ("MinProb" and "man"). In contrast, MAR data should be imputed with methods such as k-nearest neighbor ("knn") 
#or maximum likelihood ("MLE") functions. See the MSnbase vignette and more specifically the impute function description for more information.



# All possible imputation methods are printed in an error, if an invalid function name is given.
# impute(data_norm, fun = "")


# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "QRILC")



# Test all possible comparisons of samples

data_diff_imp <- test_diff(data_imp, type = "manual", test = c("FC_vs_GC", "X0g_vs_GC", "FC_vs_X0g"))


# Denote significant proteins based on user defined cutoffs

dep_imp <- add_rejections(data_diff_imp, alpha = 0.05)


# Plot the first and second principal components
plot_pca(dep_imp, x = 1, y = 2, n = nrow(dep_imp), point_size = 4)



# Generate a wide data.frame
df_wide <- get_df_wide(dep_imp)



q_value_0g_GC <- qvalue(df_wide$X0g_vs_GC_p.val)

q_value_FC_GC <- qvalue(df_wide$FC_vs_GC_p.val)

q_value_FC_0g <- qvalue(df_wide$FC_vs_X0g_p.val)


df_wide$x0g_GC_qvalue <- q_value_0g_GC$qvalues

df_wide$FC_GC_qvalue <- q_value_FC_GC$qvalues

df_wide$FC_x0g_qvalue <- q_value_FC_0g$qvalues



#Write file
#write_xlsx(df_wide, "df_wideFINAL.xlsx")