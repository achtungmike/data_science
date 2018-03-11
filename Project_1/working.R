# # library(tidyverse)
# # library(modelr)
# # 
# # nfold <-3
# # 
# # setwd("C:/Users/Mike/OneDrive/Mike/UC/Spring_2018/BMIN7054_Data_Science/hw/data_science/homework_1")
# # file <- "her2_lumA_test.txt"
# # 
# # # Load our gene sets:
# # # First 50 genes listed in her2_lumA.txt 
# # first50 = c('ARHGEF10L','HIF3A','RNF17','RNF10','RNF11','RNF13','GTF2IP1',
# #             'REM1','MTVR2','RTN4RL2','C16orf13','C16orf11','FGFR1OP2','TSKS',
# #             'ATRX','PMM2','ASS1','NCBP1','ZNF709','ZNF708','RBM14','NCBP2',
# #             'DISC1','CAMK1','RPL37','SPR','ZNF700','ZNF707','CAMK4','ZNF704',
# #             'LOC339240','GOLGA6B','RNF115','RNF112','ZC3H14','SPN','HMGCLL1',
# #             'NACAP1','LRRTM1','GRIN1','RBMY1A3P','DHX8','DHX9','LOC441204',
# #             'TCOF1','LRRTM3','NUP98','XPC','SLC12A2','GRINA')
# # 
# # # Pam50 genes provided in file PAM50.txt and force into char vector
# # pam = as.vector(as.matrix(read_tsv(file = 'PAM50.txt', col_names = F)))
# # 
# # # Read in the expression data and stuff in df, remove extraneous columns
# # input <- read_tsv(file = file, skip = 1, col_names = F)
# # colnames(input) <- as.vector(read_tsv(file = file, col_names = F, n_max = 1))
# # input <- select(input, -one_of(c('PROBE','Name_GeneSymbol','ID_geneid','DESCRIPTION')))
# # 
# # # Remove the subtype row and stuff it into another df
# # # DO WE HAVE TO REMOVE ID COLUMN IN SUBTYPE?
# # subtype <- slice(input, 1)
# # input <- slice(input, 2:n())
# # 
# # # Get the sample names for each subtype, will come in handy later.
# # lum_samples <- colnames(select(subtype, which(subtype == 'LumA')))
# # her_samples <- colnames(select(subtype, which(subtype == 'Her2')))
# # 
# # # Categorize by PAM50 call.  First select rows based on gene list, then 
# # # select columns based on PAM50 call.
# # ##################
# # # First 50 Genes #
# # ##################
# # lum <- filter(input, id %in% first50) %>%
# #   select(lum_samples)
# # 
# # her <- filter(input, id %in% first50) %>%
# #   select(her_samples)
# #              
# # 
# # 
# # lum_folds <- createFolds(lum_samples, k = nfold)
# # her_folds <- createFolds(her_samples, k = nfold)
# # 
# # for (i in 1:nfold)
# # {
# #   # Create the test sets for each subtype
# #   test_lum <- lum[,unlist(lum_folds[i])]
# #   test_her <- her[,unlist(her_folds[i])]
# #   
# #   # Create training sets for each subtype
# #   # -- choose only columns (cancer samples) NOT assigned to a selected group (test_group) 
# #   trainingA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
# #   trainingB <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
# #   
# #   
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ###############
# # # PAM50 Genes #
# # ###############
# # 
# # lum <- filter(input, id %in% pam) %>%
# #   select(which(subtype == 'LumA'))
# # 
# # her <- filter(input, id %in% pam) %>%
# #   select(which(subtype == 'Her2'))
# # 
# # 
# # 
# # 
# # ###############
# # # ALL   Genes #
# # ###############
# # lum <- select(input, which(subtype == 'LumA'))
# # 
# # her <- select(input, which(subtype == 'Her2'))
# 
# 
# 
# 
# 
# 
# library(readr)
# library(knitr)
# 
# # Code re-used from in class example
# # provided by Dr. Michal Kouril, CCHMC
# 
# # Set variables with parameter values
# file <- "her2_lumA.txt"   # this is the file we will be reading from
# 
# # Load our gene sets:
# # First 50 genes listed in her2_lumA.txt 
# first50 = c('ARHGEF10L','HIF3A','RNF17','RNF10','RNF11','RNF13','GTF2IP1',
#             'REM1','MTVR2','RTN4RL2','C16orf13','C16orf11','FGFR1OP2','TSKS',
#             'ATRX','PMM2','ASS1','NCBP1','ZNF709','ZNF708','RBM14','NCBP2',
#             'DISC1','CAMK1','RPL37','SPR','ZNF700','ZNF707','CAMK4','ZNF704',
#             'LOC339240','GOLGA6B','RNF115','RNF112','ZC3H14','SPN','HMGCLL1',
#             'NACAP1','LRRTM1','GRIN1','RBMY1A3P','DHX8','DHX9','LOC441204',
#             'TCOF1','LRRTM3','NUP98','XPC','SLC12A2','GRINA')
# 
# 
# # Pam50 genes provided in file PAM50.txt and force into char vector
# pam = as.vector(as.matrix(read_tsv(file = 'PAM50.txt', col_names = F)))
# 
# 
# # Pull in header
# header <- scan(file, nlines = 1, what = character())
# 
# # read the entire file skipping first two lines (headers)
# data <- read.table(file, skip = 1, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
# 
# # set the column names to the header
# names(data) <- header
# 
# # read the second header line -- it contains the type of cancer
# sub <- scan(file, skip = 1, nlines = 1, what = character())
# 
# 
# # Her2: 47 vars
# # LumA: 201 vars
# 
# rownames(data) <- data$id
# 
# # First 50 genes
# lum_50 <- data[data$id %in% first50, sub == 'LumA' ]
# her_50 <- data[data$id %in% first50, sub == 'Her2']
# 
# # All genes
# lum_all <- data[,sub == 'LumA']
# her_all <- data[,sub == 'Her2']
# 
# # PAM50 genes
# lum_pam <- data[data$id %in% pam, sub == 'LumA']
# her_pam <- data[data$id %in% pam, sub == 'Her2']
# 
# 
# 
# 
# 
# 
# 
# 
# 
# cv_log <- function(lum_dat, her_dat, nfolds)
# {
#   
#   #LumA will be 0
#   #Her2 will be 1
#   
#   
#   FOLD = nfolds
#   
#   # transpose dataframes and add labels 
#   lum_labs <- rownames(lum_dat)
#   lum_dat <- sapply(lum_dat, as.numeric)
#   rownames(lum_dat) <- lum_labs
#   lum_dat <- data.frame(t(lum_dat))
#   lum_dat$call <- 0
#   lum_dat$sample <- rownames(lum_dat)
#   
#   her_labs <- rownames(her_dat)
#   her_dat <- sapply(her_dat, as.numeric)
#   rownames(her_dat) <- her_labs
#   her_dat <- data.frame(t(her_dat))
#   her_dat$call <- 1
#   her_dat$sample <- rownames(her_dat)
#   
#   # Merge into one df.
#   combo <- merge.data.frame(lum_dat, her_dat, all = TRUE)
#   rownames(combo) <- combo$sample
#   combo$sample <- NULL
#   combo$call <- as.factor(combo$call)
#   
#  
#   
#   
#   # Create array for results
#   res <- array()
#   
#   groups <- split(sample(rownames(combo)), 1+(seq_along(rownames(combo)) %% FOLD))
#   
#   for (it in 1:FOLD) {
#     
#     # create test set for each cancer type
#     # -- choose only columns (cancer samples) assigned to a selected group (test_group) 
#     test <- combo[rownames(combo) %in% unlist(groups[it]),]
#     
#     # create training set for each cancer type
#     # -- choose only columns (cancer samples) NOT assigned to a selected group (test_group) 
#     train <- combo[!(rownames(combo) %in% unlist(groups[it])),]
#   
# 
#     mod <- glm(call ~., data = train, family = binomial, control = list(maxit = 50))
#     
#     preds <- predict(mod, newdata = test, type = "response")
#     
#     pred_sub <- ifelse(preds >= 0.5, 1, 0)
#     
#     mat <- table(test$call, pred_sub)
#     rownames(mat) <- c("Obs LumA", "Obs Her2")
#     colnames(mat) <- c("Pred LumA", "Pred Her2")
#     
#     mat
#     # save to the result array so we can average it across all test groups
#     res[it] <- (mat[2] + mat[3]) / (mat[1] + mat[2] + mat[3] + mat[4])
#   }
#   
#   # calculate a mean and standard deviation
#   mu = mean(res)
#   sig = sd(res)
#   
#   return(c(mu, sig))
#   
# }
# 
# 
# # Prep output dataframe
# rows = c('3-fold', '5-fold', '10-fold')
# genes = c('first50', 'All', 'PAM50')
# 
# res = data.frame(matrix(ncol = 3, nrow = 3))
# colnames(res) <- genes
# rownames(res) <- rows
# 
# 
# # Loop to change num of folds and store output in result df
# for (fold in c(3,5,10))
# {
#   curr_fold = (paste(fold,'-fold', sep = ""))
#   
#   # first50
#   out = cv(lum_50, her_50, fold)
#   res[curr_fold, 'first50'] <- paste(round(out[1], digits = 2), "\u00b1", round(out[2], digits = 2))
#   
#   # All genes
#   #out = cv(lum_all, her_all, fold)
#   #res[curr_fold, 'All'] <- paste(round(out[1], digits = 2), "\u00b1", round(out[2], digits = 2))
#   
#   # PAM50
#   out = cv(lum_pam, her_pam, fold)
#   res[curr_fold, 'PAM50'] <- paste(round(out[1], digits = 2), "\u00b1", round(out[2], digits = 2))
#   
# }
# 
# kable(res)





b_dat = data.frame("Method" = character() , "Mean" = integer(), "min" = integer(), "max" = integer(), "group" = integer(), stringsAsFactors=FALSE)
lr_dat = data.frame("Method" = character() , "Mean" = integer(), "min" = integer(), "max" = integer(), "group" = integer(),  stringsAsFactors=FALSE)
# Prep output dataframe
rows = c('3-fold', '5-fold', '10-fold')
genes = c('first50', 'PAM50', 'All')

cent_res = data.frame(matrix(ncol = 3, nrow = 3))
colnames(cent_res) <- genes
rownames(cent_res) <- rows

lr_res = data.frame(matrix(ncol = 3, nrow = 3))
colnames(lr_res) <- genes
rownames(lr_res) <- rows


## CENTROID RUNS

# Loop to change num of folds and store output in result df
for (fold in c(3,5,10))
{
  curr_fold = (paste(fold,'-fold', sep = ""))
  
  # first50
  out = cv(lum_50, her_50, fold)
  cent_res[curr_fold, 'first50'] <- paste(round((out[1] * 100), digits = 2), "\u00b1", round((out[2]*100), digits = 2))
  mean = round((out[1] * 100), digits = 2)
  sd = round((out[2]*100), digits = 2)
  min = mean - sd
  max = mean + sd
  tmp = list((paste('First50_', fold,'-fold', sep = "")), mean,min, max , 0)
  b_dat[nrow(b_dat) + 1, ] <- tmp
  
  
  # PAM50
  out = cv(lum_pam, her_pam, fold)
  cent_res[curr_fold, 'PAM50'] <- paste(round((out[1] * 100), digits = 2), "\u00b1", round((out[2]*100), digits = 2))
  mean = round((out[1] * 100), digits = 2)
  sd = round((out[2]*100), digits = 2)
  min = mean - sd
  max = mean + sd

  tmp = list((paste('PAM50_', fold,'-fold', sep = "")), mean,min, max, 1)
  b_dat[nrow(b_dat) + 1, ] <- tmp
  
  # All genes
  out = cv(lum_all, her_all, fold)
  cent_res[curr_fold, 'All'] <- paste(round((out[1] * 100), digits = 2), "\u00b1", round((out[2]*100), digits = 2))
  
  
  mean = round((out[1] * 100), digits = 2)
  sd = round((out[2]*100), digits = 2)
  min = mean - sd
  max = mean + sd
  tmp = list((paste('AllGenes_', fold,'-fold', sep = "")), mean, min, max, 2)
  b_dat[nrow(b_dat) + 1, ] <- tmp
  
  
  
  
}

## Logistic Regression RUNS

# Loop to change num of folds and store output in result df
for (fold in c(3,5,10))
{
  curr_fold = (paste(fold,'-fold', sep = ""))
  
  # first50
  out = cv_log(lum_50, her_50, fold)
  lr_res[curr_fold, 'first50'] <- paste(round((out[1] * 100), digits = 2), "\u00b1", round((out[2]*100), digits = 2))
  
  mean = round((out[1] * 100), digits = 2)
  sd = round((out[2]*100), digits = 2)
  min = mean - sd
  max = mean + sd
  tmp =  c((paste('first50_', fold,'-fold', sep = "")), mean, min, max, 0)
  lr_dat[nrow(lr_dat) + 1, ] <- tmp
  
  # PAM50
  out = cv_log(lum_pam, her_pam, fold)
  lr_res[curr_fold, 'PAM50'] <- paste(round((out[1] * 100), digits = 2), "\u00b1", round((out[2]*100), digits = 2))
  mean = round((out[1] * 100), digits = 2)
  sd = round((out[2]*100), digits = 2)
  min = mean - sd
  max = mean + sd
  tmp =  c((paste('PAM50_', fold,'-fold', sep = "")), mean, min, max, 1)
  lr_dat[nrow(lr_dat) + 1, ] <- tmp
  
  
  # All genes - not running logistic regression - stack overflow error.
  #out = cv(lum_all, her_all, fold)
  lr_res[curr_fold, 'All'] <- "N/A"
  
}

library(ggplot2)
b_dat$group <- as.factor(b_dat$group)
p <- ggplot(b_dat) + geom_crossbar(aes(ymin = min, ymax = max, x = Method, y = Mean, fill = group))
p

lr_dat$Mean <- as.numeric(lr_dat$Mean)
lr_dat$min <- as.numeric(lr_dat$min)
lr_dat$max <- as.numeric(lr_dat$max)
lr_dat$group <- as.factor(lr_dat$group)
lr <- ggplot(lr_dat) + geom_crossbar(aes(ymin = min, ymax = max, x = Method, y = Mean, fill = group))
lr
kable(cent_res, align = c('r', 'r', 'r'), caption = "Centroid Based Cross-Validation Results")
kable(lr_res, align = c('r', 'r', 'r'), caption = "Logistic Regression Based Cross-Validation Results")
