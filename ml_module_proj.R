#!/usr/bin/env Rscript

# FTE-35306 Machine Learning
# Bioinformatics Project

## === Call libraries ===================================================

cur_loc <- "C:/Users/steph/Desktop/Bioinformatics/FTE35306 - Machine Learning/Project"
setwd(cur_loc)

library(boot)
library(broom)
library(caret)
library(caretEnsemble)
library(class)
library(GGally)
library(MASS)
library(mlbench)
library(RColorBrewer)
library(tidyverse)

options(digits=6) # reduce number of significant figures shown in printing

## === Functions ========================================================

# Generic linear regression with no CV
linear_regression <- function(dataset, dependent_variable, formula_to_use,
                              ratio =0.75, list = F){
  
  inTraining <- createDataPartition(dataset[[dependent_variable]], p = ratio, list = list)
  train <- dataset[ inTraining,]
  test  <- dataset[-inTraining,]
  
  # Create LR models
  lm.fit <- lm(formula_to_use, data=train)
  
  train.mse <- mean((train[[dependent_variable]] 
                    - predict(lm.fit, train))^2
  )
  
  test.mse <- mean((test[[dependent_variable]] 
                    - predict(lm.fit, test))^2
  )
  
  output <- vector("list", 4)
  names(output) <- c("train.MSE", "test.MSE", "AIC", "fit")
  output$train.MSE <- train.mse
  output$test.MSE <- test.mse
  output$AIC <- AIC(lm.fit)
  output$fit <- lm.fit
  # output$R2 <- summary(lm.fit)$r.squared 
  
  return(output)
}


# Carries out linear regression using caret package; returns the model and the test RMSE
resampled_lm <- function(dataset, dependent_variable, control, predictors = ".", 
                         seed = 24, ratio = 0.75, list = FALSE){
  inTraining <- createDataPartition(dataset[[dependent_variable]], p = ratio, list = list)
  training <- dataset[ inTraining,]
  testing  <- dataset[-inTraining,]
  
  formula_to_use <- as.formula(paste(dependent_variable, "~ ", paste(predictors, collapse="+")))
  
  set.seed(seed)
  lin_fit <- train(formula_to_use, data=training, method="glm", trControl=control)
  
  predictions <- predict(lin_fit, testing)
  
  output <- vector("list", 3)
  output <- setNames(output, c("model", "training_index", "performance"))
  output$model <- lin_fit
  output$training_index <- inTraining
  output$performance <- postResample(pred = predictions, obs = testing[[dependent_variable]])
  return(output)
}



# Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## === Main program ===================================================

datafile <- "expr4T.dat"
dat <- read.table(datafile)

# Variables of interest
g1 <- "ENSG00000271043.1_MTRNR2L2"
g2 <- "ENSG00000229344.1_RP5.857K21.7"
g3 <- "ENSG00000237973.1_hsa.mir.6723"
genes_of_interest <- c(g1, g2, g3)
tissue <- "tissue"
seed <- 14
set.seed(seed)

## == Data cleaning for LR ============================================

#Last column contains tissue label
tissue_ind <- ncol(dat)-1
m <- sapply(dat[,1:tissue_ind],mean)

plot(log(m))
try_line <- log(mean(m) + 1*sd(m))
abline(h=try_line, col="red")

s <- sapply(dat[,1:tissue_ind],sd)

# Visualise the distribution of gene variance
try_line <- log(mean(s) + 1*sd(s))
plot(log(s))
abline(h=try_line, col="red")
try_line <- log(mean(s) - 1*sd(s))
abline(h=try_line, col="red")
sm <- s/m

#Filter genes based on minimum mean and stdev/mean
# mean(m)
# sd(m)
# min(m)
# 
# mean(s)
# sd(s)
# max(s)

minsm <- 1 #set some reasonable value here
minm <- 20 #set some reasonable value here

# In the following, we make sure that tissue label is kept;
# by giving it an artificial sufficiently high value for m and sm
m <- c(m,minm+1)
sm <- c(sm,minsm+1)
length(sm)
length(m)

dim(dat)
length(which(sm>minsm & m>minm))
# length(which(m>minm))
(which(sm>minsm & m>minm))
dat_2 <- dat[, which(sm>minsm & m>minm)]


# Seperate out the genes of interest and the tissue variable

for(gene in genes_of_interest){
  if(!(gene %in% names(dat_2))){
    dat_2 <- col_bind(dat_2, dat[[gene]])
  }
  
}

numeric_data <- dat_2 %>% dplyr::select(everything(), -contains("tissue"))

# Consider scaling data (standardise the data, mean = 0, sd ~ 1)
scaled_data <- scale(numeric_data) # Creates issue for logging data (no longer strictly positive)

linear_reg_data <- numeric_data[ , -which(names(numeric_data) %in% genes_of_interest)]


## === Feature Selection ==============================

# create control object for 10-fold cv
ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE,
                   number = 10)

# Fit the model with subset size corresponding to the following:
subsets <- c(1:5, 10, 15)

x <- linear_reg_data

gene_profile <- vector("list", 3)
gene_profile <- setNames(gene_profile, c(g1, g2, g3))


for(gene in genes_of_interest){
  y <- numeric_data[[gene]]
  
  set.seed(seed)
  
  gene_profile[[gene]] <- rfe(x, y,
                              sizes = subsets,
                              rfeControl = ctrl)
}


predictors(gene_profile[[g1]])
predictors(gene_profile[[g2]])
predictors(gene_profile[[g3]])

max_num_predictors <- 15

g1_num_predictors <- min(gene_profile[[g1]]$bestSubset,max_num_predictors)
g2_num_predictors <- min(gene_profile[[g2]]$bestSubset,max_num_predictors)
g3_num_predictors <- min(gene_profile[[g3]]$bestSubset,max_num_predictors)

g1_predictors <- gene_profile[[g1]]$optVariables[1:g1_num_predictors]
g2_predictors <- gene_profile[[g2]]$optVariables[1:g2_num_predictors]
g3_predictors <- gene_profile[[g3]]$optVariables[1:g3_num_predictors]

predictors <- unique(c(g1_predictors, g2_predictors, g3_predictors))


trellis.par.set(caretTheme())
ggplot(gene_profile[[g1]], type = c("g", "o")) +
  labs(title = "Feature Elimination",
       subtitle = g1)
ggplot(gene_profile[[g2]], type = c("g", "o")) +
  labs(title = "Feature Elimination",
       subtitle = g2)
ggplot(gene_profile[[g3]], type = c("g", "o")) +
  labs(title = "Feature Elimination",
       subtitle = g3)


data_to_use <- numeric_data %>% 
  dplyr::select(predictors, genes_of_interest)

## === Correlated Features ==============================================
# Remove highly correlated variables (removing those with the highest
# mean correlation from each pair)
# cof_df <- cor(linear_reg_data)
# cols_to_remove <- findCorrelation(cof_df, cutoff=0.9)
# cols_to_remove <- sort(cols_to_remove)
# 
# data.new <- linear_reg_data[,-c(cols_to_remove)]
# 
# # Combine the reduced dataset with the genes of interest and tissue
# data_to_use <- bind_cols(data.new, dat_2[, which(names(dat_2) %in% genes_of_interest)])
# 
# # final_data <- bind_cols(data.new, dat_2[, which(names(dat_2) %in% c(g1,g2, g3, tissue))])
# # 
# # write.csv(data_to_use, file = "cleaned_numeric_data.csv",row.names=FALSE)
# # write.csv(final_data, file = "final_data.csv",row.names=FALSE)

## === Linear Regression =============================================

control <- trainControl(method="repeatedcv", number=10, repeats=3) # 10 => ten-fold cv

# We will do linear regression across one tissue; choose that with the 
# most associated observations
tissue_most_sample <- names(summary(dat_2$tissue)[summary(dat_2$tissue) 
                                                  == max(summary(dat_2$tissue))])

tissue_data <- data_to_use[dat_2[["tissue"]] == tissue_most_sample,]
 
logged_data <- log(data_to_use + 1)
logged_tissue_data <- logged_data[dat_2[["tissue"]] == tissue_most_sample,]

scaled_data <- data_to_use %>%
  scale %>%
  as.data.frame

summary(scaled_data)
summary(data_to_use)

lm_results <- vector("list", 15)
lm_results <- setNames(lm_results, c(g1, g2, g3,
                                     paste0(g1, "_", tissue_most_sample),
                                     paste0(g2, "_", tissue_most_sample),
                                     paste0(g3, "_", tissue_most_sample),
                                     paste0(g1, "_scaled"),
                                     paste0(g2, "_scaled"),
                                     paste0(g3, "_scaled"),
                                     paste0(g1, "_logged"),
                                     paste0(g2, "_logged"),
                                     paste0(g3, "_logged"),
                                     paste0(g1, "_logged_", tissue_most_sample),
                                     paste0(g2, "_logged_", tissue_most_sample),
                                     paste0(g3, "_logged_", tissue_most_sample)))

datasets <- list(data_to_use, tissue_data, scaled_data, logged_data, logged_tissue_data)

# ggpairs(tissue_data)
summary(data_to_use)
summary(logged_data)
summary(scaled_data)

cv_effect_list <- vector("list", 2) %>% setNames(list("CV", "No_CV"))
# controls <- list(control, control_no_cv)
count <- -2

# Carry out LR for both no test set, no CV; CV; test set
for(i in 1:2){
  for(gene in genes_of_interest){
    for(dataset in datasets){
      count <- count + 3
      if(identical(dataset, tissue_data) | identical(dataset, logged_tissue_data)){
        loc_max_num_predictors <- floor((nrow(dataset) / 10) * 0.75)
        loc_predictors <- gene_profile[[gene]]$optVariables[1:min(gene_profile[[gene]]$bestSubset,loc_max_num_predictors)]
      } else{
        loc_predictors <- "."
      }
      
      if(i == 2){
        # Define function for generic LR
        formula_to_use <- as.formula(paste(gene, "~", paste(loc_predictors, collapse="+")))
        lm_results[[names(lm_results)[count]]] <- linear_regression(dataset, gene, formula_to_use)
      } else{ 
        # CV LR
        lm_results[[names(lm_results)[count]]] <- resampled_lm(dataset, gene, control,
                                                               predictors = loc_predictors,
                                                               seed = seed)
      }
    }
    count <- count - length(lm_results) + 1
  }
  cv_effect_list[[i]] <- lm_results
  count <- -2
}

# Inspect outcome
# No CV or test
# gene 1
summary(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_scaled$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_brain_cerebellum$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged_brain_cerebellum$fit)$r.squared

sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_scaled$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_brain_cerebellum$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged_brain_cerebellum$train.MSE)

sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_scaled$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_brain_cerebellum$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000271043.1_MTRNR2L2_logged_brain_cerebellum$test.MSE)

# gene 2
summary(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_scaled$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_brain_cerebellum$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged_brain_cerebellum$fit)$r.squared

sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_scaled$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_brain_cerebellum$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged_brain_cerebellum$train.MSE)

sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_scaled$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_brain_cerebellum$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000229344.1_RP5.857K21.7_logged_brain_cerebellum$test.MSE)

# gene 3
summary(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_scaled$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_brain_cerebellum$fit)$r.squared
summary(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged_brain_cerebellum$fit)$r.squared

sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_scaled$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_brain_cerebellum$train.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged_brain_cerebellum$train.MSE)

sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_scaled$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_brain_cerebellum$test.MSE)
sqrt(cv_effect_list$No_CV$ENSG00000237973.1_hsa.mir.6723_logged_brain_cerebellum$test.MSE)

# Inspect the outcome
# scaling has no impact on performace (all of each variable
# undergoes the same linear transformation, does not change the fit)
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_scaled$model
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged$model
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_brain_cerebellum$model
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged_brain_cerebellum$model

# See impact of scaling in coefficients: now the same unit (in this case were already,
# but this is generally true)
# Intercept is value when predictors are at their mean value
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model$finalModel
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_scaled$model$finalModel


cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$performance
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_scaled$performance
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged$performance
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_brain_cerebellum$performance
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged_brain_cerebellum$performance

# gene 2 results
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7$model
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_scaled$model
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_logged$model
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_brain_cerebellum$model
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_logged_brain_cerebellum$model

cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7$performance
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_scaled$performance
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_logged$performance
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_brain_cerebellum$performance
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7_logged_brain_cerebellum$performance

# gene 3 result
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_scaled$model
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged$model
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_brain_cerebellum$model
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged_brain_cerebellum$model


cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$performance
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_scaled$performance  
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged$performance
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_brain_cerebellum$performance
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged_brain_cerebellum$performance

# Plot variable importance
ggplot(varImp(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model)) +
  labs(title = "Feature Importance",
       subtitle = g1)

ggplot(varImp(cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7$model)) +
  labs(title = "Feature Importance",
       subtitle = g2)

ggplot(varImp(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model)) +
  labs(title = "Feature Importance",
       subtitle = g3)

g1_important_pred <- varImp(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model)$importance
g1_important_pred$Names <- rownames(g1_important_pred)
g1_important_pred <- g1_important_pred %>% arrange(desc(Overall)) %>% filter(Overall > 25)


g2_important_pred <- varImp(cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7$model)$importance
g2_important_pred$Names <- rownames(g2_important_pred)
g2_important_pred <- g2_important_pred %>% arrange(desc(Overall)) %>% filter(Overall > 25)


g3_important_pred <- varImp(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model)$importance
g3_important_pred$Names <- rownames(g3_important_pred)
g3_important_pred <- g3_important_pred %>% arrange(desc(Overall)) %>% filter(Overall > 25)

g1_imp_1 <- g1_important_pred$Names[1]
g1_imp_2 <- g1_important_pred$Names[2]

g2_imp_1 <- g2_important_pred$Names[1]

g3_imp_1 <- g3_important_pred$Names[1]

# important_gene <- "ENSG00000225972.1_MTND1P23" # based on above plots

g1_reduced_lm <- resampled_lm(data_to_use, g1, control, seed = seed, predictors = g1_imp_1)
g1_reduced_lm$model

g1_reduced_lm_v2 <- resampled_lm(data_to_use, g1, control, seed = seed, predictors = c(g1_imp_1, g1_imp_2))
g1_reduced_lm_v2$model

# g1_reduced_lm_v3 <- resampled_lm(data_to_use, g1, control, seed = seed, predictors = c(g2, g3,important_gene))
# g1_reduced_lm_v3$model

g2_reduced_lm <- resampled_lm(data_to_use, g2, control, seed = seed, predictors = g2_imp_1)
g2_reduced_lm$model

# g2_reduced_lm_v2 <- resampled_lm(data_to_use, g2, control, seed = seed, 
#                                  predictors = c(g3,important_gene))
# g2_reduced_lm_v2$model

g3_reduced_lm <- resampled_lm(data_to_use, g3, control, seed = seed, predictors = g3_imp_1)
g3_reduced_lm$model



# Compare reduced models to saturated
g1_reduced_lm$performance
g1_reduced_lm_v2$performance
# g1_reduced_lm_v3$performance
cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$performance

g2_reduced_lm$performance
# g2_reduced_lm_v2$performance
cv_effect_list$CV$ENSG00000229344.1_RP5.857K21.7$performance # note loss of quality in saturated model

g3_reduced_lm$performance
# g3_reduced_lm_v2$performance
# g3_reduced_lm_v3$performance
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$performance


cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model$finalModel$coefficients
cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_scaled$model$finalModel$coefficients

par(mfrow=c(2,2))
plot(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model$finalModel, which = c(1,2))
plot(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged$model$finalModel, which = c(1,2))
par(mfrow=c(1,1))


p1 <- ggplot(varImp(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723$model)) +
        labs(title = "Feature Importance",
             subtitle = g3)

p2 <- ggplot(varImp(cv_effect_list$CV$ENSG00000237973.1_hsa.mir.6723_logged$model)) +
  labs(title = "Feature Importance",
       subtitle = paste("Log(", g3, ")"))

p1
p2
# multiplot(p1, p2, cols=2)

par(mfrow=c(2,2))
plot(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model$finalModel, which = c(1,2))
plot(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged$model$finalModel, which = c(1,2))
par(mfrow=c(1,1))


p1 <- ggplot(varImp(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2$model)) +
  labs(title = "Feature Importance",
       subtitle = g3)

p2 <- ggplot(varImp(cv_effect_list$CV$ENSG00000271043.1_MTRNR2L2_logged$model)) +
  labs(title = "Feature Importance",
       subtitle = paste("Log(", g3, ")"))

p1
p2


## === Classification Functions ========================================

# Function to return the accuracy of a classification model
# along with other useful outputs
logistic_accuracy <- function(model, train, test, classes, dependent_variable,
                              type="response", threshold=0.5){
  # Inspect train
  train_probs <- predict(model, train, type=type)
  train_pred <- rep(classes[[1]] ,nrow(train))
  
  # Inspect test
  model_probs <- predict(model, test, type=type)
  model_pred <- rep(classes[[1]] ,nrow(test))
  
  train_pred[train_probs > threshold] <- classes[[2]]
  model_pred[model_probs > threshold] <- classes[[2]]
  
  train_acc <- mean(train_pred==train[[dependent_variable]])
  model_acc <- mean(model_pred==test[[dependent_variable]])
  
  output <- setNames(list(model_probs, model_pred, model_acc, train_acc), 
                     c("probabilities", "predictions", "test.accuracy", "train.accuracy"))
  
  return(output)
}

lda_basic <- function(model, training, testing, dependent_variable, classes){
  # Inspect train
  train_probs <- predict(lda.fit, training)
  train_pred <- rep(classes[[1]] ,nrow(training))
  
  # Inspect test
  model_probs <- predict(lda.fit, testing)
  model_pred <- rep(classes[[1]] ,nrow(testing))
  
  train_pred <- train_probs$class
  model_pred <- model_probs$class
  
  train_acc <- mean(train_pred==training[[dependent_variable]])
  model_acc <- mean(model_pred==testing[[dependent_variable]])
  
  output <- setNames(list(model_probs, model_pred, model_acc, train_acc), 
                     c("probabilities", "predictions", "test.accuracy", "train.accuracy"))
}

# Carry out KNN
knn_basic <- function(training, testing, dependent_variable, k = 5){
  
  # Create index to access dependent_variable
  dependent_variable_ind <- match(dependent_variable, names(training))
  
  knn_train <- training %>% dplyr::select(-dependent_variable_ind)
  knn_test <- testing %>% dplyr::select(-dependent_variable_ind)
  
  knn.fit.train <- knn(knn_train, knn_train, 
                       training[[dependent_variable]], k = k)
  # table(knn.fit.train, training$tissue)
  knn.train.accuracy <- mean(knn.fit.train == training[[dependent_variable]])
  
  knn.fit.test <- knn(knn_train, knn_test, training[[dependent_variable]], k = k)
  # table(knn.fit.test,testing$tissue)
  knn.test.accuracy <- mean(knn.fit.test == testing[[dependent_variable]])
  knn_out <- list(model = knn.fit.test,
                  train.accuracy = knn.train.accuracy,
                  test.accuracy = knn.test.accuracy,
                  k = k)
  return(knn_out)
}


# Finds levels or declares error if appropriate
find_levels_raise_error <- function(data, dependent_variable, pair){
  # Check if pair is declared correctly
  if(is.null(pair)){
    if(length(levels(data[[dependent_variable]])) > 2){
      stop(paste("Please declare pair.", dependent_variable, "has more than  2 levels. "))
    } else{
      pair <- levels(data[[dependent_variable]])
    }
  }
  return(pair)
}


# Creates the data objects for rfe
create_x_y_input_for_rfe <- function(data, dependent_variable, classifier, pair=NULL){
  
  dependent_variable_ind <- match(dependent_variable, names(data))
  
  if(classifier){
    
    pair <- find_levels_raise_error(data, dependent_variable, pair)
    
    y <- data %>%
      filter(get(dependent_variable) %in% pair) %>%
      droplevels
    
    y <- y[[dependent_variable]]
    
    # Create data frame (or tibble) of predictors
    
    x <- data %>%
      filter(get(dependent_variable) %in% pair) %>% 
      dplyr::select(-dependent_variable_ind)
    
  } else{
    y <- data[[dependent_variable]]
    
    # Create data frame (or tibble) of predictors
    x <- data %>%
      dplyr::select(-dependent_variable_ind)
  }
  
  output <- vector("list", 2) %>% setNames(list("x", "y"))
  output$x <- x
  output$y <- y
  return(output)
}


# Returns the set of useful predictors and the rfe object used to fid them
find_predictors <- function(data, dependent_variable, 
                            subsets, control, max_num_predictors=15,
                            seed=24, classifier=TRUE, pair=NULL){
  
  # Create vector of dependent variable
  inputs <- create_x_y_input_for_rfe(data, dependent_variable, classifier, pair=pair)
  x <- inputs$x
  y <- inputs$y
  
  # Set seed for recreateable analysis
  set.seed(seed)
  
  loc_name <- paste(pair, collapse="_v_")
  
  # Carry out recursive feature selection using caret
  feature_selection <- rfe(x, y, sizes = subsets, rfeControl = ctrl)
  
  # Define the number of features to select
  num_predictors <- min(feature_selection$bestSubset, max_num_predictors)
  
  # Select the best features (limited by the input max_num_predictors)
  useful_predictors <- feature_selection$optVariables[1:num_predictors]
  
  # Declare output vector
  output <- vector("list", 2) %>% setNames(c("rfe", "predictors"))
  output$rfe <- feature_selection
  output$predictors <- useful_predictors
  
  return(output)
}


## === Classification Feature Selection ==================================

# Declare the pairs of interest
pair1 <- c("brain_amygdala", "brain_hippocampus")
pair2 <- c("brain_hippocampus",  "brain_nucleusaccumbens")
pair3 <- c("brain_spinalcord", "brain_substantianigra")
pair4 <- c("brain_cerebellarhemisphere", "brain_cerebellum")
pair5 <- c("brain_cerebellum", "brain_amygdala")

plot(dat_2[,ncol(dat_2)]) # check frequencies, range of [55, 120] (approx)
# This will limit the number of predictors we can use

pairs_of_interest <- list(pair1, pair2, pair3, pair4, pair5)

## Feature selection based on logistic regression
# create control object for 10-fold cv

# Consider trying rfFuncs too
ctrl <- rfeControl(functions = lrFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE,
                   number = 10)


# Fit the model with subset size corresponding to the following:
subsets <- c(1:5, 7, 10)

# Create the names for each list for the 5 pairings
pairing_names <- rep(0, length(pairs_of_interest))
for(i in 1:length(pairs_of_interest)){
  pairing_names[[i]] <- paste(pairs_of_interest[[i]], collapse="_v_")
}

# Declare the objects to hold the rfe object and useful predictors
# for each pair.
tissue_profile <- vector("list", length(pairs_of_interest)) %>% 
  setNames(pairing_names)

useful_predictors <- vector("list", length(pairs_of_interest)) %>% 
  setNames(pairing_names)

# Set the max number of predictors for the classification problems
# due to the limited number of samples this is lower than for regression
max_num_predictors_class <- 12

# Create the data frame
tissue_data <- numeric_data %>% 
  bind_cols(dplyr::select(dat_2, tissue))

# Find all the relevant predictors across each pair of tissues (based on LR)
for(pair in pairs_of_interest){
  
  feature_selection <- find_predictors(tissue_data, "tissue", 
                                       subsets, control, 
                                       max_num_predictors=max_num_predictors_class,
                                       seed=seed, classifier=TRUE, pair=pair)
  
  
  loc_name <- paste(pair, collapse="_v_")
  
  tissue_profile[[loc_name]] <- feature_selection$rfe
  useful_predictors[[loc_name]] <- feature_selection$predictors
}

# Create a list of all the relevant predictors
tissue_predictors <- unique(unlist(useful_predictors))

useful_predictors[[paste(pair5, collapse="_v_")]]

# This is probably unnecessary
tissue_data <- numeric_data %>% 
  dplyr::select(tissue_predictors) %>% 
  bind_cols(dplyr::select(dat_2, tissue))

# still will want to select the relevant subset for each pair.

## === Classification ======================================================

pair_logistic_info <- vector("list", length(pairing_names)) %>%
  setNames(pairing_names)

pair_lda_info <- vector("list", length(pairing_names)) %>%
  setNames(pairing_names)

pair_knn_info <- vector("list", length(pairing_names)) %>%
  setNames(pairing_names)

# Idea here is to inspect the logistic regression model for each pair
for(pair in pairs_of_interest){
  
  # Create name to access apropriate level of list
  loc_name <- paste(pair, collapse="_v_")
  
  pair_data <- tissue_data %>% filter(tissue %in% pair) %>% droplevels %>%
    dplyr::select(useful_predictors[[loc_name]], tissue)
  
  # Create the test and training data
  set.seed(seed)
  inTraining <- createDataPartition(pair_data[["tissue"]], p = 0.75, list = FALSE)
  training <- pair_data[ inTraining,]
  testing  <- pair_data[-inTraining,]
  
  # Create and record models
  formula_to_use <- as.formula(paste(tissue, "~ ", 
                                     paste(useful_predictors[[loc_name]], collapse="+")
  )
  )
  
  log_fit <- glm(formula_to_use, family="binomial", data=training)
  
  pair_logistic_info[[loc_name]] <- logistic_accuracy(log_fit, training, 
                                                      testing, pair, "tissue",
                                                      type="response", threshold=0.5)
 
  lda.fit <- lda(formula_to_use, data = pair_data, subset = inTraining)
  
  pair_lda_info[[loc_name]] <- lda_basic(lda.fit, training, testing, "tissue", pair)
  
  pair_knn_info[[loc_name]] <- knn_basic(training, testing, "tissue", k = 5)
}


pair_lda_info$brain_amygdala_v_brain_hippocampus$test.accuracy
pair_lda_info$brain_hippocampus_v_brain_nucleusaccumbens$test.accuracy
pair_lda_info$brain_spinalcord_v_brain_substantianigra$test.accuracy
pair_lda_info$brain_cerebellarhemisphere_v_brain_cerebellum$test.accuracy
pair_lda_info$brain_cerebellum_v_brain_amygdala$test.accuracy

pair_lda_info$brain_amygdala_v_brain_hippocampus$train.accuracy
pair_lda_info$brain_hippocampus_v_brain_nucleusaccumbens$train.accuracy
pair_lda_info$brain_spinalcord_v_brain_substantianigra$train.accuracy
pair_lda_info$brain_cerebellarhemisphere_v_brain_cerebellum$train.accuracy
pair_lda_info$brain_cerebellum_v_brain_amygdala$train.accuracy

inspect_pair5_data <- tissue_data %>%
  filter(tissue %in% pair5) %>%
  select(useful_predictors[[5]], tissue)
ggpairs(inspect_pair5_data)
# labs(title = "Pairwise plots",
#      subtitle = "Amygdala, Cerebellum")

inspect_pair5_data_specific_preds <- tissue_data %>%
  filter(tissue %in% pair5) %>%
  select(useful_predictors[[5]][c(2,4,6,8,9)], tissue)

ggpairs(inspect_pair5_data_specific_preds) + 
  labs(title = "Pairwise plots",
       subtitle = "Amygdala, Cerebellum")

inspect_pair2_data <- tissue_data %>%
  filter(tissue %in% pair2) %>%
  select(useful_predictors[[2]], tissue)
ggpairs(inspect_pair2_data)
# labs(title = "Pairwise plots",
#      subtitle = "Amygdala, Cerebellum")

inspect_pair2_data_specific_preds <- tissue_data %>%
  filter(tissue %in% pair2) %>%
  select(useful_predictors[[2]][c(2,6,7)], tissue)

ggpairs(inspect_pair2_data_specific_preds) + 
  labs(title = "Pairwise plots",
       subtitle = "Hippocampus, Nucleus accumbens")

## === Classification with Cross Validation ================================

metric <- "Accuracy"

# Declare methods of interest for classification
classifier_methods <- c("lda", "knn", "rpart", "adaboost", "rf", "svmRadial", "glm") # glm needs special handling
classifier_methods_names <- c("lda", "knn", "tree", "boost", "rf", "svm", "log_reg")
classifier_methods <- setNames(classifier_methods, classifier_methods_names)

# For each model we want to save certain associated items
items_of_interest <- c("model", "confustion.table", "test.resample")

# Exclude LR as uses a different style to others (so needs to be outside the loop)
similar_methods <- classifier_methods
similar_methods$log_reg <- NULL

classification_results <- vector("list", length(pairs_of_interest)) %>% 
  setNames(pairing_names)

for(pair in names(classification_results)){
  classification_results[[pair]] <- vector("list", length(classifier_methods)) %>% 
    setNames(classifier_methods_names)
  
  for(method in names(classification_results[[pair]])){
    classification_results[[pair]][[method]] <- vector("list", length(items_of_interest)) %>% 
      setNames(items_of_interest)
  }
}

# For each pair:
# This takes ages. Consider saving the object and reading in.
for(pair in pairs_of_interest){
  loc_name <- paste(pair, collapse="_v_")
  
  pair_data <- tissue_data %>% filter(tissue %in% pair) %>% droplevels %>%
    dplyr::select(useful_predictors[[loc_name]], tissue)
  
  set.seed(seed)
  inTraining <- createDataPartition(pair_data[["tissue"]], p = 0.75, list = FALSE)
  training <- pair_data[ inTraining,]
  testing  <- pair_data[-inTraining,]
  
  formula_to_use <- as.formula(paste(tissue, "~ ", 
                                     paste(useful_predictors[[loc_name]], collapse="+")
                                     )
                               )
  
  set.seed(seed)
  
  classification_results[[loc_name]][["log_reg"]][["model"]] <- train(formula_to_use, 
                                                                      data=training,
                                                                      method="glm", 
                                                                      trControl=control, 
                                                                      family="binomial", 
                                                                      metric=metric)
  
  predictions <- predict(classification_results[[loc_name]][["log_reg"]][["model"]], testing)
  classification_results[[loc_name]][["log_reg"]][["test.resample"]] <- postResample(pred = predictions, 
                                                                                     obs = testing[[tissue]])

  classification_results[[loc_name]][["log_reg"]][["confustion.table"]] <- confusionMatrix(data = predictions,
                                                                                           reference = testing$tissue)

  for(class_method in names(similar_methods)){
    
    curr_method <- similar_methods[[class_method]]
    
    set.seed(seed)
    classification_results[[loc_name]][[class_method]][["model"]] <- train(formula_to_use,
                                                                           data=training,
                                                                           method=curr_method, 
                                                                           trControl=control,
                                                                           metric=metric)
                          
    predictions <- predict(classification_results[[loc_name]][[class_method]][["model"]], testing)
    classification_results[[loc_name]][[class_method]][["test.resample"]] <- postResample(pred = predictions,
                                                                                          obs = testing[[tissue]])

    classification_results[[loc_name]][[class_method]][["confustion.table"]] <-  confusionMatrix(data = predictions, 
                                                                                                 reference = testing$tissue)
  }
}

# Save the relevant plots (note this is from CV not test data)
graph_folder <- "Classification Graphs"

for(i in 1:length(classification_results)){
  
  pair <- classification_results[[i]]
  
  # collect the resampling results
  resamps <- resamples(list(LOG = pair[["log_reg"]][["model"]],
                            LDA = pair[["lda"]][["model"]],
                            KNN = pair[["knn"]][["model"]],
                            SVM = pair[["svm"]][["model"]],
                            Tree = pair[["tree"]][["model"]],
                            RF = pair[["rf"]][["model"]],
                            Boost = pair[["boost"]][["model"]]))
  
  # Investigate comparison
  
  # First a box plot for the accuracy and Kappa value of each model
  theme1 <- trellis.par.get()
  theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
  theme1$plot.symbol$pch = 16
  theme1$plot.line$col = rgb(1, 0, 0, .7)
  theme1$plot.line$lwd <- 2
  trellis.par.set(theme1)
  
  plot_title <- names(classification_results)[i]
  plot_name <- paste0("boxplot_", plot_title, ".jpg")
  mypath <- file.path(cur_loc, graph_folder, plot_name)
  
  jpeg(file=mypath)
  print(bwplot(resamps, layout = c(2, 1), main=plot_title))
  dev.off()
  
  # Then a plot of the mean accuracy of each model with error bars
  trellis.par.set(caretTheme())
  
  plot_name <- paste0("dotplot_", plot_title, ".jpg")
  
  mypath <- file.path(cur_loc, graph_folder, plot_name)
  jpeg(file=mypath)
  print(dotplot(resamps, metric = "Accuracy", main=plot_title))

  dev.off()
 
}

# We notice pair5 has several perfect predictors; investigate by inspecting data
pair5_data <- tissue_data %>%
  dplyr::select(useful_predictors[[paste(pair5, collapse="_v_")]],
                tissue) %>%
  dplyr::filter(tissue %in% pair5) %>%
  droplevels

pairs(pair5_data, main="Pair-wise scatterplots for pair 5 data")


# classification_results$brain_amygdala_v_brain_hippocampus$svm$confustion.table
# 
# classification_results$brain_cerebellum_v_brain_amygdala$rf$confustion.table
# classification_results$brain_cerebellarhemisphere_v_brain_cerebellum$rf$confustion.table

# classification_results$brain_amygdala_v_brain_hippocampus$svm$confustion.table



# write.table(classification_results$brain_amygdala_v_brain_hippocampus$svm$confustion.table$table,
            # file="sample_confusion_table.csv",sep=",")

tissue_choice <- classification_results$brain_hippocampus_v_brain_nucleusaccumbens
model_boost <- tissue_choice$boost$model
model_svm <-tissue_choice$svm$model
model_log <- tissue_choice$log_reg$model

importance_threshold <- 25

# Investigate variable importace within the models
plot(varImp(model_boost))
boost_important_pred <- varImp(model_boost)$importance
boost_important_pred$Names <- rownames(boost_important_pred)
boost_important_pred <- boost_important_pred %>% 
  arrange(desc(brain_nucleusaccumbens), desc(brain_hippocampus)) %>% 
  filter(brain_nucleusaccumbens > importance_threshold 
         | brain_hippocampus > importance_threshold)

# Investigate variable importace within the models
plot(varImp(model_svm))
svm_important_pred <- varImp(model_svm)$importance
svm_important_pred$Names <- rownames(svm_important_pred)
svm_important_pred <- svm_important_pred %>% 
  arrange(desc(brain_nucleusaccumbens), desc(brain_hippocampus)) %>% 
  filter(brain_nucleusaccumbens > importance_threshold
         | brain_hippocampus > importance_threshold)

# Investigate variable importace within the models
plot(varImp(model_log))
log_important_pred <- varImp(model_log)$importance
log_important_pred$Names <- rownames(log_important_pred)
log_important_pred <- log_important_pred %>% 
  arrange(desc(Overall)) %>%
  filter(Overall > importance_threshold)

important_class_predictors <- unique(c(log_important_pred$Names,
                                       svm_important_pred$Names,
                                       boost_important_pred$Names))

log_important_pred$Names[!(log_important_pred$Names %in% important_class_predictors)]
svm_important_pred$Names[!(svm_important_pred$Names %in% important_class_predictors)]
boost_important_pred$Names[!(boost_important_pred$Names %in% important_class_predictors)]

# The models have the same predictros as important, but the LR model has a different
# order and weighting





# 
# plot(varImp(lda_fit))
# lda_important_pred <- varImp(lda_fit)$importance
# lda_important_pred$Names <- rownames(lda_important_pred)
# lda_important_pred <- lda_important_pred %>% 
#   arrange(desc(brain_amygdala), desc(brain_hippocampus)) %>% 
#   filter(brain_amygdala > 25 | brain_hippocampus > 25) # though both classes have the same associated importance
# 
# 
# plot(varImp(knn_fit)) # identical to LDA for importance
# knn_important_pred <- varImp(knn_fit)$importance
# knn_important_pred$Names <- rownames(knn_important_pred)
# knn_important_pred <- knn_important_pred %>% 
#   arrange(desc(brain_amygdala), desc(brain_hippocampus)) %>% 
#   filter(brain_amygdala > 25 | brain_hippocampus > 25)

## === Unsupervised methods ========================================

## === PCA =========================================================

pca_function <- function(data_mat, scale=F)
  # Function for PCA that also returns the variance explained by the first two components
{
  prcomp_out<-prcomp(data_mat, scale=scale)
  var<-(prcomp_out$sdev)^2
  sumvar<-sum(var)
  expl1<-100*var[1]/sumvar
  expl2<-100*var[2]/sumvar
  explboth<-sum(expl1,expl2)
  output <- setNames(list(prcomp_out, sumvar, expl1, expl2, explboth), 
                     c("prcomp", "sumvar", "expl1", "expl2", "explboth"))
  return(output)
}

# return num_components variance explained by the components or prcomp_obj
pca_variance_expl <- function(prcomp_obj, num_components){
  var<-(prcomp_obj$sdev)^2
  sumvar<-sum(var)
  
  output <- setNames(rep(0, num_components), paste0("expl_", 1:num_components))
  
  for(i in 1:num_components){
    curr_level <- paste0("expl_", i)
    output[[curr_level]] <- 100*var[i]/sumvar
  }
  return(output)
}


# Biplot using output of above function
biplot_labelled <- function(pca_function_output)
{
  xlab<-paste('PC1 ',format(pca_function_output$expl1,digits=3),'%')
  ylab<-paste('PC2 ',format(pca_function_output$expl2,digits=3),'%')
  mainlab<-paste('PCA biplot genes (',format(
    pca_function_output$explboth,digits=3),'%)',sep="")
  
  biplot(pca_function_output$prcomp,xlab=xlab,ylab=ylab, main=mainlab)
}


heatmap_genes <- function(gene_data, mypalette){
  # Make a matrix and also the transpose of the data matrix
  data_mat <- as.matrix(gene_data)
  tdata_mat<-t(data_mat)
  
  # Heatmap of the data matrix, in original order
  heatmap(data_mat,Rowv=NA, Colv=NA, col=mypalette, scale='none')
  # Heatmap of the data matrix, rows and columns ordered according to a hierarchical clustering procedure
  heatmap(data_mat, col=mypalette, scale='none')
  
  # Heatmap of the correlations among the different conditions
  heatmap(cor(data_mat),symm=TRUE,col=mypalette)
  # Heatmap of the correlations among the genes
  heatmap(cor(tdata_mat), symm=TRUE,col=mypalette)
}

# Some heatmapping

# Choose palette going from Red to Blue for heatmaps
mypalette<-brewer.pal(10,"RdBu")
# Reverse the palette so that red is for high values, blue for low
mypalette<-mypalette[10:1]

# pca_matrix <- linear_reg_data %>% 
#   as.matrix


# Inspect PCA for the pre- and post- feature selection data
pca_matrix_post_fe <- data_to_use %>% 
  as.matrix

pca_matrix_pre_fe <- linear_reg_data %>% 
  as.matrix

pca_matrix_all <- dat %>% dplyr::select(-tissue) %>% 
  as.matrix

# summary(pca_matrix_all)
# str(pca_matrix_all)

heatmap_genes(cor(pca_matrix_post_fe), mypalette)

unscaled_pca_post_fe <- pca_function(pca_matrix_post_fe)
biplot_labelled(unscaled_pca_post_fe)

unscaled_pca_pre_fe <- pca_function(pca_matrix_pre_fe)

unscaled_pca_all <- pca_function(pca_matrix_all)
scaled_pca_all <- pca_function(pca_matrix_all, scale=T)

unscaled_pca_post_fe$explboth
unscaled_pca_pre_fe$explboth
unscaled_pca_all$explboth


pca_var_expl <- pca_variance_expl(unscaled_pca_all$prcomp, 5)
sum(pca_var_expl)

pca_var_expl <- pca_variance_expl(scaled_pca_all$prcomp, 5)
sum(pca_var_expl)

screeplot(unscaled_pca_all$prcomp, type="lines", main="PCA for raw data (unscaled)")
screeplot(scaled_pca_all$prcomp, type="lines", main="PCA for raw data (scaled)")

pr.var <- scaled_pca_all$prcomp$sdev ^2
pve <- pr.var/sum(pr.var )

plot(pve[1:15] , xlab=" Principal Component ", ylab=" Proportion of
     Variance Explained ", ylim=c(0,1) ,type="b")
plot(cumsum (pve[1:15]), xlab=" Principal Component ", ylab ="
       Cumulative Proportion of Variance Explained ", ylim=c(0,1) ,
       type="b")

pr.var <- unscaled_pca_all$prcomp$sdev ^2
pve <- pr.var/sum(pr.var )


plot(pve[1:15] , xlab=" Principal Component ", ylab=" Proportion of
     Variance Explained ", ylim=c(0,1) ,type="b")
plot(cumsum (pve[1:15]), xlab=" Principal Component ", ylab ="
     Cumulative Proportion of Variance Explained ", ylim=c(0,1) ,
     type="b")

# biplot(unscaled_pca_all$prcomp , scale =0)

library(ISLR)

Cols <- function (vec){
  cols <- rainbow (length(unique(vec)))
  return (cols[as.numeric (as.factor(vec))])
}



example_tissue <- dat %>%
  dplyr::select(tissue) %>%
  dplyr::filter(tissue %in% pair1) %>%
  as.matrix

example_gene <- dat %>%
  dplyr::filter(tissue %in% pair1) %>%
  dplyr::select(-tissue) %>%
  as.matrix


pr.out <- prcomp(example_gene, scale=T)
summary (pr.out)

# 

par(mfrow =c(1,2))

plot(pr.out$x [,1:2], col =Cols(example_tissue), pch =19,
       xlab ="Z1",ylab="Z2")
plot(pr.out$x[,c(1,3) ], col =Cols(example_tissue), pch =19,
       xlab ="Z1",ylab="Z3")


plot(unscaled_pca_all$prcomp$x[,1:2], col=Cols(dat$tissue), pch =19,
     xlab ="Z1",ylab="Z2")
plot(unscaled_pca_all$prcomp$x[,c(1,3) ], col=Cols(dat$tissue), pch =19,
     xlab ="Z1",ylab="Z3")
par(mfrow= c(1,1))





## === K Means Clustering ==========================================

# Enables investigating across multiple values of k
# Gives warnings if k is not a vector of numbers, but still works
k_means <- function(data.matrix, k_vec, nstart=50){
  kclusts <- data.frame(k=k_vec) %>%
    group_by(k) %>%
    do(kclust=kmeans(data.matrix, .$k, nstart=nstart))
  
  output <- vector("list", 3) %>% 
    setNames(c("clusters", "assignments", "clusterings"))
  
  # tidy summarises on a per-cluster level
  output$clusters <- kclusts %>% group_by(k) %>% do(tidy(.$kclust[[1]]))
  
  # augment adds the point classifications to the original dataset
  output$assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], points.matrix))
  
  output$clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))
  
  return(output)
}

# Returns a plot with <shape> at the cluster centre
k_means_plot <- function(k_means_obj, x_var, y_var,
                         palette_set = "Set1", shape="x", size=10){
  # Set up a larger palette of colours (isn't working?)
  getPalette <- colorRampPalette(brewer.pal(9, palette_set))
  
  # Plot the original points, with each point colored according to the original cluster:
  p1 <- ggplot(k_means_obj$assignments, 
               mapping = aes_string(x = x_var, y = y_var)) +
    geom_point(mapping = aes(colour=.cluster)) +
    facet_wrap(~ k) +
    scale_fill_manual(values = getPalette(max(k_means_obj$clusterings$k)))
  
  names(k_means_obj$clusters)[2:45] <- names(k_means_obj$assignments)[2:45]
  
  # We can then add the centers of the cluster using the data from tidy:
  p2 <- p1 + geom_point(data=k_means_obj$clusters, size=size, shape=shape)
  return(p2)
  
}

set.seed(seed)
num_tissues <- length(unique(tissue_data$tissue))

# Seperate the data into a matrix
points.matrix <- tissue_data %>% select(-tissue) %>% as.matrix

# Carry out k means clustering
k_means_tissue <- k_means(points.matrix, 1:(num_tissues+5))

# Investigate the optimal number of clusters
plot(k_means_tissue$clusterings$tot.withinss, type='l')

# This plot is misleading as k=1 is so awful. Exclude this value:
k_means_tissue$clusterings$totss

plot_var <- k_means_tissue$clusterings$tot.withinss / k_means_tissue$clusterings$totss

# Set up a plot
par(mar = c(2, 3, 2, 1),
    mgp = c(2, 0.4, 0), 
    las = 1,
    tck = -.01,
    xaxs = "i", yaxs = "i")


plot(2:18, plot_var[2:18], type='l', lwd=2, col="blue",
     xlab = "Number of clusters", ylab = "Within SS",
     axes = FALSE, # Don't plot the axes
     frame.plot = FALSE, # Remove the frame 
     panel.first = abline(h = seq(0.05, 0.2, 0.05), col = "grey80"))
abline(h = seq(0.025, 0.2, 0.05), lty=2, col = "grey80")

title("K Means Clustering: Within SS as a function of k", adj = 1, 
      cex.main = 0.8, font.main = 2, col.main = "black")

at <- pretty(2:18)
mtext(side = 1, text = at, at = at, 
      col = "grey20", line = 1, cex = 0.9)

at <- pretty(plot_var[2:18])
mtext(side = 2, text = at, at = at, col = "grey20", line = 1, cex = 0.9)

abline(v = c(10, 12, 14), lty=2, col = "red") # three possible choices for k

k_means_tissue$clusterings$tot.withinss


# Graph the various combinations of variables and the associated clusters
# name_combinations <- combn(names(tissue_data), 2) # too many; I did this, but probably less interesting to you
graph_folder <- "K Means Graphs"
name_combinations <- list(c(g1, g2), c(g1, g3), c(g2, g3))

# If correcting this, consider using only a subsection of the combinations as it takes awhile
for(i in 1:ncol(name_combinations)){
  x_var <- name_combinations[1,i]
  y_var <- name_combinations[2,i]
  title <- paste("K-means clustering")
  plot_name <- paste0("kmeans_", "_", x_var, "_", y_var, ".jpg")
  mypath <- file.path(cur_loc, graph_folder, plot_name)
  
  jpeg(file=mypath)
  plot(k_means_plot(k_means_tissue, x_var,y_var))
  dev.off()
  
}

## === Hierarchical Clustering =====================================

# Cluster based on the Euclidean distance between the conditions

hier_clust <- function(distance, num_clusters, linkage="complete"){
  clusters <- hclust(distance, method = linkage)
  clusterCut <- cutree(clusters, num_clusters)
  return(clusterCut)
}

distEucl<-dist(points.matrix)
distcorr<-as.dist(1-cor(points.matrix))

library(factoextra)

# Elbow method
fviz_nbclust(points.matrix, hcut, method = "wss", k=18, nstart=50) +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(points.matrix, hcut, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(seed)
fviz_nbclust(points.matrix, hcut, k.max = 18,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")
abline(v=c(7,13,16), lty=2, col="blue")

install.packages("NbClust", dep=T)
library(NbClust)



clusters <- hclust(distEucl, method = "average")
clusterCut <- cutree(clusters, 11)
plot(clusters)
plot(as.dendrogram(clusters), ylim = c(0,10000), xlim=c(1163, 1300))

# abline(v=1163, lty=2, col="red")
table(clusterCut, tissue_data$tissue)


distEucl<-dist(points.matrix)
hclust(distEucl)
plot(hclust(distEucl))
abline(h=14.5)

# Consider cutting tree using pre-defined number of clusters
clusterCut <- cutree(hclust(distEucl), 13)

table(clusterCut, tissue_data$tissue)


# Cluster based on the distance of correlations between the genes
distcorr<-as.dist(1-cor(points.matrix))
hclust(distcorr)
plot(hclust(distcorr))

clusterCut <- cutree(hclust(distcorr), 13)



# Add line of cutoff (this was decided for the low number of clusters of sufficient size)
abline(h=1.275)

hclust(distcorr, method='single')
plot(hclust(distcorr))
abline(h=1.275)

hclust(distcorr, method='average')
plot(hclust(distcorr))
abline(h=1.275)
