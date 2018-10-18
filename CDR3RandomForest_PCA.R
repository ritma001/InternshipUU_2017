# fig 8 #
# random forest #
################################################################################
vdjRF <- function(df, fold, plot = F){
  # load necessary R packages
  lapply(c("randomForest", "pROC", "caret", "dplyr", "ROCR", "RColorBrewer"), 
         library, character.only = TRUE)
  
  # ensure caterogical outcome
  df$dominant <- factor(df$dominant)
  
  # create k-fold cross validated subsets
  set.seed(fold)
  train_index <- createFolds(df$dominant, k = fold, list = TRUE, returnTrain = TRUE)
  
  # initialise the color vector
  cols <- colorRampPalette(brewer.pal(11, "Set3"))(fold)
  # random sampling test set and training set
  for(i in 1:fold){
    train <- df[unlist(train_index[i]), ]
    test <- df[-unlist(train_index[i]), ]
    
    # randomly select ID of train data (so, 25% ID and 75% SD for training)
    train_sd <- sample_n(train[train$dominant == "SD", ], sum(train$dominant == "ID")*3)
    train_i <- rbind(train_sd, train[train$dominant == "ID", ])
    
    test_sd <- sample_n(test[test$dominant == "SD", ], sum(test$dominant == "ID")*3)
    test_i <- rbind(test_sd, test[test$dominant == "ID", ])
    
    # train model
    rf.fit <- randomForest(formula = dominant ~ ., data = train_i,
                           strata = dominant, importance = TRUE,proximity = TRUE)
    # print variable importance
    print(rf.fit)
    
    #obatain accuracy
    print(rf.fit$confusion[, 'class.error'])
    
    
    # print number of ID and SD in train dataset
    cat("train:", summary(factor(train_i$dominant)))
    # estimate accuracy of the test set
  
    if(plot == F){
      ## insert the roc plot or varimp plot ##
      varImpPlot(rf.fit)
    } else {
      testSet <- test_i[,!colnames(df) == "dominant"]
      
      # calculate accuracy
      pred.class <- predict(rf.fit, newdata = testSet)
      test_i[,"dominant"]
      rightPred <- as.character(pred.class) == test_i[,"dominant"]
      print(rightPred)
      accuracy <- sum(rightPred)/nrow(testSet)
      print(accuracy)
      
      # plot ROC curve
      rf.pred <- predict(rf.fit, newdata = testSet, type = "prob")[, 1]
      if(i == 1){
        plot.roc(unlist(test_i[,"dominant"]), rf.pred, grid = TRUE, print.auc = TRUE, 
                 print.auc.y = 0.1, print.auc.x = 0, col = cols[i], lwd = 4,
                 cex.lab = 1.5, cex.axis=1.5, print.auc.cex = 1.5)
      } else {
        plot.roc(unlist(test_i[,"dominant"]), rf.pred, grid = TRUE, print.auc = TRUE, 
                 print.auc.y = 0.1+(i-1)/20, print.auc.x = 0, col = cols[i], 
                 lwd = 4, print.auc.cex = 1.5, add = TRUE)
      } 
    }
  }
}

################################################################################
# train on 6 original feature from previous analysis; V, J, length of CDR3 and 
# nonVJ, hydropathy, affinity and aa content
# train on 6 features
vdjRF(a2NLV[, c(18:19, 23, 28, 34, 46, 57)], 5, T)
aucNLV <- mean(c(0.724,0.612,0.602,0.653,0.625))
# NLV classifier on 24 features #
vdjRF(a2NLV[, c(18:19, 23, 28, 34, 36:53, 57)], 5)
aucNLV<- meanc(c(0.748, 0.612, 0.517, 0.714, 0.639))
accNLV.id <- mean(1-c(0.85, 0.74, 0.78, 0.85, 0.68))*100
accNLV.sd <- mean(1-c(0.09, 0.09, 0.09, 0.04, 0.06))*100

################################################################################
vdjRF(a2GLC[, c(18:19, 23, 28, 34, 46, 57)], 5, T)
aucGLC <- mean(c(0.74, 0.57,0.447,0.555,0.516))

# GLC classifers on 24 features #
vdjRF(a2GLC[, c(18:19, 23, 28, 34, 36:53, 57)], 5, T)
aucGLC <- meanc(c(0.693, 0.442, 0.558, 0.620, 0.555))
accGLC.id <- mean(1-c(0.94, 0.82, 0.88, 0.85, 0.77))*100
accGLC.sd <- mean(1-c(0.11, 0.07, 0.10, 0.07, 0.04))*100

################################################################################
vdjRF(a2GIL[, c(18:19, 23, 28, 34, 46, 57)], 5, T)
aucGIL <- mean(c(0.396,0.73,0.693,0.646,0.854))
# GLC classifers on 24 features #
vdjRF(a2GIL[, c(18:19, 23, 28, 34, 36:53, 57)], 5, T)
aucGIL <- mean(c(0.688,0.520,0.767,0.604,0.458))
accGIL.id <- mean(1-c(0.83, 0.65, 0.88, 0.72, 0.72))*100
accGIL.sd <- mean(1-c(0.07, 0.06, 0.10, 0.07, 0.07))*100


################################################################################
vdjRF.epitope <- function(df, fold, plot = F){
  # load necessary R packages
  lapply(c("randomForest", "pROC", "caret", "dplyr", "ROCR", "RColorBrewer"), 
         library, character.only = TRUE)
  
  # ensure caterogical outcome
  df$Epitope <- factor(df$Epitope)
  
  # create k-fold cross validated subsets
  set.seed(fold)
  train_index <- createFolds(df$Epitope, k = fold, list = TRUE, returnTrain = TRUE)
  
  # initialise the color vector
  cols <- colorRampPalette(brewer.pal(11, "Set3"))(fold)
  # random sampling test set and training set
  for(i in 1:fold){
    train <- df[unlist(train_index[i]), ]
    test <- df[-unlist(train_index[i]), ]
    
    # train model
    rf.fit <- randomForest(formula = Epitope ~ ., data = train,
                           strata = Epitope, importance = TRUE,proximity = TRUE)
    # print variable importance
    print(rf.fit)
    
    #obatain accuracy
    print(rf.fit$confusion[, 'class.error'])
    
    accuracies <- c()
    
    if(plot == F){
      varImpPlot(rf.fit)
    } else {
      ## insert the roc plot or varimp plot ##
      rf.pred <- predict(rf.fit, newdata = test[,!colnames(df) == "Epitope"], type = "prob")[, 1]
      
      # calculate accuracy
      pred.class <- predict(rf.fit, newdata = test[,!colnames(df) == "Epitope"])
      rightPred <- as.character(pred.class) == test[,"Epitope"]
      print(rightPred)
      accuracy <- sum(rightPred)/nrow(test)
      print(accuracy)
      accuracies <- c(accuracies, accuracy)
   
      # plot ROC curve
      if(i == 1){
        plot.roc(unlist(test[,"Epitope"]), rf.pred, grid = TRUE, print.auc = TRUE, 
                 print.auc.y = 0.1, print.auc.x = 0, col = cols[i], lwd = 4,
                 cex.lab = 1.5, cex.axis=1.5, print.auc.cex = 1.5)
      } else {
        plot.roc(unlist(test[,"Epitope"]), rf.pred, grid = TRUE, print.auc = TRUE, 
                 print.auc.y = 0.1+(i-1)/20, print.auc.x = 0, col = cols[i], 
                 lwd = 4, print.auc.cex = 1.5, add = TRUE)
      }
    }
    print(accuracies)
  }
}

x <- droplevels(do.call(rbind, list(a2GIL, a2GLC, a2NLV)))
# 
vdjRF.epitope(x[, c(10, 18:19, 23, 28, 46, 57)], 5)
# 24 features
colnames(x)["CDR3_AA_GRAVY", "V.new", "J.new"] <-  c("Hydropathy", "V.gene", "J.gene")
vdjRF.epitope(x[!duplicated(x[, "CDR3"]), c(18:19, 23, 28, 10, 36:53, 57)], 5, T)
aucCD8 <- mean(c(0.901,0.836,0.840,0.865, 0.828))
################################################################################

# vdj of kidera factor
vdj.pca <- function(vdj, color.factor){
  # import library #
  library(ggfortify)
  library(pca3d)
  
  # convert type factor to character and interger to numeric
  vdj <- convertType(is.integer, as.numeric, vdj)
  vdj <- convertType(is.factor, as.character, vdj)
  
  # pca #
  num.factors <- sapply(vdj, is.numeric)
  vdj.pca <- prcomp(na.omit(vdj[,num.factors]), center = T, scale. = T)
  
  # 3dplot
  ##pca3d(vdj.pca, group=factor(vdj$dominant))
 
  # plot and print the result
  p <- autoplot(vdj.pca, 
                data = na.omit(vdj), 
                colour = color.factor) + custom_theme
  
  print(vdj.pca)
  print(summary(vdj.pca))
  return(p)
}













