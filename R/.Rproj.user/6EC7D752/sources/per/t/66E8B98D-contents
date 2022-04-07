
######################
#### LIBRARY LOAD ####
######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")
install.packages("RcppCNPy")
install.packages("reticulate")

library(ropls)
library(caret)
library(reticulate)

#use_condaenv("r-reticulate", required = TRUE)


######################
#### DATA LOADING ####
######################


reactome_labels <- read.csv('../python/data/reactome_labels.csv', header=FALSE)
reactome_neglogp <- read.csv('../python/data/metrics/neglogpvalue_reactome.csv', row.names = 1)
reactome_closeness <- read.csv('../python/data/metrics/process_closeness.csv', row.names = 1)
reactome_betweenness <- read.csv('../python/data/metrics/process_betweenness.csv', row.names = 1)
reactome_rwr<- read.csv('../python/data/metrics/process_rwr.csv', row.names = 1)

reactome_neglogp_80 <- read.csv('../python/data/metrics/process_ppi80_hyper.csv')
reactome_closeness_80 <- read.csv('../python/data/metrics/process_ppi80_closeness.csv')
reactome_betweenness_80 <- read.csv('../python/data/metrics/process_ppi80_betweenness.csv')
reactome_rwr_80 <- read.csv('../python/data/metrics/process_ppi80_rwr.csv')

rownames(reactome_labels) <- row.names(reactome_neglogp)

##########################
#### ADJACENCY MATRIX ####
##########################

adj_matrix <- np$load('../python/data/adjacency_matrix.npy')

##########################
#### OPLS-DA FUNCTION ####
##########################

Matt_Coef <- function (conf_matrix)
{
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  return(mcc_final)
}


opls_da <- function(metrics, labels){
  columns <- c()
  tp <- c()
  fp <- c()
  tn <- c()
  fn <- c()
  recall <- c()
  precision <- c()
  f_measure <- c()
  accuracy <- c()
  mcc <- c()
  significant <- c()
  for(column in 1:ncol(labels)){
    cat("Current working on process number", column)
    cat("\n")
    columns[column] <- names(metrics)[column]
    column_labels <- factor(labels[,column], levels=c(1,0))
    train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
    process.oplsda <- try(opls(metrics, column_labels, orthoI = NA, predI=1, subset=train.index, crossvalI = 10), silent=TRUE)
    significant[column] <- TRUE
    if("try-error" %in% class(process.oplsda)) {
      significant[column] <- FALSE
      process.oplsda <- try(opls(metrics, column_labels, orthoI = 9, predI=1, subset=train.index, crossvalI = 10), silent=TRUE)
    } 
    if("try-error" %in% class(process.oplsda)) {
      significant[column] <- FALSE
      process.oplsda <- opls(metrics, column_labels, orthoI = 1, predI=1, subset=train.index, crossvalI = 10)}
    
    x_pred <- table(column_labels[-train.index], predict(process.oplsda, metrics[-train.index,]))
    result <- confusionMatrix(x_pred)
    tp[column] <- result$table[1,1]
    tn[column] <- result$table[2,2]
    fp[column] <- result$table[1,2]
    fn[column] <- result$table[2,1]
    accuracy[column] <- result$byClass['Prevalence']
    precision[column] <- result$byClass['Pos Pred Value']
    recall[column] <- result$byClass['Sensitivity']
    f_measure[column] <- 2 * ((precision[column] * recall[column]) / (precision[column] + recall[column]))
    mcc[column] <- Matt_Coef(result)
  }
  opls_da.resuls <- data.frame(columns, tp, fp, tn, fn, recall, precision, accuracy, f_measure, mcc, significant, row.names = 1)
  opls_da.resuls
}

# Yet to be implemented

opls_da_multiple_df <- function(metrics, labels){
  opls_da.resuls <- list()
  for(n in 1:(length(unique(metrics[,1]))-1)){
    columns <- c()
    recall <- c()
    precision <- c()
    f_measure <- c()
    significant <- c()
    metrics_values <- metrics[metrics$X == n,3:435]
    proteins <- metrics[metrics$X == n,2]
    cat("Currently working on reduction number", n)
    cat("\n\n")
    for(column in 1:ncol(labels)){
      cat("Currently working on process number", column)
      cat("\n")
      metrics_values_row <- data.frame(metrics_values[,column])
      columns[column] <- names(metrics_values)[column]
      column_labels <- factor(labels[proteins,column], levels=c(1,0))
      train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
      process.oplsda <- try(opls(metrics_values_row, column_labels, orthoI = NA, predI=1, subset=train.index, crossvalI = 10), silent=TRUE)
      significant[column] <- TRUE
      if("try-error" %in% class(process.oplsda)) {
        significant[column] <- FALSE
        process.oplsda <- try(opls(metrics_values_row, column_labels, orthoI = 9, predI=1, subset=train.index, crossvalI = 10), silent=TRUE)
      } 
      if("try-error" %in% class(process.oplsda)) {
        significant[column] <- FALSE
        process.oplsda <- opls(metrics_values_row, column_labels, orthoI = 1, predI=1, subset=train.index, crossvalI = 10)}
      test_values <- data.frame(metrics_values_row[-train.index,1])
      x_pred <- table(column_labels[-train.index], predict(process.oplsda, test_values))
      result <- confusionMatrix(x_pred)
      precision[column] <- result$byClass['Pos Pred Value']
      recall[column] <- result$byClass['Sensitivity']
      f_measure[column] <- 2 * ((precision[column] * recall[column]) / (precision[column] + recall[column]))
    } 
    opls_da.resuls[[n+1]] <- data.frame(columns, recall, precision, f_measure, significant, row.names = 1)
  }
  opls_da.resuls <- bind_rows(opls_da.resuls)
  opls_da.resuls
}


#################################
#### REACTOME HYPERGEOMETRIC ####
#################################


opls_da_hyper.results <- opls_da(reactome_neglogp, reactome_labels)
write.csv(opls_da_hyper.results,"../python/data/clf/reactome_hyper_oplsda.csv")

opls_da_hyper.results <- read.csv('../python/data/clf/reactome_hyper_oplsda.csv')
reactome_hyper.recall <- opls_da_hyper.results$recall
reactome_hyper.precision <- opls_da_hyper.results$precision
reactome_hyper.f_measure <- opls_da_hyper.results$f_measure

jpeg("hyper_bp.jpg")
reactome_hyper.bp <- boxplot(reactome_hyper.recall, reactome_hyper.precision, reactome_hyper.f_measure,
                       main = "Metrics",
                       names = c("recall", "precision", "f_measure"))
dev.off()
reactome_hyper.bp$stats[3,]

# Yet to be implemented

opls_da_hyper_80.results <- opls_da_multiple_df(reactome_neglogp_80, reactome_labels)

write.csv(opls_da_hyper_80.results,"../python/data/clf/reactome_hyper_80_oplsda.csv", row.names = FALSE)

opls_da_hyper_80.results <- read.csv('../python/data/clf/reactome_hyper_80_oplsda.csv')
reactome_hyper_80.recall <- opls_da_hyper_80.results$recall
reactome_hyper_80.precision <- opls_da_hyper_80.results$precision
reactome_hyper_80.f_measure <- opls_da_hyper_80.results$f_measure

jpeg("hyper_80_bp.jpg")
reactome_hyper_80.bp <- boxplot(reactome_hyper_80.recall, reactome_hyper_80.precision, reactome_hyper_80.f_measure,
                             main = "Metrics",
                             names = c("recall", "precision", "f_measure"))
dev.off()
reactome_hyper_80.bp$stats[3,]


############################
#### REACTOME CLOSENESS ####
############################


opls_da_closeness.results <- opls_da(reactome_closeness, reactome_labels)
write.csv(opls_da_closeness.results,"../python/data/clf/reactome_closeness_oplsda.csv", row.names = TRUE)

opls_da_closeness.results <- read.csv('../python/data/clf/reactome_closeness_oplsda.csv')
reactome_closeness.recall <- opls_da_closeness.results$recall
reactome_closeness.precision <- opls_da_closeness.results$precision
reactome_closeness.f_measure <- opls_da_closeness.results$f_measure

jpeg("closeness_bp.jpg")
reactome_closeness.bp <- boxplot(reactome_closeness.recall, reactome_closeness.precision, reactome_closeness.f_measure,
                                 main = "Metrics",
                                 names = c("recall", "precision", "f_measure"))
dev.off()
reactome_closeness.bp$stats[3,]

##############################
#### REACTOME BETWEENNESS ####
##############################

opls_da_betweenness.results <- opls_da(reactome_betweenness, reactome_labels)
write.csv(opls_da_betweenness.results,"../python/data/clf/reactome_betweenness_oplsda.csv", row.names = TRUE)

opls_da_betweenness.results <- read.csv('../python/data/clf/reactome_betweenness_oplsda.csv')
reactome_betweenness.recall <- opls_da_betweenness.results$recall
reactome_betweenness.precision <- opls_da_betweenness.results$precision
reactome_betweenness.f_measure <- opls_da_betweenness.results$f_measure

jpeg("betweenness_bp.jpg")
reactome_betweenness.bp <- boxplot(reactome_betweenness.recall, reactome_betweenness.precision, reactome_betweenness.f_measure,
                                 main = "Metrics",
                                 names = c("recall", "precision", "f_measure"))
dev.off()
reactome_betweenness.bp$stats[3,]

######################
#### REACTOME RWR ####
######################

opls_da_rwr.results <- opls_da(reactome_rwr, reactome_labels)
write.csv(opls_da_rwr.results,"../python/data/clf/reactome_rwr_oplsda.csv", row.names = TRUE)

opls_da_rwr.results <- read.csv('../python/data/clf/reactome_rwr_oplsda.csv')
reactome_rwr.recall <- opls_da_rwr.results$recall
reactome_rwr.precision <- opls_da_rwr.results$precision
reactome_rwr.f_measure <- opls_da_rwr.results$f_measure

jpeg("rwr_bp.jpg")
reactome_rwr.bp <- boxplot(reactome_rwr.recall, reactome_rwr.precision, reactome_rwr.f_measure,
                                 main = "Metrics",
                                 names = c("recall", "precision", "f_measure"))
dev.off()
reactome_rwr.bp$stats[3,]


# Yet to be implemented

opls_da_rwr_80.results <- opls_da_multiple_df(reactome_rwr_80, reactome_labels)
write.csv(opls_da_hyper_80.results,"../python/data/clf/reactome_rwr_80_oplsda.csv", row.names = FALSE)

opls_da_rwr_80.results <- read.csv('../python/data/clf/reactome_rwr_80_oplsda.csv')
reactome_rwr_80.recall <- opls_da_rwr_80.results$recall
reactome_rwr_80.precision <- opls_da_rwr_80.results$precision
reactome_rwr_80.f_measure <- opls_da_rwr_80.results$f_measure

jpeg("rwr_80_bp.jpg")
reactome_rwr_80.bp <- boxplot(reactome_rwr_80.recall, reactome_rwr_80.precision, reactome_rwr_80.f_measure,
                                main = "Metrics",
                                names = c("recall", "precision", "f_measure"))
dev.off()
reactome_rwr_80.bp$stats[3,]




