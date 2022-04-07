library(ropls)
library(caret)
library(foreach)
library(doParallel)
library(reticulate)
library(igraph)

source_python('../python/src/features/disease_process_proteins.py')
source_python('../python/src/features/metrics_functions.py')
np <- import("numpy")

####CUSTOM FUNCTIONS####

oplsda.fs <- function(metrics, labels){
  myCluster <- makeCluster(12, # number of cores to use
                             type = "PSOCK") # type of cluster
  registerDoParallel(myCluster)
  fs <- foreach(column = 1:ncol(labels), .combine = 'rbind') %dopar% {
    library(ropls)
    library(caret)
    column_labels <- factor(labels[,column], levels=c(1,0))
    train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
    all_labels  <- seq(1, nrow(labels))
    test.index <- setdiff(all_labels, train.index)
    if(length(test.index) == 3440){
      test.index <- append(test.index, NA)
    }
    process.oplsda <- try(opls(metrics, column_labels, orthoI = NA, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
    if("try-error" %in% class(process.oplsda)) {
      process.oplsda <- try(opls(metrics, column_labels, orthoI = 9, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
    } 
    if("try-error" %in% class(process.oplsda)) {
      process.oplsda <- opls(metrics, column_labels, orthoI = 1, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none')}
    VIP <- getVipVn(process.oplsda)

  return(list(VIP, as.vector(test.index)))
  }
  stopCluster(myCluster)
  fs
}

oplsda.fs.multiple <- function(metrics, labels){
  results <- list()
  o.df <- NULL
  test.df <- NULL
  for(n in 1:(length(unique(metrics[,1])))){
    metrics_values <- metrics[metrics$X == n-1,3:431]
    proteins <- metrics[metrics$X == n-1,2]
    myCluster <- makeCluster(12) # type of cluster
    registerDoParallel(myCluster)
    fs <- foreach(column = 1:ncol(labels), .combine = 'rbind') %dopar% {
      library(ropls)
      library(caret)
      column_labels <- factor(labels[proteins,column], levels=c(1,0))
      train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
      all_labels  <- seq(1, length(column_labels))
      test.index <- setdiff(all_labels, train.index)
      if(length(test.index) == 3440){
        test.index <- append(test.index, NA)
      }
      process.oplsda <- try(opls(metrics_values, column_labels, orthoI = NA, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
      if("try-error" %in% class(process.oplsda)) {
        process.oplsda <- try(opls(metrics_values, column_labels, orthoI = 9, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
      } 
      if("try-error" %in% class(process.oplsda)) {
        process.oplsda <- opls(metrics_values, column_labels, orthoI = 1, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none')}
      VIP <- getVipVn(process.oplsda)
      o <- order(-VIP)
      o <- o[1:11]
      if (column %in% o){}
      else {
        o <- o[1:10]
        o[11] <- column
      }
      return(list(o, as.vector(test.index)))
    }
    o.df <- rbind(o.df, fs[,1])
    test.df <- rbind(test.df, fs[,2])
    if(n==10){
      results[[1]] <- o.df
      results[[2]] <- test.df
      }
    stopCluster(myCluster)
  }
  fs.df <- as.data.frame(results[[1]][,1])
  colnames(fs.df) <- c(1,2,3,4,5,6,7,8,9,10)
  for (i in 2:ncol(labels)){
    df <- as.data.frame(results[[1]][,i])
    colnames(df) <- colnames(fs.df)
    fs.df <- rbind(fs.df, df)
  }
  test.df <- as.data.frame(lapply(results[[2]][,1], 'length<-', 17000))
  colnames(test.df) <- c(1,2,3,4,5,6,7,8,9,10)
  #for (i in 2:length(colnames(reactome_labels))){
  for (i in 2:ncol(labels)){
    df <-  as.data.frame(lapply(results[[2]][,i], 'length<-', 17000))
    colnames(df) <- colnames(test.df)
    test.df <- rbind(test.df, df)
  }
  return_results <- list(fs.df, test.df)
  return_results
}

oplsda_related.fs <- function(adj_matrix, protein_ppis, labels, metric='bridge'){
  myCluster <- makeCluster(8, # number of cores to use
                           ) # type of cluster
  registerDoParallel(myCluster)
  #registerDoSNOW(myCluster)
  #iterations <- ncol(reactome_labels)
  #pb <- txtProgressBar(max = iterations, style = 3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress = progress)
  #fs <- foreach(column = 1:ncol(labels), .combine = 'cbind', .options.snow = opts) %dopar% {
  fs <- foreach(column = 1:ncol(labels), .combine = 'cbind') %dopar% {
    library(ropls)
    library(caret)
    library(reticulate)
    source_python('../python/lib/metrics_functions.py')
    process <- colnames(labels)[column]
    column_labels <- factor(labels[,column], levels=c(1,0))
    train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
    all_labels  <- seq(1, nrow(labels))
    test.index <- setdiff(all_labels, train.index)
    if(metric=='bridge'){metrics <- bridge_scores(adj_matrix, protein_ppis, process)}
    if(metric=='scp'){metrics <- scp(protein_ppis,process)}
    if(metric=='slb'){metrics <- slb(process)}
    if(length(test.index) == 3440){
      test.index <- append(test.index, NA)
    }
    process.oplsda <- try(opls(metrics, column_labels, orthoI = NA, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
    if("try-error" %in% class(process.oplsda)) {
      process.oplsda <- try(opls(metrics, column_labels, orthoI = 9, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
    } 
    if("try-error" %in% class(process.oplsda)) {
      process.oplsda <- opls(metrics, column_labels, orthoI = 1, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none')}
    VIP <- getVipVn(process.oplsda)
    o <- order(-VIP)
    o <- o[1:10]
    o <- metrics[,o]
    return(list(o, as.vector(test.index)))
  }
  stopCluster(myCluster)
  fs
}

oplsda.fs.multiple_related <- function(adjacency_matrix, labels_, ppi_80, proteins_by_process, metric){
  results <- list()
  o.df <- NULL
  test.df <- NULL
  #for(n in 1:(length(unique(metrics[,1])))){
  for(n in 10:10){
    myCluster <- makeCluster(8, type='PSOCK')
    registerDoParallel(myCluster)
    if(metric=='bridge'){
      cat('Computing Bridge')
      cat('\n')
      ppis <- as.data.frame(ppi_80[n,,])
      protein_ppi_adj_matrix <- protein_ppi_by_process(adjacency_matrix, ppis, TRUE)
      protein_ppi <- protein_ppi_adj_matrix[[1]]
      adj_matrix <- protein_ppi_adj_matrix[[2]]
      proteins_red <- protein_ppi_adj_matrix[[3]]
      labels <- labels_[rownames(labels_) %in% proteins_red,]}
    if(metric=='scp'){
      cat('Computing SCP')
      cat('\n')
      ppis <- adjacency_matrix[adjacency_matrix$X == n-1,3:431]
      proteins_red <- adjacency_matrix[adjacency_matrix$X == n-1,2]
      labels <- labels_[rownames(labels_) %in% proteins_red,]}
    if(metric=='slb'){
      cat('Computing SLB')
      cat('\n')
      ppis <- as.data.frame(ppi_80[n,,])
      proteins_red <- graph_reduction(ppis)
      labels <- labels_[rownames(labels_) %in% proteins_red,]}
    
    cat('Initiating Foreach')
    cat('\n')
    fs <- foreach(column = 1:ncol(labels), .combine = 'rbind') %dopar% {
    #fs <- foreach(column = 1:2, .combine = 'rbind') %dopar% {
    #for(column in 1:2){
      gc()
      library(ropls)
      library(caret)
      library(reticulate)
      source_python('../python/lib/metrics_functions.py')
      column_labels <- factor(labels[,column], levels=c(1,0))
      process <- colnames(labels)[column]
      train.index <- unlist(createDataPartition(column_labels, p = .8, list = TRUE))
      all_labels  <- seq(1, length(column_labels))
      test.index <- setdiff(all_labels, train.index)
      if(metric=='bridge'){metrics <- bridge_scores(adj_matrix, protein_ppi, process)}
      if(metric=='scp'){metrics <- scp(ppis,process)}
      if(metric=='slb'){
        metrics <- slb(process, proteins_by_process, reduced=TRUE)}
      if(length(test.index) == 3440){
        test.index <- append(test.index, NA)
      }
      process.oplsda <- try(opls(metrics, column_labels, orthoI = NA, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
      if("try-error" %in% class(process.oplsda)) {
        process.oplsda <- try(opls(metrics, column_labels, orthoI = 9, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none'), silent=TRUE)
      } 
      if("try-error" %in% class(process.oplsda)) {
        process.oplsda <- opls(metrics, column_labels, orthoI = 1, predI=1, crossvalI = 10, subset=train.index, fig.pdfC = "none", info.txtC='none')}
      VIP <- getVipVn(process.oplsda)
      o <- order(-VIP)
      o <- o[1:10]
      o <- metrics[,o]
      cat(column)
      return(list(o, as.vector(test.index)))
    }
  fs.df <- as.data.frame(fs[,1])
  while(nrow(fs.df)<=16920){fs.df <- rbind(fs.df, rep(NA, 4290))}
  #while(nrow(fs.df)<=16920){fs.df <- rbind(fs.df, rep(NA, 20))}
  colnames(fs.df) <- seq(1,4290)
  #colnames(fs.df) <- seq(1,20)
  test.df <- as.data.frame(lapply(fs[,2], 'length<-', 17000))
  colnames(test.df) <- c(1,2,3,4,5,6,7,8,9,10)
  #colnames(test.df) <- c(1,2)
  if(n==10){
    np$savez_compressed('../python/data/processed/fs/reactome_slb_80_fs_apid_huri4.npz', fs.df ,allow_pickle = TRUE)
    write.csv(test.df,"../python/data/processed/fs/reactome_slb_80_test_apid_huri4.csv", row.names = FALSE)
  }
  else{
    fs.df_load <- np$load('../python/data/processed/fs/reactome_slb_80_fs_apid_huri3.npz')['arr_0']
    colnames(fs.df_load) <- seq(1,4290)
    np$savez_compressed('../python/data/processed/fs/reactome_slb_80_fs_apid_huri3.npz', rbind(fs.df_load, fs.df) ,allow_pickle = TRUE)
    test.df_load = read.csv('../python/data/processed/fs/reactome_slb_80_test_apid_huri3.csv')
    colnames(test.df_load) <- c(1,2,3,4,5,6,7,8,9,10)
    #colnames(test.df_load) <- c(1,2)
    write.csv(rbind(test.df_load, test.df),"../python/data/processed/fs/reactome_slb_80_test_apid_huri3.csv", row.names = FALSE)
    rm(fs.df_load)
    rm(test.df_load)
  }
  rm(fs.df)
  rm(test.df)
  gc()
  stopCluster(myCluster)
  }
  
}


####DATA LOAD####

#####PROCESS#####

adjacency_matrix <- np$load('../python/data/processed/adjacency_matrix_apid_huri.npy')
ppi_80 <-np$load('../python/data/processed/ppis/ppis_red80_apid_huri.npy', allow_pickle = TRUE)
reactome_labels <- read.csv('../python/data/processed/reactome_labels_apid_huri.csv', header=FALSE)
reactome_proteins_indexes <- read.csv('../python/data/processed/reactome_proteins_indexes_apid_huri.csv')
reactome_protein_ppi_by_process <- read.csv('../python/data/processed/protein_ppi_by_process_apid_huri.csv')

hypergeometric <- read.csv('../python/data/processed/metrics/process_hypergeometric_apid_huri.csv', row.names = 1)
closeness <- read.csv('../python/data/processed/metrics/process_closeness_apid_huri.csv', row.names = 1)
betweenness <- read.csv('../python/data/processed/metrics/process_betweenness_apid_huri.csv', row.names = 1)
fraction_betweenness <- read.csv('../python/data/processed/metrics/process_fraction_betweenness_apid_huri.csv', row.names = 1)
rwr <- read.csv('../python/data/processed/metrics/process_rwr_apid_huri.csv', row.names = 1)

row.names(reactome_labels) <- row.names(hypergeometric)
colnames(reactome_labels) <- colnames(hypergeometric)

hypergeometric_80 <- read.csv('../python/data/processed/metrics/process_ppi80_hyper_apid_huri.csv')
closeness_80 <- read.csv('../python/data/processed/metrics/process_ppi80_closeness_apid_huri.csv')
betweenness_80 <- read.csv('../python/data/processed/metrics/process_ppi80_betweenness_apid_huri.csv')
fraction_betweenness_80 <- read.csv('../python/data/processed/metrics/process_ppi80_fraction_betweenness_apid_huri.csv')
rwr_80 <- read.csv('../python/data/processed/metrics/process_ppi80_rwr_apid_huri.csv')
graph_apid_huri <- read_graph("../python/data/processed/graph_apid_huri", format = "gml")

hypergeometric_protein80 <- read.csv('../python/data/processed/metrics/process_protein80_hyper_apid_huri.csv')
closeness_protein80 <- read.csv('../python/data/processed/metrics/process_protein80_closeness_apid_huri.csv')
betweenness_protein80 <- read.csv('../python/data/processed/metrics/process_protein80_betweenness_apid_huri.csv')
rwr_protein80 <- read.csv('../python/data/processed/metrics/process_protein80_rwr_apid_huri.csv')
fraction_betweenness_protein80 <- read.csv('../python/data/processed/metrics/process_protein80_fraction_betweenness_apid_huri.csv')

#####DISEASE#####
######STEINER TREE######

disgenet_labels <- read.csv('../python/data/processed/disgenet_filtered_labels_apid_huri.csv', header=FALSE)
disgenet_proteins_indexes <- read.csv('../python/data/processed/disgenet_proteins_indexes_apid_huri_filtered.csv')

disease_hypergeometric <- read.csv('../python/data/processed/metrics/disease_hypergeometric_apid_huri.csv', row.names = 1)
disease_closeness <- read.csv('../python/data/processed/metrics/disease_closeness_apid_huri.csv', row.names = 1)
disease_betweenness <- read.csv('../python/data/processed/metrics/disease_betweenness_apid_huri.csv', row.names = 1)
disease_fraction_betweenness <- read.csv('../python/data/processed/metrics/disease_fraction_betweenness_apid_huri.csv', row.names = 1)
disease_rwr <- read.csv('../python/data/processed/metrics/disease_rwr_apid_huri.csv', row.names = 1)

row.names(disgenet_labels) <- row.names(disease_hypergeometric)
colnames(disgenet_labels) <- colnames(disease_hypergeometric)

disease_hypergeometric_protein80 <- read.csv('../python/data/processed/metrics/disease_protein80_hyper_apid_huri.csv')
disease_closeness_protein80 <- read.csv('../python/data/processed/metrics/disease_protein80_closeness_apid_huri.csv')
disease_betweenness_protein80 <- read.csv('../python/data/processed/metrics/disease_protein80_betweenness_apid_huri.csv')
disease_rwr_protein80 <- read.csv('../python/data/processed/metrics/disease_protein80_rwr_apid_huri.csv')
disease_fraction_betweenness_protein80 <- read.csv('../python/data/processed/metrics/disease_protein80_fraction_betweenness_apid_huri.csv')

disease_hypergeometric_ppi80 <- read.csv('../python/data/processed/metrics/disease_ppi80_hyper_apid_huri.csv')
disease_closeness_ppi80 <- read.csv('../python/data/processed/metrics/disease_ppi80_closeness_apid_huri.csv')
disease_betweenness_ppi80 <- read.csv('../python/data/processed/metrics/disease_ppi80_betweenness_apid_huri.csv')
disease_rwr_ppi80 <- read.csv('../python/data/processed/metrics/disease_ppi80_rwr_apid_huri.csv')
disease_fraction_betweenness_ppi80 <- read.csv('../python/data/processed/metrics/disease_ppi80_fraction_betweenness_apid_huri.csv')

######CONSERVATIVE#######

disgenet_labels_conservative <- read.csv('../python/data/processed/disgenet_conservative_labels_apid_huri.csv', header=FALSE)
disgenet_proteins_indexes_conservative <- read.csv('../python/data/processed/disgenet_prot_index_main_comp.csv')

disease_hypergeometric_conservative <- read.csv('../python/data/processed/metrics/disease_hypergeometric_conservative_apid_huri.csv', row.names = 1)
disease_closeness_conservative <- read.csv('../python/data/processed/metrics/disease_closeness_conservative_apid_huri.csv', row.names = 1)
disease_betweenness_conservative <- read.csv('../python/data/processed/metrics/disease_betweenness_conservative_apid_huri.csv', row.names = 1)
disease_fraction_betweenness_conservative <- read.csv('../python/data/processed/metrics/disease_fraction_betweenness_conservative_apid_huri.csv', row.names = 1)
disease_rwr_conservative <- read.csv('../python/data/processed/metrics/disease_rwr_conservative_apid_huri.csv', row.names = 1)

row.names(disgenet_labels_conservative) <- row.names(disease_hypergeometric_conservative)
colnames(disgenet_labels_conservative) <- colnames(disease_hypergeometric_conservative)

disease_hypergeometric_protein80_conservative <- read.csv('../python/data/processed/metrics/disease_protein80_hyper_conservative_apid_huri.csv')
disease_closeness_protein80_conservative <- read.csv('../python/data/processed/metrics/disease_protein80_closeness_conservative_apid_huri.csv')
disease_betweenness_protein80_conservative <- read.csv('../python/data/processed/metrics/disease_protein80_betweenness_conservative_apid_huri.csv')
disease_rwr_protein80_conservative <- read.csv('../python/data/processed/metrics/disease_protein80_rwr_conservative_apid_huri.csv')
disease_fraction_betweenness_protein80_conservative <- read.csv('../python/data/processed/metrics/disease_protein80_fraction_betweenness_conservative_apid_huri.csv')

disease_hypergeometric_ppi80_conservative <- read.csv('../python/data/processed/metrics/disease_ppi80_hyper_conservative_apid_huri.csv')
disease_closeness_ppi80_conservative <- read.csv('../python/data/processed/metrics/disease_ppi80_closeness_conservative_apid_huri.csv')
disease_betweenness_ppi80_conservative <- read.csv('../python/data/processed/metrics/disease_ppi80_betweenness_conservative_apid_huri.csv')
disease_rwr_ppi80_conservative <- read.csv('../python/data/processed/metrics/disease_ppi80_rwr_conservative_apid_huri.csv')
disease_fraction_betweenness_ppi80_conservative <- read.csv('../python/data/processed/metrics/disease_ppi80_fraction_betweenness_conservative_apid_huri.csv')

####FEATURE SELECTION PROCESS####

hyper_fs <- oplsda.fs(hypergeometric, reactome_labels)
write.csv(as.data.frame(hyper_fs[,1]),"../python/data/processed/fs/reactome_hyper_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(hyper_fs[,2]),"../python/data/processed/fs/reactome_hyper_test_apid_huri.csv", row.names = FALSE)

closeness_fs <- oplsda.fs(closeness, reactome_labels)
write.csv(as.data.frame(closeness_fs[,1]),"../python/data/processed/fs/reactome_closeness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(closeness_fs[,2]),"../python/data/processed/fs/reactome_closeness_test_apid_huri.csv", row.names = FALSE)

betweenness_fs <- oplsda.fs(betweenness, reactome_labels)
write.csv(as.data.frame(betweenness_fs[,1]),"../python/data/processed/fs/reactome_betweenness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(betweenness_fs[,2]),"../python/data/processed/fs/reactome_betweenness_test_apid_huri.csv", row.names = FALSE)

rwr_fs <- oplsda.fs(rwr, reactome_labels)
write.csv(as.data.frame(rwr_fs[,1]),"../python/data/processed/fs/reactome_rwr_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(rwr_fs[,2]),"../python/data/processed/fs/reactome_rwr_test_apid_huri.csv", row.names = FALSE)

fraction_betweenness_fs <- oplsda.fs(fraction_betweenness, reactome_labels)
write.csv(as.data.frame(fraction_betweenness_fs[,1]),"../python/data/processed/fs/reactome_fraction_betweenness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(fraction_betweenness_fs[,2]),"../python/data/processed/fs/reactome_fraction_betweenness_test_apid_huri.csv", row.names = FALSE)



hyper_80_fs <- oplsda.fs.multiple(hypergeometric_80, reactome_labels)
write.csv(as.data.frame(hyper_80_fs[[1]]),"../python/data/processed/fs/reactome_hyper_80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(hyper_80_fs[[2]]),"../python/data/processed/fs/reactome_hyper_80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

betweenness_80_fs <- oplsda.fs.multiple(betweenness_80, reactome_labels)
write.csv(as.data.frame(betweenness_80_fs[[1]]),"../python/data/processed/fs/reactome_betweenness_80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(betweenness_80_fs[[2]]),"../python/data/processed/fs/reactome_betweenness_80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

fraction_betweenness_80_fs <- oplsda.fs.multiple(fraction_betweenness_80, reactome_labels)
write.csv(as.data.frame(fraction_betweenness_80_fs[[1]]),"../python/data/processed/fs/reactome_fraction_betweenness_80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(fraction_betweenness_80_fs[[2]]),"../python/data/processed/fs/reactome_fraction_betweenness_80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

closeness_80_fs <- oplsda.fs.multiple(closeness_80, reactome_labels)
write.csv(as.data.frame(closeness_80_fs[[1]]),"../python/data/processed/fs/reactome_closeness_80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(closeness_80_fs[[2]]),"../python/data/processed/fs/reactome_closeness_80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

rwr_80_fs <- oplsda.fs.multiple(rwr_80, reactome_labels)
write.csv(as.data.frame(rwr_80_fs[[1]]),"../python/data/processed/fs/reactome_rwr_80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(rwr_80_fs[[2]]),"../python/data/processed/fs/reactome_rwr_80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)



hyper_protein80_fs <- oplsda.fs.multiple(hypergeometric_protein80, reactome_labels)
write.csv(as.data.frame(hyper_protein80_fs[[1]]),"../python/data/processed/fs/reactome_hyper_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(hyper_protein80_fs[[2]]),"../python/data/processed/fs/reactome_hyper_protein80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

betweenness_protein80_fs <- oplsda.fs.multiple(betweenness_protein80, reactome_labels)
write.csv(as.data.frame(betweenness_protein80_fs[[1]]),"../python/data/processed/fs/reactome_betweenness_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(betweenness_protein80_fs[[2]]),"../python/data/processed/fs/reactome_betweenness_protein80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

closeness_protein80_fs <- oplsda.fs.multiple(closeness_protein80, reactome_labels)
write.csv(as.data.frame(closeness_protein80_fs[[1]]),"../python/data/processed/fs/reactome_closeness_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(closeness_protein80_fs[[2]]),"../python/data/processed/fs/reactome_closeness_protein80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

rwr_protein80_fs <- oplsda.fs.multiple(rwr_protein80, reactome_labels)
write.csv(as.data.frame(rwr_protein80_fs[[1]]),"../python/data/processed/fs/reactome_rwr_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(rwr_protein80_fs[[2]]),"../python/data/processed/fs/reactome_rwr_protein80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)

fraction_betweenness_protein80_fs <- oplsda.fs.multiple(fraction_betweenness_protein80, reactome_labels)
write.csv(as.data.frame(fraction_betweenness_protein80_fs[[1]]),"../python/data/processed/fs/reactome_fraction_betweenness_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(fraction_betweenness_protein80_fs[[2]]),"../python/data/processed/fs/reactome_fraction_betweenness_protein80_test_apid_huri.csv", row.names = FALSE)
stopCluster(myCluster)


#bridge_fs <- oplsda_related.fs(adjacency_matrix,reactome_protein_ppi_by_process, reactome_labels)
#write.csv(as.data.frame(bridge_fs[1,]), "../python/data/processed/fs/reactome_bridge_fs_apid_huri.csv", row.names = FALSE)
#write.csv(as.data.frame(bridge_fs[2,]), "../python/data/processed/fs/reactome_bridge_test_apid_huri.csv", row.names = FALSE)

#scp_fs <- oplsda_related.fs(adjacency_matrix,closeness, reactome_labels, metric = 'scp')
#write.csv(as.data.frame(scp_fs[1,]), "../python/data/processed/fs/reactome_scp_fs_apid_huri.csv", row.names = FALSE)
#write.csv(as.data.frame(scp_fs[2,]), "../python/data/processed/fs/reactome_scp_test_apid_huri.csv", row.names = FALSE)

#slb_fs <-oplsda_related.fs(FALSE, FALSE, reactome_labels, metric = 'slb')
#write.csv(as.data.frame(slb_fs[1,]), "../python/data/processed/fs/reactome_slb_fs_apid_huri.csv", row.names = FALSE)
#write.csv(as.data.frame(slb_fs[2,]), "../python/data/processed/fs/reactome_slb_test_apid_huri.csv", row.names = FALSE)

#bridge_80_fs <- oplsda.fs.multiple_related(reactome_proteins_indexes, reactome_labels, ppi_80, metric='bridge')
#np$savez_compressed('../python/data/processed/fs/reactome_bridge_80_fs_apid_huri.npz', as.data.frame(bridge_80_fs[[1]]) ,allow_pickle = TRUE)
#write.csv(as.data.frame(bridge_80_fs[[2]]),"../python/data/processed/fs/reactome_bridge_80_test_apid_huri.csv", row.names = FALSE)

#scp_80_fs <- oplsda.fs.multiple_related(closeness_80, reactome_labels, FALSE, metric='scp')
#np$savez_compressed('../python/data/processed/fs/reactome_scp_80_fs_apid_huri.npz', as.data.frame(scp_80_fs[[1]]) ,allow_pickle = TRUE)
#write.csv(as.data.frame(scp_80_fs[[2]]),"../python/data/processed/fs/reactome_scp_80_test_apid_huri.csv", row.names = FALSE)

#slb_80_fs <- oplsda.fs.multiple_related(betweenness_80, reactome_labels, ppi_80, reactome_proteins_indexes, metric='slb')


####FEATURE SELECTION DISEASE####

#####STEINER TREE#####

disease_hyper_fs <- oplsda.fs(disease_hypergeometric, disgenet_labels)
write.csv(as.data.frame(disease_hyper_fs[,1]),"../python/data/processed/fs/disease/disease_hyper_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_fs[,2]),"../python/data/processed/fs/disease/disease_hyper_test_apid_huri.csv", row.names = FALSE)

disease_closeness_fs <- oplsda.fs(disease_closeness, disgenet_labels)
write.csv(as.data.frame(disease_closeness_fs[,1]),"../python/data/processed/fs/disease/disease_closeness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_fs[,2]),"../python/data/processed/fs/disease/disease_closeness_test_apid_huri.csv", row.names = FALSE)

disease_betweenness_fs <- oplsda.fs(disease_betweenness, disgenet_labels)
write.csv(as.data.frame(disease_betweenness_fs[,1]),"../python/data/processed/fs/disease/disease_betweenness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_fs[,2]),"../python/data/processed/fs/disease/disease_betweenness_test_apid_huri.csv", row.names = FALSE)

disease_fraction_betweenness_fs <- oplsda.fs(disease_fraction_betweenness, disgenet_labels)
write.csv(as.data.frame(disease_fraction_betweenness_fs[,1]),"../python/data/processed/fs/disease/disease_fraction_betweenness_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_fraction_betweenness_fs[,2]),"../python/data/processed/fs/disease/disease_fraction_betweenness_test_apid_huri.csv", row.names = FALSE)

disease_rwr_fs <- oplsda.fs(disease_rwr, disgenet_labels)
write.csv(as.data.frame(disease_rwr_fs[,1]),"../python/data/processed/fs/disease/disease_rwr_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_fs[,2]),"../python/data/processed/fs/disease/disease_rwr_test_apid_huri.csv", row.names = FALSE)



disease_hyper_ppi80_fs <- oplsda.fs.multiple(disease_hypergeometric_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_hyper_ppi80_fs[[1]]),"../python/data/processed/fs/disease/disgenet_hyper_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_ppi80_fs[[2]]),"../python/data/processed/fs/disease/disgenet_hyper_ppi80_test_apid_huri.csv", row.names = FALSE)

disease_betweenness_ppi80_fs <- oplsda.fs.multiple(disease_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_betweenness_ppi80_fs[[1]]),"../python/data/processed/fs/disease/disgenet_betweenness_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_ppi80_fs[[2]]),"../python/data/processed/fs/disease/disgenet_betweenness_ppi80_test_apid_huri.csv", row.names = FALSE)

disease_closeness_ppi80_fs <- oplsda.fs.multiple(disease_closeness_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_closeness_ppi80_fs[[1]]),"../python/data/processed/fs/disease/disgenet_closeness_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_ppi80_fs[[2]]),"../python/data/processed/fs/disease/disgenet_closeness_ppi80_test_apid_huri.csv", row.names = FALSE)

disease_rwr_ppi80_fs <- oplsda.fs.multiple(disease_rwr_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_rwr_ppi80_fs[[1]]),"../python/data/processed/fs/disease/disgenet_rwr_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_ppi80_fs[[2]]),"../python/data/processed/fs/disease/disgenet_rwr_ppi80_test_apid_huri.csv", row.names = FALSE)

disgenet_fraction_betweenness_ppi80_fs <- oplsda.fs.multiple(disease_fraction_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs[[1]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs[[2]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_test_apid_huri.csv", row.names = FALSE)



disease_hyper_protein80_fs <- oplsda.fs.multiple(disease_hypergeometric_protein80, disgenet_labels)
write.csv(as.data.frame(disease_hyper_protein80_fs[[1]]),"../python/data/fs/disease/disgenet_hyper_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_protein80_fs[[2]]),"../python/data/fs/disease/disgenet_hyper_protein80_test_apid_huri.csv", row.names = FALSE)

disease_betweenness_protein80_fs <- oplsda.fs.multiple(disease_betweenness_protein80, disgenet_labels)
write.csv(as.data.frame(disease_betweenness_protein80_fs[[1]]),"../python/data/fs/disease/disgenet_betweenness_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_protein80_fs[[2]]),"../python/data/fs/disease/disgenet_betweenness_protein80_test_apid_huri.csv", row.names = FALSE)

disease_closeness_protein80_fs <- oplsda.fs.multiple(disease_closeness_protein80, disgenet_labels)
write.csv(as.data.frame(disease_closeness_protein80_fs[[1]]),"../python/data/fs/disease/disgenet_closeness_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_protein80_fs[[2]]),"../python/data/fs/disease/disgenet_closeness_protein80_test_apid_huri.csv", row.names = FALSE)

disease_rwr_protein80_fs <- oplsda.fs.multiple(disease_rwr_protein80, disgenet_labels)
write.csv(as.data.frame(disease_rwr_protein80_fs[[1]]),"../python/data/fs/disease/disgenet_rwr_protein80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_protein80_fs[[2]]),"../python/data/fs/disease/disgenet_rwr_protein80_test_apid_huri.csv", row.names = FALSE)

disgenet_fraction_betweenness_ppi80_fs <- oplsda.fs.multiple(disease_fraction_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs[[1]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_fs_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs[[2]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_test_apid_huri.csv", row.names = FALSE)

#####CONSERVATIVE MODULE#####
disease_hyper_fs_conservative <- oplsda.fs(disease_hypergeometric_conservative, disgenet_labels_conservative)
write.csv(as.data.frame(disease_hyper_fs_conservative[,1]),"../python/data/processed/fs/disease/disease_hyper_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_fs_conservative[,2]),"../python/data/processed/fs/disease/disease_hyper_test_conservative_apid_huri.csv", row.names = FALSE)

disease_closeness_fs_conservative <- oplsda.fs(disease_closeness_conservative, disgenet_labels_conservative)
write.csv(as.data.frame(disease_closeness_fs_conservative[,1]),"../python/data/processed/fs/disease/disease_closeness_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_fs_conservative[,2]),"../python/data/processed/fs/disease/disease_closeness_test_conservative_apid_huri.csv", row.names = FALSE)

disease_betweenness_fs_conservative <- oplsda.fs(disease_betweenness_conservative, disgenet_labels_conservative)
write.csv(as.data.frame(disease_betweenness_fs_conservative[,1]),"../python/data/processed/fs/disease/disease_betweenness_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_fs_conservative[,2]),"../python/data/processed/fs/disease/disease_betweenness_test_conservative_apid_huri.csv", row.names = FALSE)

disease_fraction_betweenness_fs_conservative <- oplsda.fs(disease_fraction_betweenness_conservative, disgenet_labels_conservative)
write.csv(as.data.frame(disease_fraction_betweenness_fs_conservative[,1]),"../python/data/processed/fs/disease/disease_fraction_betweenness_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_fraction_betweenness_fs_conservative[,2]),"../python/data/processed/fs/disease/disease_fraction_betweenness_test_conservative_apid_huri.csv", row.names = FALSE)

disease_rwr_fs_conservative <- oplsda.fs(disease_rwr_conservative, disgenet_labels_conservative)
write.csv(as.data.frame(disease_rwr_fs_conservative[,1]),"../python/data/processed/fs/disease/disease_rwr_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_fs_conservative[,2]),"../python/data/processed/fs/disease/disease_rwr_test_conservative_apid_huri.csv", row.names = FALSE)



disease_hyper_ppi80_fs_conservative <- oplsda.fs.multiple(disease_hypergeometric_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_hyper_ppi80_fs_conservative[[1]]),"../python/data/processed/fs/disease/disgenet_hyper_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_ppi80_fs_conservative[[2]]),"../python/data/processed/fs/disease/disgenet_hyper_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_betweenness_ppi80_fs_conservative <- oplsda.fs.multiple(disease_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_betweenness_ppi80_fs_conservative[[1]]),"../python/data/processed/fs/disease/disgenet_betweenness_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_ppi80_fs_conservative[[2]]),"../python/data/processed/fs/disease/disgenet_betweenness_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_closeness_ppi80_fs_conservative <- oplsda.fs.multiple(disease_closeness_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_closeness_ppi80_fs_conservative[[1]]),"../python/data/processed/fs/disease/disgenet_closeness_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_ppi80_fs_conservative[[2]]),"../python/data/processed/fs/disease/disgenet_closeness_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_rwr_ppi80_fs_conservative <- oplsda.fs.multiple(disease_rwr_ppi80, disgenet_labels)
write.csv(as.data.frame(disease_rwr_ppi80_fs_conservative[[1]]),"../python/data/processed/fs/disease/disgenet_rwr_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_ppi80_fs_conservative[[2]]),"../python/data/processed/fs/disease/disgenet_rwr_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)

disgenet_fraction_betweenness_ppi80_fs_conservative <- oplsda.fs.multiple(disease_fraction_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)



disease_hyper_protein80_fs_conservative <- oplsda.fs.multiple(disease_hypergeometric_protein80, disgenet_labels)
write.csv(as.data.frame(disease_hyper_protein80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_hyper_protein80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_hyper_protein80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_hyper_protein80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_betweenness_protein80_fs_conservative <- oplsda.fs.multiple(disease_betweenness_protein80, disgenet_labels)
write.csv(as.data.frame(disease_betweenness_protein80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_betweenness_protein80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_betweenness_protein80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_betweenness_protein80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_closeness_protein80_fs_conservative <- oplsda.fs.multiple(disease_closeness_protein80, disgenet_labels)
write.csv(as.data.frame(disease_closeness_protein80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_closeness_protein80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_closeness_protein80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_closeness_protein80_test_conservative_apid_huri.csv", row.names = FALSE)

disease_rwr_protein80_fs_conservative <- oplsda.fs.multiple(disease_rwr_protein80, disgenet_labels)
write.csv(as.data.frame(disease_rwr_protein80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_rwr_protein80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disease_rwr_protein80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_rwr_protein80_test_conservative_apid_huri.csv", row.names = FALSE)

disgenet_fraction_betweenness_ppi80_fs_conservative <- oplsda.fs.multiple(disease_fraction_betweenness_ppi80, disgenet_labels)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs_conservative[[1]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_fs_conservative_apid_huri.csv", row.names = FALSE)
write.csv(as.data.frame(disgenet_fraction_betweenness_ppi80_fs_conservative[[2]]),"../python/data/fs/disease/disgenet_fraction_betweenness_ppi80_test_conservative_apid_huri.csv", row.names = FALSE)