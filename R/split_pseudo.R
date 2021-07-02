#' Split Pseudo Data
#'
#' This function prepares the training and testing data for training of the pseudo spots. The training and testing data are portions of gene
#' expression matrix and the labels are the portions normalized of the cell type of interest. 
#' 
#' @param reduced_mat This takes in a traditional gene expression matrix that should have been normalized and reduced from "important_gene" function
#' @param  count_metadata This refers to the cell type count found in metadata and should include all cell types. Should have same cell numbers as reduced_mat
#' @param coi This is the cell of interest, takes in the name of the cell type you are deconvoluting. Should be name from metadata
#' @param ratio This is the percentage split of training/test. Default to 70/30 training/testing data. Input range from 0 - 1
#' @keywords training
#' @return Returns a list of training/testing data and its respective labels
#' @export
#' 
split_pseudo <- function(reduced_mat, count_metadata, coi, ratio = .70){
  
  df_sim <- count_metadata
  sim_cell_count <- df_sim %>% as.matrix()
  sim_labels <- apply(sim_cell_count,1,function(x) {  x/(x %>% sum)} ) %>% t()
  cell_num <- which(str_locate(colnames(sim_cell_count), coi) == 1)
  
  #Assigning training and testing data
  sample <- sample( dim(reduced_mat)[2], size = floor(ratio*dim(reduced_mat)[2]), replace = F)
  
  train_labels <- sim_labels[,cell_num][sample]
  test_labels <- sim_labels[,cell_num][-sample]
  
  train_cell_data <- reduced_mat[,sample] %>% as.matrix() %>% t()
  test_cell_data <- reduced_mat[,-sample] %>% as.matrix() %>% t()
  
  #Randomization
  shuffle_train <- sample( c(1:dim(train_cell_data)[1]) ,dim(train_cell_data)[1] )
  shuffle_test <- sample( c(1:dim(test_cell_data)[1]) ,dim(test_cell_data)[1] )
  
  train_cell_data <- train_cell_data[shuffle_train,]
  test_cell_data <- test_cell_data[shuffle_test,]
  
  train_labels <- train_labels[shuffle_train]
  test_labels <- test_labels[shuffle_test]
  
  return(list("training data" = train_cell_data, "training label" = train_labels, "test data" = test_cell_data, "test label" = test_labels))
  
  
}
