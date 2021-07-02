
#' Train Neural Network
#'
#' The purpose of this function is to traing the pseudo spot data on a neural network. This function uses an external depedency for setting
#' up the neural network: Keras. Can freely change parameters of neural network. Feel free to refer to Keras guide on NN.
#' 
#' @param train_cell_data Transpose of regular gene expression matrix of data that requires training. split_pseudo does the transposing automatically
#' @param train_labels This is the training data which is the portion of cell type to predict.
#' @keywords training
#' @return Returns a trained network. Can be saved as an .h5 file using TensorFlow
#' @export
#' 
train_network <- function(train_cell_data, train_labels){
  
 network <- keras::keras_model_sequential() 
  network %>% 
    #layer_dropout(rate = 0.5) %>% 
    # layer_dense(units = 128, activation = "relu", input_shape = dim(train_cell_data)[2] ) %>%
    layer_dense(units = 64, activation = "relu", input_shape = dim(train_cell_data)[2] ) %>%
    # layer_dropout(rate = 0.1) %>%
    layer_dense(units = 32, activation = "relu" ) %>%
    # layer_dense(units = 4, activation = "relu") %>%
    layer_dropout(rate = 0.1)    %>%
    layer_dense(units = 1, activation = "relu")
  
  network %>% keras::compile(
    optimizer = 'adam',
    # optimizer = optimizer_adam(lr = 0.8),
    # optimizer = optimizer_rmsprop(),
    # optimizer = 'sgd',
    # loss = "hinge",
    loss = "mae",
    # metrics = list("mean_absolute_error")
    metrics = c("accuracy")
  )
  
  
  history <- network %>% keras::fit(
    train_cell_data, train_labels, 
    epochs = 15, batch_size = 50, 
    validation_split = 0.2
  )
  
  return(network)
  
}

