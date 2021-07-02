

#' Sum Normalization
#'
#' This function normalizes the columns of the gene expression matrix, so that each column sums to one. In other words all 
#' the genes in each samples sums to one
#' 
#' @param sparse_mat Takes in a sparse gene expression matrix
#' @keywords normalization
#' @return Returns a sparse matrix
#' @export
#' 
sum_norm <- function(sparse_mat){
  
  gene_exp_mat <- sparse_mat
  gene_exp_mat@x = gene_exp_mat@x/rep.int(Matrix::colSums(gene_exp_mat),diff(gene_exp_mat@p))
  
  return(gene_exp_mat)

  
}


