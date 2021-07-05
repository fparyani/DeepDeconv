
#' Quantile Normalization
#'
#' This function performs the quantile normalization on a gene expression matrix. Some warning of this function: May run into memory issue
#' on personal machine if matrix is too large. If needed subsample original matrix and continue operation.
#' 
#' @param gene_exp_mat This gene matrix is assumed to be sum normalized but otherwise a regular sparse matrix
#' @keywords gene selection
#' @return Returns a sparse matrix of the same dimension
#' @export
#' 
quant_norm <- function(gene_exp_mat){
  
  #Finds the 95th percentile for each gene across all the samples
  quant_95 <-  apply(gene_exp_mat,1,function(y) {quant <- quantile(y,probs = seq(0,1,0.05))[20] })
  
  #Transposes the given matrix and makes sure it is sparse
  gene_exp_mat_t <- Matrix::t(gene_exp_mat)
  gene_exp_mat_t <- as(gene_exp_mat_t, 'sparseMatrix')
  
  #Quant renormalization procedure#
  nz=diff(gene_exp_mat_t@p)
  quant_95=rep(0,length(nz))
  num_elements=floor(0.05*nrow(gene_exp_mat_t))
  nonzero_cols=which(nz>=num_elements)
  temp=apply(gene_exp_mat_t[,nonzero_cols],2,function(x){y=x[x>0];return(y[order(-y)[num_elements]])})
  quant_95[nonzero_cols]=temp
  
  #Assigns a value of 1 to values above 95th percentile
  gene_exp_mat_t@x[gene_exp_mat_t@x > quant_95[gene_exp_mat@i + 1] ] <- 1
  #Finds which values are below the 95th percentile
  gene_exp_true <- which(gene_exp_mat_t@x < quant_95[gene_exp_mat_t@i + 1])
  #Divides the values below the 95th percentile by the 95th percentile
  gene_exp_mat_t@x[gene_exp_true] <- gene_exp_mat_t@x[gene_exp_true] / quant_95[gene_exp_mat_t@i[gene_exp_true] + 1]
  
  gene_exp_mat_t <- Matrix::t(gene_exp_mat_t)
  
  return(gene_exp_mat_t)
  
}
