#This takes in snRNA data that has already been quantile normalized

#Takes in a factored group and the quantile_normalized matrix



#' Find Important Gene
#'
#' The purpose of this function is to reduce the dimensions of a gene expression matrix by finding the most 
#' relevant genes of a particular cell type using pairwise Wilcox test. This function assumes the distribution of your cell types in quant_mat
#' reflects the set of genes you are searching for. 
#' 
#' @param quant_mat A gene expression matrix that has already been quantile normalized
#' @param groups The groups entered refer to the variouscell types of your dataset and are assumed to be in factored form when entered. 
#' @param cell_type This is the name of the cell whose gene signature you are looking for, should be one of the name from "groups"
#' @param st_gene If you are working with spatial transcriptomic data, add its gene list to ensure feasibility when running the model
#' @param num_gene Hyper-parameter that approximately determines number of genes to sample from all permutation of pairwise groups
#' @keywords gene selection
#' @return Returns a vector of the location of the gene on the matrix inputted
#' @export
#' 
important_gene <- function(quant_mat, factor_group, cell_type, st_gene = NA, num_gene = 500 ){
  
  pwt <- scran::pairwiseWilcox(quant_mat, groups=factor_group, direction = 'up')
  result_pwt <- stats::setNames(lapply(pwt$statistics,  function(x) {as.data.table(x, keep.rownames=TRUE) }), apply(pwt$pairs, 1, paste, collapse = '_'))     
  gene_up_exp <- lapply(result_pwt,function(x){ x$p.value %>% order() } )
  important_gene <- lapply(gene_up_exp,function(x){ x[1:num_gene] })
  
  coi = which(pwt$pairs$first == cell_type)
  
  cell_gene <- important_gene[coi] %>% unlist() %>% unique() %>% sample()
  
  if (is.na(st_gene) == TRUE){
    
    return (cell_gene)
    
  } else {
    
    final_gene_s <- match(rownames(quant_mat)[cell_gene],st_gene)
    final_gene_s <- final_gene_s[!is.na(final_gene_s)]
    cell_gene <- match(st_gene[final_gene_s],rownames(quant_mat))
    
    return(cell_gene)
    
  }
  
  
}
