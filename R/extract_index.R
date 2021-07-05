#' @export
#' 
#' 


extract_index <- function(cell_data){
  
  index <- lapply(
    cell_data,
    function(x)
    {
      extract_row <- rownames(x) %>% as.numeric()
    }
  )
  
  return(index)
}