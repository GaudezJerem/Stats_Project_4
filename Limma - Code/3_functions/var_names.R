
var_names <- function(colData, pattern, begin, end) {

  
  
  if(begin == FALSE){
    start <- 1
  } else {
   
    start <- regexpr(pattern = pattern, 
                     text = colData)+begin %>%
      unlist()
    
  }
  
  stop <- regexpr(pattern = pattern, 
                  text = colData)+end %>%
    unlist()
  
  
  substr(colData, start, stop)
  
}








