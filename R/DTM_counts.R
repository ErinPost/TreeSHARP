#' Convert leaf counts to node or branch counts for DTM variable.
#'
#' @export
#' @param Y Numeric matrix of DTM observations
#' @param tree An object of class `ape` specificying the DTM structure of Y.
#' @param branches A logical (default = `FALSE`) indicating if counts should be in branch order.
#' @param prop A logical (default = `FALSE`) indicating if counts should be converted to proportions.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
DTM_counts <-  function(Y, tree, branches = FALSE, prop = FALSE){


  V <- unique(tree$edge[,1])
  n <- length(V)
  N <- nrow(Y)


  # output <- matrix(nrow = nrow(Y), ncol = V[n])
  # output[,c(1:ncol(Y))] <- Y


  output <- cbind(Y,  matrix(nrow = N,
                             ncol = V[n] - ncol(Y)))
  prop_output <- cbind(Y,  matrix(nrow = N,
                                  ncol = V[n] - ncol(Y)))
  names(output) <- NULL
  names(prop_output) <- NULL

  for(ii in n:1){

    child_index <- tree$edge[which(tree$edge[,1] == V[ii]),2]
    if(length(child_index) ==1){
      output[,V[ii]] <- output[,child_index]
    }else{
      output[,V[ii]] <- apply(output[,child_index], 1, sum)
    }
  }

  if(prop){
    for(ii in n:1){
      child_index <- tree$edge[which(tree$edge[,1] == V[ii]),2]
      if(length(child_index) ==1){
        prop_output[,child_index] <- 1
      }else{
        prop_output[,child_index] <- output[,child_index] / output[, V[ii]]
      }

    }
    prop_output<- ifelse(is.nan(prop_output), 0, prop_output)
    output <- prop_output
  }

  if(branches){
    # organize output by BRANCH index instead of NODE index
    output <- output[,tree$edge[,2]]
  }

  output
}

