#' Produces expected counts for each node for specfied DTM variables.
#'
#' @export
#' @param gamma Numeric matrix of Dirichlet parameter values.
#' @param tree An object of class `ape` specificying the DTM structure of Y.
#' @param N_i Optional numeric vector of sample sizes. If `NULL` (default), will return proportions.
#' @return A numeric `matrix`
#'
predict_DTM <- function(gamma, tree, N_i = NULL){
  # produces expected values for each NODE
  if(is.null(N_i)) N_i <- rep(1, nrow(gamma))
  V <- unique(tree$edge[,1])
  n <- nrow(gamma)

  y_VK <- matrix(nrow = n, ncol = max(tree$edge))
  y_VK[, V[1]] <- N_i
  for(vv in V){
    child_index <- tree$edge[which(tree$edge[,1] == vv),2]
    branches <- which(tree$edge[,1] == vv)
    N <- y_VK[, vv]
    shape <- gamma[, branches]
    if(length(branches) > 1){
      gamma_sum <- apply(shape, 1, sum)
      y_VK[, child_index] <- N * shape / gamma_sum
    }else{
      y_VK[, child_index]
    }
  }

  y_VK

}
