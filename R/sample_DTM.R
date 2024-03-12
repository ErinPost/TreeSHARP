#' Sample from a DTM distribution
#'
#' @export
#' @param gamma Numeric matrix of dirichlet parameters for each branch.
#' @param N Numeric vector of sample sizes for each observation.
#' @param tree An object of class `ape` specificying the DTM structure of Y.
#' @return A numeric `matrix` object with samples'
sample_DTM <- function(gamma, N, tree){

  nodes <- max(tree$edge)
  V <- unique(tree$edge[,1])
  I <- length(N)

  samples <- matrix(nrow = I, ncol = nodes)

  samples[,V[1]] <- N # root node has fixed count

  for(vv in V){
    #print(vv)
    child_index <- tree$edge[which(tree$edge[,1] == vv),2]
    branches <- which(tree$edge[,1] == vv)

    N_vv <- samples[, vv]
    shape <- gamma[, branches]

    input <- if(nrow(gamma) == 1){
      matrix(c(shape, N_vv),nrow=1)
    }else{
      cbind(shape, N_vv)
    }
    hold <- apply(input, 1,
                  function(x){
                    n <- length(x)
                    p <- NA
                    while(all(is.na(p))){
                      suppressWarnings(p <- MCMCpack::rdirichlet(1,
                                                                 x[1:(n-1)]))                         }
                    stats::rmultinom(n = 1, size = x[n], prob = p)
                  } )

    samples[, child_index] <- t(hold)

  }
  return(samples)
}
