#' Calculate likelihood for DTM random variable
#'
#' @export
#' @param y Numeric matrix of DTM count values.
#' @param gamma Numeric matrix of Dirichlet parameters for DTM variable
#' @param tree An object of class `ape` specificying the DTM structure of y.
#' @param log Logical indicating if log-likelihood should be returned (Default = `TRUE`)
#' @param simplify Logical. TRUE (default): return overall likelihood; FALSE: return observation-wise likelihoods
#' @return A numeric vector
#'
DTM_likelihood <- function(y, gamma, tree, log = TRUE, simplify = TRUE){

  if(is.data.frame(y)) y <- as.matrix(y)
  # y is a matrix... I x B -- need to modify the "actuals" to
  # gamma is also a matrix ... I x B
  V <- unique(tree$edge[,1])
  lprod_v <- matrix(nrow = nrow(y), ncol = length(V))
  index <- 1
  for(vv in V){
    child_index <- tree$edge[which(tree$edge[,1] == vv),2]
    branches <- which(tree$edge[,1] == vv)

    sum_cv <- lgamma(y[,branches] + gamma[, branches]) -
      (lgamma(y[,branches] + 1) + lgamma(gamma[,branches]))

    if(length(branches) == 1){
      lnum <- lgamma( y[,branches] +1 ) +
        lgamma(gamma[,branches])

      ldenom <- lgamma(y[,branches]+ gamma[,branches])

      lprod_v[,index] <- (lnum - ldenom) + sum_cv
    }else{
      lnum <- lgamma( apply(y[,branches],1,sum) +1 ) +
        lgamma(apply(gamma[,branches],1,sum))

      ldenom <- lgamma(apply(y[,branches],1,sum)
                       + apply(gamma[,branches],1,sum))

      lprod_v[,index] <- matrix((lnum - ldenom) + apply(sum_cv, 1, sum),
                                nrow = nrow(y))
    }

    index <- index + 1
  }

  out <- apply(lprod_v, 1, sum)

  #print(sum(log(out)))
  if(simplify){
    if(log){
      return(sum(out))
    }else{
      return(exp(sum(out)))
    }}else{
      if(log){
        return(out)
      }else{
        return(exp(out))
      }
    }

}
