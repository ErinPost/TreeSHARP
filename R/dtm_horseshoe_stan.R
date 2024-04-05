#' Bayesian Dirichlet Tree Multinomial regression with Stan
#'
#' @export
#' @param X Numeric matrix of input values.
#' @param Y Numeric matrix of output values.
#' @param tree An object of class `ape` specificying the DTM structure of Y.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
dtm_horseshoe_stan <- function(X, Y, tree,...){

  # Prepare stan inputs
  Y_branches <- DTM_counts(Y = Y,
                           tree = tree,
                           branches = FALSE,
                           prop = FALSE)

  edge_inputs <- treesharp_stan_tree_inputs(tree)

  stan_inputs <- list(
    B = nrow(tree$edge), # number of branches
    I = nrow(X),
    Y = Y_branches,
    P = ncol(X),
    X = X,
    V = length(edge_inputs$V),
    parents = array(edge_inputs$V),
    children = array(edge_inputs$C_v),
    n_children = array(edge_inputs$n_C_v),
    parent_rep = array(edge_inputs$V_rep),
    branches = array(edge_inputs$B_v)
  )

  # Execute stan fit

  out <- rstan::sampling(stanmodels$dtm_horseshoe, data = stan_inputs, ...)

  # Return results
  return(out)

}
