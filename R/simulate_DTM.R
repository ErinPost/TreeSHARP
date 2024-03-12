#' Simulate a DTM for use with treesharp_stan function.
#'
#' @export
#' @param subject_sim number of subjects/groups/observations to simulate
#' @param subject_N vector of length(subject_sim with group sizes. NULL (default) will randomly sample from Unif(7000, 10000)
#' @param tree object of type `ape`
#' @param num_leaf if tree is NULL, number of leaves on generated tree
#' @param covariates_sim number of covariates to generate (i.e. ncol(X))
#' @param which_continuous indexes of continuous covariates
#' @param which_binary indexes of binary covariates
#' @param binary_p success probability for binary draws
#' @param zeta matrix indicating nonzero variables for each branch
#' @param rho correlation for continuous covariates at lag 1. Either rho or Sigma must be specified.
#' @param Sigma correlation matrix for continuous covariates. Either rho or Sigma must be specified.
#' @param num_cov number of nonzero regression coefficients per branch
#' @param beta_min minimum absolute value for regression coefficients
#' @param beta_max maximum absolute value for regression coefficients
#' @param seed random seed
#' @return A list with, X, Y, alpha, beta, zeta, tree
#'
simulate_DTM <- function (subject_sim = 100,
                          subject_N = NULL,
                          tree = NULL,
                          num_leaf = 5,
                          covariates_sim = 50,
                          which_coninuous = NULL,
                          which_binary = NULL,
                          binary_p = NULL,
                          zeta = NULL,
                          rho = 0.3,
                          Sigma = NULL,
                          num_cov = 5,
                          beta_min = 0.8,
                          beta_max = 1.2,
                          seed = 1234)
{

  ############################################################################
  # This code has been adapted from the function simulate_DTM found in the   #
  # MicroBVS package created by Matthew Koslovsky.                           #
  # *    Title: MicroBVS::simulate_DTM() source code                         #
  # *    Author: Koslovsky, Matthew                                          #
  # *    Date: 2023                                                          #
  # *    Availability: https://github.com/mkoslovsky/MicroBVS                #
  ############################################################################

  library(mvtnorm)
  library(MCMCpack)
  library(ape)

  set.seed(seed)

  if (!is.null(rho)) {
    if(is.null(which_coninuous)){
      Sigma <- matrix(0, covariates_sim, covariates_sim)
      Sigma = rho^abs(row(Sigma) - col(Sigma))
    }else{
      Sigma <- matrix(0, length(which_coninuous),
                      length(which_coninuous))
      Sigma = rho^abs(row(Sigma) - col(Sigma))
    }
  }
  tree.ex <- if (is.null(tree)) {
    rtree(n = num_leaf)
  }else{
    tree
  }
  V <- tree.ex$Nnode
  Cv <- table(tree.ex$edge[, 1])
  K <- length(tree.ex$tip.label)
  B_sim <- sum(Cv)

  # generate the X matrix as a mixture of continuous and binary

  if(!is.null(which_coninuous) & !is.null(which_binary)){
    X <- matrix(0, nrow = subject_sim, ncol = covariates_sim)
    X[,which_coninuous] <- scale(rmvnorm(subject_sim,
                                         rep(0, length(which_coninuous)), Sigma))

    X[,which_binary] <- scale(matrix(sample(c(1,0),
                                            length(which_coninuous) * subject_sim,
                                            replace = TRUE,
                                            prob = c(binary_p, 1 - binary_p))),
                              center = TRUE,scale=TRUE)

  }else{
    X <- scale(rmvnorm(subject_sim, rep(0, covariates_sim), Sigma))
  }

  if(!is.null(zeta)){

    zeta_sim <- zeta

  }else{
    zeta_sim <- matrix(0, B_sim, covariates_sim)
    for (i in sample(seq(1, B_sim), B_sim)) {
      select <- sample(1:covariates_sim, num_cov)
      zeta_sim[i, select] <- 1
    }
    # for (i in sample(seq(1, B_sim), num_branch)) {
    #     select <- sample(1:covariates_sim, num_cov)
    #     zeta_sim[i, select] <- 1
    # }
  }

  alpha_sim <- matrix(1, nrow = subject_sim, ncol = 1) %*%
    (runif(n = B_sim, -1.3, 1.3))
  true_cov <- which(zeta_sim == 1)

  beta_sim <- matrix(0, B_sim, covariates_sim)
  beta_sim[true_cov] <- runif(sum(zeta_sim), beta_min, beta_max) *
    sample(c(-1, 1), sum(zeta_sim), replace = TRUE)
  inside_sim <- exp(alpha_sim + X %*% t(beta_sim))

  node_counts <- matrix(0, nrow = subject_sim, ncol = (V + K))
  if(!is.null(subject_N)){
    node_counts[, (K + 1)] <- subject_N
  }else{
    node_counts[, (K + 1)] <- sample(seq(7500, 10000), subject_sim)
  }
  for (b in (K + 1):(sum(Cv) + 1)) {

    node <- which(tree.ex$edge[, 1] == b)
    inside_branches <- inside_sim[, node]
    prob_sim <- apply(inside_branches, 1, function(x) {
      out <-rdirichlet(1, x)
      if(any(is.na(out))){
        out <- x/sum(x)
      }else{
        out
      }
    })
    for (i in 1:subject_sim) {
      y <- t(rmultinom(1, node_counts[i, b], t(prob_sim)[i,
      ]))
      node_counts[i, (tree.ex$edge[node, 2])] <- y
    }
  }
  Y <- node_counts[, 1:K]

  return(list(Y = Y, X = X,
              alpha_sim = alpha_sim,
              beta_sim = beta_sim,
              zeta_sim = zeta_sim,
              tree = tree.ex,
              Sigma = Sigma))
}

