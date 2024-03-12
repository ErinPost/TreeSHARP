#' Support function for treesharp_stan
#'
#' @export
#' @param tree An object of class `ape` specificying the DTM structure of Y.
#' @return A `list` with 6 element
#'
treesharp_stan_tree_inputs <- function(tree){
  df <- tree$edge |>
    data.frame() |>
    group_by(X1) |>
    dplyr::summarize(children = n()) |>
    dplyr::select(node = X1, children)

  edge_order <- tree$edge |>
    data.frame() |>
    mutate(branch = 1:n()) |>
    arrange(X1)

  out <- list(B = nrow(tree$edge),
       C_v = edge_order$X2,
       B_v = edge_order$branch,
       V = df$node,
       n_C_v = df$children,
       V_rep = rep((1:length(df$node)), df$children))

  return(out)
}
