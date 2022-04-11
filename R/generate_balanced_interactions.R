#' Generate balanced interactions
#'
#' @param K number of communities
#' @param L number of types
#' @param K_each number of communities each type should be assigned to
#' @param max_ite number of maximum attempts to generate assignments
#'
#' @return A L-by-K sparse matrix of type community assignments.
#' @export
generate_balanced_assignments = function(K, L, K_each = NA, max_ite = 1000) {
  if (is.na(K_each)) {
    K_each <- min(K - 1, round(max(K / L, L / K )))
  }
  assignments <- matrix(0, nrow = L, ncol = K)
  ite <- 1
  is_connected <- FALSE
  while (any(Matrix::colSums(assignments)==0) | !is_connected |
         anyDuplicated.matrix(assignments, MARGIN = 1))
  {
    assignments <- Matrix::sparseMatrix(i = rep(1:L, rep(K_each, L)),
                                        j = sapply(rep(K, L), sample, size = K_each),
                                        x = 1,
                                        dims = c(L, K))
    H <- assignments %*% Matrix::t(assignments)
    is_connected <- igraph::is.connected(igraph::graph_from_adjacency_matrix(
      as.matrix(H>0)))
    if (ite == max_ite)
      stop("Maximum iterations reached. No assignments were generated.")
    else ite <- ite + 1
  }
  assignments <- Matrix::Diagonal(x = apply(assignments, 1,
                                            function(x) sqrt(sum(x^2))^(-1))) %*%
    assignments
  return(assignments)
}
