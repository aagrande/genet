test_that("test generate_balanced_assignments", {
  K <- 100
  L <- 200
  K_each <- 5
  Z <- generate_balanced_assignments(K, L, K_each)

  testthat::expect_equal(nrow(Z), L)
  testthat::expect_equal(ncol(Z), K)
  testthat::expect_equal(Matrix::rowSums(Z > 0), rep(K_each, L))

  H <- Z %*% Matrix::t(Z)
  is_connected <- igraph::is.connected(igraph::graph_from_adjacency_matrix(
    as.matrix(H>0)))
  testthat::expect_identical(is_connected, TRUE)
})
