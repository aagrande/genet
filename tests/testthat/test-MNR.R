test_that("MNR recovers a positive root of the degree centrality problem", {
  compute_system <- function(theta, weighted_interactions, d) {
    L <- length(theta)
    f <- Matrix::rowSums(Matrix::t(weighted_interactions) *
                   matrix(theta, nrow = L, ncol = L, byrow = T))
    f <-  f * theta - d
    return(f)
  }

  K <- 100
  L <- 200
  K_each <- 5
  Z <- generate_balanced_assignments(K, L, K_each)
  H <- Z %*% Matrix::t(Z)
  theta_ini <- runif(L)
  target <- runif(L)
  cards <- rbinom(L, 100, .5) + 2

  cards_interactions <- compute_weighted_interactions(H, cards, dc = T, target)
  root <- MNR(cards_interactions, theta_ini, target, verbose = F)
  value_at_root <- compute_system(root, cards_interactions, target)

  testthat::expect_identical(all(root > 0), TRUE)
  testthat::expect_equal(value_at_root, rep(0, L))
})
