test_that("generate_networks produces networks with target degree centrality", {
  K <- 10
  L <- 20
  K_each <- 3
  Z <- generate_balanced_assignments(K, L, K_each)
  cards <- rbinom(L, 100, .5) + 2
  average_degree <- 10
  target <- runif(L)
  target <- target * average_degree / (sum(cards * target) / sum(cards))

  g <- genet(1, Z, cards, dc = T, target)[[1]]
  testthat::expect_equal(igraph::vcount(g), sum(cards))

  empirical_type_degrees <- sapply(split(igraph::degree(g), rep(1:L, cards)),
                                   mean)
  plot(1:L, empirical_type_degrees, col = "red",
       xlab = "types", ylab = "degree", col.lab = "grey50", xaxt = "n")
  title(expression("Target " * phantom("vs empirical degree centrality")),
        col.main = "blue")
  title(expression(phantom("Target vs ") * "empirical " *
                     phantom("degree centrality")),
        col.main = "red")
  title(expression(phantom("Target ") * "vs " * phantom("empirical ") *
                     "degree centrality"),
        col.main = "grey50")
  segments(1:L, rep(0, L), 1:L, target, col = "blue")
  points(1:L, target, col = "blue", pch = 19, cex = .5)
  axis(2, col = "grey50", col.axis = "grey50")
  box(col = "grey50")
})
