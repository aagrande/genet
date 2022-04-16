rBKN_wrapper <- function(n, cards, H, root, verbose = T, way_too_verbose = F,
                         parallel = F) {
  if (!methods::is(H, "dgCMatrix")) stop("H is not a dgCMatrix")
  if (parallel & !requireNamespace("parallel", quietly = TRUE)) {
    warning("Package \"parallel\" must be installed to generate networks
              parallelly. The generator will run on a single-thread.")
  }
  null_list <- vector("list", n)
  if (parallel & requireNamespace("parallel", quietly = TRUE)) {
    graphs <- parallel::mclapply(null_list,
                                 function(x) rBKN(cards, H, root, verbose = F,
                                                  way_too_verbose = F))
  }
  if (!parallel) {
    graphs <- lapply(null_list,
                     function(x) rBKN(cards, H, root, verbose, way_too_verbose))
  }
  # note that the vertices max(x) + 1, ..., sum(cards) may be missing
  n_vertices_to_add <- sum(cards) - sapply(graphs, function(x) max(x))
  graphs <- lapply(graphs, igraph::graph_from_edgelist, directed = FALSE)
  if (any(n_vertices_to_add > 0)) {
    for (g in which(n_vertices_to_add > 0)) {
      graphs[[g]] <- igraph::add_vertices(graphs[[g]], nv = n_vertices_to_add[g])
    }
  }
  return(graphs)
}

compute_expected_degrees <- function(root, interactions, target_eigenc = NA) {
  if (is.na(target_eigenc[1])) {
    return(as.vector(root * (root %*% interactions)))
  } else {
    interactions <- diag(1 / target_eigenc) %*% interactions
    return(as.vector(root * (root %*% interactions)))
  }
}

#' Generate Networks
#'
#' Generate networks with overlapping communities and target node centralities.
#'
#' @param n number of random networks to generate.
#' @param multiedges logical. If \code{FALSE}, omit multiedges.
#' @param assignments matrix with non-negative entries and whose rows encode the
#'  membership assignments of types.
#' @param cards vector of type cardinalities. Must be whole numbers.
#' @param dc logical. If \code{TRUE}, generate networks with target degree
#'   centrality; if \code{FALSE}, generate networks with target eigencentrality.
#' @param target_centrality vector of target centralities.
#' @param target_degree target average expected degree. This argument is
#'   required only when generating networks with target eigencentrality (i.e.
#'   \code{dc = FALSE}).
#' @param verbose logical. If \code{TRUE}, give verbose output.
#' @param return_root logical. Should the degree correction terms be returned?
#'
#' @details The rows of the matrix argument `assignments` are normalized so as
#' to have unit norm.
#'
#' @return A list with \code{n} \code{igraph} graph objects. If 
#' \code{return_root = TRUE}, a list with the following two objects:
#' \item{graphs}{a list with \code{n} \code{igraph} graph objects;}
#' \item{root}{a vector with the degree correction terms.}
#' 
#' @examples
#'assignments <- Matrix::sparseMatrix(i = c(1, 2, 3, 3, 3),
#' j = c(1, 2, 1, 2, 3),
#' x = c(1, 1, 1, 2, 1))
#' row_norms <- apply(assignments, 1, function(x) sqrt(sum(x ^ 2)))
#' assignments <- Matrix::diag(row_norms ^ -1) %*% assignments
#' 
#' cards <- c(4, 4, 6)
#' 
#' target_degree <- c(4, 6, 10)
#' 
#' g <- genet(n = 4, multiedges = F, assignments = assignments, cards = cards, dc = T,
#'            target_centrality = target_degree)
#'            
#' par(mfrow = c(2, 2), mar = rep(0, 4))
#' invisible(lapply(X = g, FUN = igraph::plot.igraph, vertex.size = 4,
#'                  vertex.label = NA, vertex.color = "blue", edge.arrow.mode = 0))
#'
#' @export
genet <- function(n, multiedges = T, assignments, cards, dc = TRUE, 
                  target_centrality, target_degree = NA, verbose = FALSE,
                  return_root = FALSE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  L <- length(cards)
  assignments <- methods::as(assignments, "sparseMatrix")
  assignments <- Matrix::Diagonal(x = apply(assignments, 1,
                                            function(x) sqrt(sum(x^2))^(-1))) %*%
    assignments
  H <- assignments %*% Matrix::t(assignments)
  interactions <- compute_weighted_interactions(H, cards, dc, target_centrality)
  theta_ini <- initialize_root(interactions, target_centrality) # runif(L)
  root <- MNR(interactions, theta_ini, target_centrality, verbose = verbose)
  if (!dc) { # rescale root to comply with target average degree
    if (is.na(target_degree[1])) {
      stop("In the eigencentrality problem, the target average degree should be
           specified.")
    }
    avg_expected_degree <-
      sum(cards * compute_expected_degrees(root, interactions, 
                                           target_centrality)) / sum(cards)
    root <- root * sqrt(target_degree / avg_expected_degree)
  }
  out <- rBKN_wrapper(n, cards, H, root, verbose, F)
  if (!multiedges) out <- lapply(out, igraph::simplify)
  if (return_root) {
    out <- list(graphs = out, root = root)
  }
  return(out)
}
