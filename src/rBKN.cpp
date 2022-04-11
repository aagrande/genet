#include "genet_types.h"
#include <stdio.h>

// [[Rcpp::export]]
void add_edges(MatrixXi &edge_list, R_xlen_t last_edge, int n_edges, MatrixXi &endpoints) {
  R_xlen_t i_edge = last_edge;
  for (int i = 0; i < n_edges; i++) {
    edge_list(i_edge, 0) = endpoints(i, 0);
    edge_list(i_edge, 1) = endpoints(i, 1);
    i_edge ++;
  }
}

// [[Rcpp::export]]
MatrixXi sample_endpoints(int i, int j, int n_edges, VectorXi &cards, VectorXi &leaders) {
  VectorXi left(n_edges), right(n_edges);
  left = Rcpp::as<VectorXi>(Rcpp::sample(cards(i), n_edges, true));
  if (i==j) {
    right = Rcpp::as<VectorXi>(Rcpp::sample(cards(j) - 1, n_edges, true));
    for (int i = 0; i < right.size(); i++) {
      if (left(i) <= right(i)) {
        right(i) = right(i) + 1;
      }
    }
  } else {
    right = Rcpp::as<VectorXi>(Rcpp::sample(cards(j), n_edges, true));
  }
  left = left.array() + leaders(i);
  right = right.array() + leaders(j);
  MatrixXi m(left.size(),2);
  m << left, right;
  return m;
}

// [[Rcpp::export]]
SEXP rBKN(
    VectorXi cards,
    const MSpMat H,
    const VectorXd root,
    bool verbose = false,
    bool way_too_verbose = false
) {
  int n_interactions = (H.cwiseSign().sum() + cards.size())/2;
  if (verbose) Rcpp::Rcout << "Number of interactions " << n_interactions << "\n";

  Rcpp::Function rpois("rpois"), sample("sample"), readline("readline");
  long long i_interaction = 0, cards_interaction;
  VectorXd expected_edges(n_interactions);
  for (int k=0; k<H.outerSize(); ++k){
    for (MSpMat::InnerIterator it(H,k); it; ++it)
    {
      if (it.row()>it.col()) continue;
      if (it.row()==it.col()) {
        long long int cards_i = cards(it.row()), cards_j = cards_i - 1;
        cards_interaction = cards_i * cards_j / 2;
      }
      if (it.row()<it.col()) {
        long long int cards_i = cards(it.row()), cards_j = cards(it.col());
        cards_interaction = cards_i * cards_j;
      }
      expected_edges(i_interaction) = cards_interaction * it.value() * root(it.row()) * root(it.col());
      i_interaction += 1;
    }
  }
  if (verbose) {
    // Rcpp::Rcout << "Max. number of connected pairs " << cards_interactions.sum() << "\n";
    Rcpp::Rcout << "Expected number of edges  " << expected_edges.sum() << "\n";
    Rcpp::Rcout << "Sampling edges for each interaction.. \n";
  }

  VectorXi edges(expected_edges.size());
  for (int i_interaction = 0; i_interaction < expected_edges.size(); i_interaction++) {
    edges(i_interaction) = Rcpp::as<int>(rpois(1, expected_edges(i_interaction)));
    // Rcpp::Rcout << "Edges for int. " << i_interaction << ": " << edges(i_interaction) << "\n";
  }
  double n_edges_total_double = 0;
  for (int i = 0; i < edges.size(); i++) {
    n_edges_total_double += edges(i);
  }
  long long int n_edges_total = n_edges_total_double;
  // int n_edges_total = edges.sum();
  if (verbose) Rcpp::Rcout << "We have sampled " << n_edges_total << " edges! \n";

  VectorXi leaders(cards.size()); // index of first node member of a type (0-based numbering)
  leaders(0) = 0;
  for (int i = 1; i < cards.size(); i++) {
    leaders(i) = leaders(i-1) + cards(i-1);
  }

  // readline("Initialize edge list?");
  MatrixXi edge_list(n_edges_total, 2);
  // readline("Edge list initialized!\n Generate edges?");
  R_xlen_t last_edge = 0;
  MatrixXi endpoints;
  i_interaction = 0;
  for (int k=0; k<H.outerSize(); ++k){
    for (MSpMat::InnerIterator it(H,k); it; ++it)
    {
      if (it.row()>it.col()) continue;
      if (edges(i_interaction)>0){
        endpoints = sample_endpoints(it.row(), it.col(), edges(i_interaction), cards, leaders);
        add_edges(edge_list, last_edge, edges(i_interaction), endpoints);
        // Rcpp::Rcout << last_edge << ", ";
        last_edge += edges(i_interaction);
        // if (last_edge > large_int) Rcpp::Rcout << last_edge << ", ";
        Rcpp::checkUserInterrupt();
      }
      i_interaction += 1;
    }
  }
  // readline("Edge generation is complete! \n");
  R_xlen_t n_edges = edge_list.rows();
  long long int max_rows = std::pow(2,31);
  if (n_edges < max_rows) {
    return Rcpp::wrap(edge_list);
  } else { // Rcpp::wrap doesn't work with long vectors
    readline("Initiate return?");
    Rcpp::IntegerVector edge_list_stacked(2 * n_edges);
    memcpy(INTEGER(edge_list_stacked), edge_list.array().data(), 2 * n_edges * sizeof(int));
    if (verbose) Rcpp::Rcout << " Matrix columns have been stacked into a single vector.\n";
    return edge_list_stacked;
  }
}
