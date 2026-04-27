#include <Rcpp.h>

//' @title Convert node-block indices into groups of node indices
//' @description Converts a vector of block indices into a list of
//' node index vectors.
//' @param block An integer vector of 1-based block indices, one for each node.
//' @param max_block Precomputed maximal value of `block`. Must be at least
//' as large as `max(block)`. Block indices larger than `max_block` are ignored.
//' @param per_node Logical; if `FALSE`, the function returns a list with
//' one element for each block. If `TRUE`, the function returns a list with
//' one element for each node.
//' @returns If `per_node` is `false`, a list of integer vectors, where each
//' vector contains the 1-based node indices of nodes belonging to the
//' corresponding block. If `per_node` is `true`, a list of integer vectors,
//' where each vector contains the 1-based node indices of nodes belonging to
//' the same block as the corresponding node; `b[i] <- a[block[j]]`
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List inlabru_group_cv_block_conversion(Rcpp::IntegerVector block,
                                             int max_block,
                                             bool per_node) {
  int n = block.size();
  std::vector< std::vector<int> > node_vectors(max_block);

  for (int i = 0; i < n; i++) {
    int b = block[i];
    if (b >= 1 && b <= max_block) {
      node_vectors[b - 1].push_back(i + 1);
    }
  }

  // Convert to an R List
  Rcpp::List block_list(max_block);
  for (int j = 0; j < max_block; j++) {
    block_list[j] = Rcpp::wrap(node_vectors[j]);
  }
  if (per_node) {
    Rcpp::List node_list(n);
    for (int i = 0; i < n; i++) {
      node_list[i] = block_list[block[i] - 1];
    }
    return node_list;
  }
  return block_list;
}
