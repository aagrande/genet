#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

typedef MappedSparseMatrix<double> MSpMat;
typedef SparseMatrix<double> SpMat;
typedef Eigen::Matrix<int64_t,Dynamic,1> VectorXi64;
