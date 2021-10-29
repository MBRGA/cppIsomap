#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::VectorXi;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Ref;

// Calculates pairwise Euclidean distances between rows of a matrix
MatrixXd cppPairwiseDistances(Map<MatrixXd> data);
VectorXi cppOrder(const Ref<const VectorXd> sortvec);