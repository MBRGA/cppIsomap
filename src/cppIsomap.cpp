#include "cppIsomap.h"

//' Computes the Isomap embedding as introduced in 2000 by Tenenbaum, de Silva and Langford.
//' 
//' Isomap is a nonlinear dimension reduction technique, that preserves global properties of the data. 
//' That means, that geodesic distances between all samples are captured best in the low dimensional embedding. 
//' This C++ version is based R version by Christoph Bartenhagen, which is in turn based on the Matlab implementation by Tenenbaum and uses Floyd's Algorithm to compute the neighbourhood graph of shortest distances, when calculating the geodesic distances.
//' A modified version of the original Isomap algorithm is included. It respects nearest and farthest neighbours. 
//' 
//' @param data N x D matrix (N samples, D features)
//' @param dims vector containing the target space dimension(s)
//' @param k number of neighbours
//' @param mod use modified Isomap algorithm
//' @param verbose show a summary of the embedding procedure at the end
//' 
//' @return It returns a N x dim matrix (N samples, dim features) with the reduced input data (list of several matrices if more than one dimension was specified)
//'  
//' @examples
//' \dontrun{
//' ## two dimensional Isomap embedding of a 1.000 dimensional dataset using k=5 neighbours
//' ## two dimensional Isomap embedding of a 1.000 dimensional dataset using k=5 neighbours
//' d = RDRToolbox::generateData(samples=20, genes=1000, diffgenes=100, blocksize=10)
//' d_low = Isomap(data=d[[1]], dims=2, k=5)
//' ## Isomap residuals for target dimensions 1-10
//' d_low = Isomap(data=d[[1]], dims=1:10, k=5, plotResiduals=TRUE)	
//'     
//' ## three dimensional Isomap embedding of a 1.000 dimensional dataset using k=10 (nearest and farthest) neighbours
//' d = generateData(samples=20, genes=1000, diffgenes=100, blocksize=10)
//' d_low = Isomap(data=d[[1]], dims=3, mod=TRUE, k=10)
//' }
//' @seealso \link[RDRToolbox]{Isomap}
//' @name cppIsomap
//' @export
// [[Rcpp::export]]
List cppIsomap(NumericMatrix data, IntegerVector dims = IntegerVector::create(2), Nullable<IntegerVector> k = R_NilValue, bool mod = false, bool verbose = true) {
  
  for (auto& x: dims) {
    if (x < 1) stop("invalid argument: target dimension is required to be a (vector of) positive integer value(s)");
  }
  
  auto num_samples = data.nrow();
  auto num_features = data.ncol();
  
  int k_ = -1;
  
  if (k.isNull() || k.as().size() < 1) {
    k_ = num_samples > 5 ? 5 : num_samples;
    Rcout << "neighbours argument missing; setting to " << k_ << "\n";
  } else if (k.as().size() > 1) {
    k_ = k.as()(0);
    Rcout << "multiple values given for neighbours argument; only using first\n";
  } else {
    k_ = k.as()(0);
  }
  
  if (k_ >= num_samples)
    stop("invalid argument: more neighbours than samples");
  if (k_ <= 1) 
    stop("invalid argument: neighbour parameter is required to be an integer value >= 2");
  
  Rcout << "Computing distance matrix ... ";
  
  Map<MatrixXd> dataEigen(data.begin(), num_samples, num_features);
  
  MatrixXd d = cppPairwiseDistances(dataEigen);
  
  Rcout << "done\n";
  
  int first_k, last_k;
  
  if (mod) {
    first_k = round((double)k_ / 2.0);
    last_k = k_ - first_k;
    Rcout << "Building graph with shortest paths (using " << first_k << " nearest and " << last_k << " farthest neighbours) ... ";
  } else {
    first_k = k_;
    last_k = 0;
    Rcout << "Building graph with shortest paths (using " << first_k << " nearest neighbours) ... ";
  }

  for (auto i = 0; i < num_samples; i++) {
    // original Isomap gets column ordering within column, but given matrix is symmetric doing it by row works 
    // and doesn't need us to save a whole matrix of orderings in advance
    VectorXi rowsort = cppOrder(d.row(i)); 
    for (auto j = first_k + 1; j < num_samples - last_k; j++) 
      d(i, rowsort[j]) = std::numeric_limits<double>::infinity();
  }
  
  int infoint = 1000;
  
  // d = pmin(d, t(d))
  for (auto i = 0; i < num_samples; i++) {
    for (auto j = i + 1; j < num_samples; j++) {
      if (d(i, j) < d(j, i)) {
        d(j, i) = d(i, j);
      } else if (d(j, i) < d(i, j)) {
        d(i, j) = d(j, i);
      }
    }
  }
  
  for (auto i = 0; i < num_samples; i++) {
    if (--infoint == 0) {
      Rcout << "\n" << i << " samples processed ... ";
      infoint = 1000;
      checkUserInterrupt();
    }
    
    VectorXd samplecol = d.col(i);
    
    for (auto j = 0; j < num_samples; j++) {
      for (auto k = j + 1; k < num_samples; k++) {
        auto newd = samplecol(j) + samplecol(k);
        d(j, k) = (d(j, k) < newd ? d(j, k) : newd);
        d(k, j) = (d(k, j) < newd ? d(k, j) : newd);
      }
    }
  }
  
  std::vector<int> num_connections(num_samples);
  std::vector<int> first_connections(num_samples);
  std::unordered_set<int> components_set = {};
  
  for (auto i = 0; i < num_samples; i++) {
    auto fc = -1;
    for (auto j = 0; j < num_samples; j++) {
      if (d(i, j) != std::numeric_limits<double>::infinity()) {
        num_connections[i] += 1;
        if (fc == -1) fc = j;
      }
    }
    first_connections[i] = fc;
    components_set.insert(fc);
  }
  
  // We can't access the unordered set nicely so copy to a vector
  std::vector<int> components(components_set.begin(), components_set.end());
  
  auto num_components = components.size();
  
  Rcout << "done\n";
  Rcout << "Computing low dimensional embedding ... ";
  
  List all_Y = List::create();
  
  for (auto& dim: dims) {
    NumericMatrix Yout(num_samples, dim);
    Map<MatrixXd> Y(Yout.begin(), num_samples, dim);
    
    for(auto c = 0; c < num_components; c++) {
      
      checkUserInterrupt();
      
      std::vector<int> comp_indices = {};
      
      for (auto i = 0; i < first_connections.size(); i++)
        if(first_connections[i] == components[c]) comp_indices.push_back(i);
      
      auto N = comp_indices.size();
      MatrixXd D(N, N);
      
      for(auto i = 0; i < N; i++)
        for(auto j = 0; j < N; j++)
          D(i,j) = d(comp_indices[i], comp_indices[j]);
      
      MatrixXd T;
      VectorXd D2rowsums;
      
      VectorXd ones;
      ones = VectorXd::Ones(N);
      
      D2rowsums = D.array().square().matrix().rowwise().sum();
      
      T = (- 0.5 * (D.array().square().matrix() - D2rowsums * ones.transpose() / (double)N 
                     - ones * D2rowsums.transpose() / (double)N).array() -0.5 * D.array().square().sum() / (double)(N * N)).matrix();
      
      if (T.rows() < dim)
        stop("Problem occured while computing embedding: Connected component %i"
               " consists of too few samples (%i) and cannot be reduced to %i dimensions."
               " Try Isomap with a neighbourhood parameter higher than %i"
               " or with dimension parameters lower than %i.", c, T.rows(), dim, k, dim);
      
      SelfAdjointEigenSolver<MatrixXd> es(T);
      
      MatrixXd eigvecs = es.eigenvectors();
      MatrixXcd eigvalsqrt = es.eigenvalues()
        .transpose()
        .cast<std::complex<double>>()
        .array()
        .sqrt(); // transpose because more useful to us as row vector. Need complex because of negative eigenvalues
      
      // Note we reverse eigenvalues and eigenvectors
      // Eigen returns eigenvalues in increasing order, R gives DECREASING order
      for (auto i = 0; i < N; i++)
        Y.row(comp_indices[i]) = (eigvecs.rightCols(dim).row(i).reverse().array() * 
          eigvalsqrt.rightCols(dim).reverse().array()).matrix().real();
      
      all_Y.push_back(Y, "dim" + std::to_string(dim));
    }
  }
  
  Rcout << "done\n";
  if (verbose) {
    Rcout << "number of samples: " << num_samples << "\n";
    Rcout << "reduction from " << num_features << " to " << dims << " dimensions\n";
    Rcout << "number of connected components in graph: " << num_components << "\n";
  }
  
  return all_Y;
  
}

VectorXi cppOrder(const Ref<const VectorXd> sortvec) {
  
  VectorXi retval(sortvec.size());
  
  for(auto i = 0; i < sortvec.size(); i++)
    retval(i) = i;
  
  std::sort(retval.data(),
            retval.data() + retval.size(),
            [&sortvec](long a, long b) {
              return sortvec[a] < sortvec[b];
            });
  
  return retval;
}

MatrixXd cppPairwiseDistances(Map<MatrixXd> data) {
  auto num_samples = data.rows();
  auto num_features = data.cols();
  
  MatrixXd pd(num_samples, num_samples);
  //NumericMatrix pd(num_samples, num_samples);
  
  for (auto i = 0; i < num_samples; i++) {
    
    for (auto j = i + 1; j < num_samples; j++) {
      
      double sqd = 0;
      
      for (auto k = 0; k < num_features; k++) {
        sqd += (data(i,k) - data(j,k)) * (data(i,k) - data(j,k));
      }
      
      pd(i, j) = pd(j, i) = sqrt(sqd);
    }
    
    pd(i,i) = 0; // We want diagonal to be 0
  }
  
  return pd;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

cppIsomap(as.matrix(wdbc[1:10,-1]), dims = c(3,6), k=3)

#randmat = matrix(runif(12),nrow=3, ncol=4)

#RDRToolbox:::pairwiseDistances(randmat)

#cppPairwiseDistances(randmat)

#cppIsomap(matrix(0, nrow = 3, ncol = 3), k = 1)

*/
