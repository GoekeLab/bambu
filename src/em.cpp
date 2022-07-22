#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;


//' EM algorithm
//' @noRd
// [[Rcpp::export]]
List em_theta (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const int maxiter,
               const double minvalue,
               const double conv // , const int nThr = 1

) {
  int M = X.n_rows; //number of isoforms

  // omp_set_num_threads(nThr) ; // using multiple threads

  // containers
  arma::rowvec theta(M) ; // define theta
  theta.fill(1) ;

  arma::mat theta_trace(M,maxiter);
  int iter = 0; // iterator
  theta_trace.col(iter) = theta.t();

  // initialize Xb

  arma::mat Xb = X ;
  arma::mat lmat = X;


  arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(),0);

  arma::mat adjP_Xt = X;   // adjusted p-values = sampling probability X theta X b, with bias parameter include
  arma::rowvec summed_by_col_adjP_Xt = arma::sum(adjP_Xt,0);   // will be used in logLikelihood estimation as well as probMat update

  double deltaTheta = 1;
  arma::vec t_before = theta_trace.col(iter);
  arma::vec t_after = theta_trace.col(iter);
  arma::vec deltaVec(M);
  while(deltaTheta > conv && iter < (maxiter-1)){

    //update iterator
    iter++;
    adjP_Xt =  (X.t() * diagmat(theta)).t();
    summed_by_col_adjP_Xt = arma::sum(adjP_Xt,0);
    lmat = (adjP_Xt * diagmat(Y / summed_by_col_adjP_Xt));
    lmat.replace(arma::datum::nan, 0);
    theta = arma::sum(lmat.t(),0)/summed_by_row_Xb; //sum by row
    //theta.replace(arma::datum::nan, 0);
    theta_trace.col(iter) = theta.t();
    //Rcout << "The value of 1-est : " << theta << "\n";
    t_before = theta_trace.col(iter-1);
    t_after = theta_trace.col(iter);
    deltaVec = abs(t_after.elem( find(t_after > minvalue)) -
          t_before.elem( find(t_after > minvalue)))/t_after.elem( find(t_after > minvalue));
    if(deltaVec.is_empty()){
          deltaTheta = 0;
    }
    if(!deltaVec.is_empty()){
         deltaTheta = max(deltaVec);
    }
  }
  // returns
  List ret ;
  ret["theta"] = theta ;
  ret["theta_trace"] = theta_trace ;
  ret["iter"] = iter;
  return(ret) ;
}


//' L1-penalized likelihood estimation
//' @noRd
// [[Rcpp::export]]
List emWithL1 (const arma::cube A, // alignment compatibility matrix array, last dimension 1 is all, 2 is full, 3 is partial, 4 is unique
               const arma::rowvec Y, // observed number of reads for each read class j
               const double K, //total read count
               const double lambda,  // the tuning parameter for bias estimation, take as
               const int maxiter,
               const double minvalue,
               const double conv  // , const int nThr = 1
){

  // initialization
  arma::mat X = A.slice(0); 
  int M = X.n_rows; //number of isoforms


  List theta_out(3); // create a empty list of size 3
  arma::rowvec theta(M);
  
  theta_out = em_theta(X, Y, lambda, maxiter, minvalue, conv) ;
  theta = Rcpp::as<arma::rowvec>(theta_out[0]) ;

  // post-process outputs
  arma::mat estMat(4,M);
  estMat.row(0) = theta;
  arma::rowvec baseSum = Y / arma::sum((X.t()*diagmat(theta)).t(),0);
  baseSum.replace(arma::datum::nan, 0);
  estMat.row(1) = arma::sum(((X.t()*diagmat(theta)).t() * diagmat(baseSum)).t(), 0) * K;
  estMat.row(2) = arma::sum(((A.slice(1).t()*diagmat(theta)).t() * diagmat(baseSum)).t(), 0) * K;
  estMat.row(3) = arma::sum(((A.slice(2).t()*diagmat(theta)).t() * diagmat(baseSum)).t(), 0) * K;
  // returns
  List ret ;
  ret["theta"] = estMat;

  return(ret) ;


}




