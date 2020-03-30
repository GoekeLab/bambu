#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;


//' EM algorithm
//'
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param conv convergence threshold of likelihoods
//' @export
// [[Rcpp::export]]
List em_theta (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const arma::rowvec b, // bias parameter treated as fixed while estimating theta
               const bool d, // indicator of whether bias parameter will be used
               const int maxiter,
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
  if(d){
    Xb = X * diagmat(exp(b));
  }


  arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(),0);

  Xb = (Xb.t() * diagmat(1.0/summed_by_row_Xb)).t();

  arma::mat adjP_Xb = Xb;   // adjusted p-values = sampling probability X theta X b, with bias parameter include
  arma::rowvec summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0);   // will be used in logLikelihood estimation as well as probMat update

  double deltaTheta = 1;
  while(deltaTheta > conv && iter < (maxiter-1)){

    //update iterator
    iter++;
    adjP_Xb =  (Xb.t() * diagmat(theta)).t();
    summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0);
    lmat = (adjP_Xb * diagmat(Y / summed_by_col_adjP_Xb));
    lmat.replace(arma::datum::nan, 0);
    theta = arma::sum(lmat.t(),0); //sum by row
    theta_trace.col(iter) = theta.t();

    deltaTheta = dot(theta_trace.col(iter) - theta_trace.col(iter-1), theta_trace.col(iter) - theta_trace.col(iter-1));


  }

  // returns
  List ret ;
  ret["theta"] = theta ;
  ret["theta_trace"] = theta_trace ;
  ret["iter"] = iter;
  return(ret) ;
}




//' L1-penalized likelihood estimation
//'
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param conv convergence threshold of likelihoods
//' @export
// [[Rcpp::export]]
List emWithL1 (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const bool d, // indicator of whether bias parameter will be used
               const int maxiter,
               const double conv  // , const int nThr = 1
){

  // initialization
  int J = X.n_cols; //number of equivalent read class
  int M = X.n_rows; //number of isoforms


  arma::rowvec est(M+J); // bias parameter
  est.fill(1);

  arma::mat est_trace(M+J,maxiter);
  //arma::mat theta_trace(M,maxiter);



  int iter = 0; // iterator
  est_trace.col(iter) = est.t();

  List theta_out(3); // create a empty list of size 5
  arma::rowvec theta = est.head(M);
  //theta_trace.col(iter) = theta.t();

  List b_out(5); // create a empty list of size 5
  arma::rowvec b(J);
  b.fill(1);


  if(d){

    arma::mat Xb = X * arma::diagmat(exp(b));
    arma::rowvec summed_by_col_Xb = arma::sum(Xb,0);
    arma::rowvec log_column_sum_Xb(J);


    arma::rowvec summed_by_col_adjP_X = arma::sum(X,0);
    arma::rowvec signvec(J);

    arma::rowvec signinputvec(J);
    b = b - diagvec(arma::log(diagmat(summed_by_col_Xb))).t();
    est.tail(J) = b;

    int maxiter_b = 500;

    double deltaTheta = 1;

    while(deltaTheta > conv && iter < (maxiter_b-1)){

      //update iterator
      iter++;
      theta = est.head(M) ;
      b = est.tail(J) ;


      // at each estimation, run EM with b fixed
      theta_out = em_theta(X, Y, lambda,b, d, maxiter, conv) ;
      theta = Rcpp::as<arma::rowvec>(theta_out[0]) ;
      est.head(M) = theta ;


      // at each estimation, update b with new theta
      summed_by_col_adjP_X = arma::sum((X.t() * diagmat(theta)).t(),0);
      signinputvec = Y- summed_by_col_adjP_X;
      signvec = diagvec(arma::sign(diagmat(signinputvec))).t() % ((abs(signinputvec)-lambda) % ((abs(signinputvec)-lambda) > 0));
      b = diagvec(arma::log1p(diagmat(signvec / summed_by_col_adjP_X))).t();
      b = b - median(b);
      Xb = X * diagmat(exp(b));
      summed_by_col_Xb = arma::sum(Xb,0);
      log_column_sum_Xb = arma::diagvec(arma::log(arma::diagmat(summed_by_col_Xb))).t();
      b = b - log_column_sum_Xb;


      est.tail(J) = b ;

      est_trace.col(iter) = est.t();
      deltaTheta = dot(est_trace.col(iter) - est_trace.col(iter-1), est_trace.col(iter) - est_trace.col(iter-1));

    }
  }
  if(!d){

    theta_out = em_theta(X, Y, lambda,b, d, maxiter, conv) ;
    theta = Rcpp::as<arma::rowvec>(theta_out[0]) ;
    est.head(M) = theta ;
    est.tail(J).zeros();
  }


  // returns
  List ret ;
  ret["theta"] = est.head(M);
  ret["b"] = est.tail(J) ;
  ret["nzindex"] = find(est.tail(J) > 0) ;
  return(ret) ;


}


