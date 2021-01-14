#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;


//' EM algorithm
//' @noRd
// [[Rcpp::export]]
List em_theta (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const arma::rowvec b, // bias parameter treated as fixed while estimating theta
               const bool d, // indicator of whether bias parameter will be used
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
  if(d){
    Xb = X * diagmat(exp(b));
  }

  arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(),0);

  //Xb = (Xb.t() * diagmat(1.0/summed_by_row_Xb)).t();
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
List emWithL1 (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const bool d, // indicator of whether bias parameter will be used
               const int maxiter,
               const double minvalue,
               const double conv  // , const int nThr = 1
){

  // initialization
  int J = X.n_cols; //number of equivalent read class
  int M = X.n_rows; //number of isoforms

  arma::rowvec est(M+J); // bias parameter
  est.fill(1);

  arma::mat est_trace(M+J,maxiter);


  int iter = 0; // iterator
  est_trace.col(iter) = est.t();

  List theta_out(3); // create a empty list of size 5
  arma::rowvec theta = est.head(M);
  //theta_trace.col(iter) = theta.t();

  List b_out(5); // create a empty list of size 5
  arma::rowvec b(J);
  b.fill(1);

  if(d){


    arma::rowvec summed_by_col_adjP_Xt = arma::sum(X,0);
    arma::rowvec signvec(J);

    arma::rowvec signinputvec(J);
    arma::rowvec abssignvec(J);
    arma::rowvec maxVec(J);
    est.tail(J) = b;

    int maxiter_b = 500;

    double deltaTheta = 1;
    arma::vec e_before = est_trace.col(iter);
    arma::vec e_after = est_trace.col(iter);
    arma::vec deltaVec(M);
    while(deltaTheta > conv && iter < (maxiter_b-1)){

      //update iterator
      iter++;
      theta = est.head(M) ;
      b = est.tail(J) ;

      // at each estimation, run EM with b fixed
      theta_out = em_theta(X, Y, lambda,b, d, 2, minvalue, conv) ;
      theta = Rcpp::as<arma::rowvec>(theta_out[0]) ;
      est.head(M) = theta ;

      // at each estimation, update b with new theta
      summed_by_col_adjP_Xt = arma::sum((X.t() * diagmat(theta)).t(),0);
      signinputvec = Y - summed_by_col_adjP_Xt;
      abssignvec = abs(signinputvec)-lambda;
      maxVec = (abssignvec % (abssignvec > 0));
      signvec = diagvec(arma::sign(diagmat(signinputvec))).t() % maxVec;
      //abssignvec.elem( find(abssignvec < 0) ).zeros();
      b = diagvec(arma::log1p(diagmat(signvec / summed_by_col_adjP_Xt))).t();
      b = b - median(b);
      b.replace(arma::datum::nan, 0);
      est.tail(J) = b;
      //Rcout << "The value of 2-b : " << b << "\n";
      est_trace.col(iter) = est.t();
      //Rcout << "The value of 2-est : " << theta << "\n";
      e_before = est_trace.col(iter-1).head(M);
      e_after = est_trace.col(iter).head(M);
      deltaVec = abs(e_after.elem( find(e_after > minvalue)) - 
          e_before.elem( find(e_after > minvalue)))/e_after.elem( find(e_after > minvalue));
      if(deltaVec.is_empty()){
          deltaTheta = 0; 
      }
      if(!deltaVec.is_empty()){
          deltaTheta = max(deltaVec);
      }

    }
  }
  if(!d){
    theta_out = em_theta(X, Y, lambda,b, d, maxiter, minvalue, conv) ;
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




