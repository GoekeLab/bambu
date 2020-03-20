#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;



//' Create likelihood objective function to estimate bias parameters
//'
//' @param par A row vector of bias parameters for each read class
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param theta estimates for isoform expression
double obj_fun_b_rcpp(const arma::rowvec par,
                         const arma::mat X,
                         const arma::rowvec Y,
                         const double lambda,
                         const arma::rowvec theta
){


  arma::mat a_adj_mat = X * diagmat(exp(par)); // define tmp variables

  arma::rowvec a_adj_mat_rowSum = arma::sum(a_adj_mat.t(),0);

  a_adj_mat = (a_adj_mat.t() * diagmat(a_adj_mat_rowSum % theta)).t();

  arma::rowvec column_agg = arma::sum(a_adj_mat,0);
  arma::rowvec log_column_agg = diagvec(arma::log(diagmat(column_agg))).t();
  log_column_agg.replace(arma::datum::nan, 0); // replace each NaN with 0
  // Compute objective value
  double obj_val = -(sum(Y%log_column_agg-column_agg)-lambda*sum(abs(par)));

  // Return a single value
  return obj_val;
}







//' Optimization function for bias parameters
//'
//' @param par A row vector of bias parameters for each read class
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param theta estimates for isoform expression
// [[Rcpp::export]]
List optim_b_rcpp(const arma::rowvec theta,
                const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
                const arma::rowvec Y, // observed number of reads for each read class j
                const double lambda,  // the tuning parameter for bias estimation, take as
                const arma::rowvec b){

  // Extract Rs optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = b,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_b_rcpp),
                                 Rcpp::_["method"] = "SANN", // objective function is non-differentiable
                                 Rcpp::_["hessian"] = "FALSE",
                                 Rcpp::_["X"] = X,
                                 Rcpp::_["Y"] = Y,
                                 Rcpp::_["lambda"] = lambda,
                                 Rcpp::_["theta"] = theta);


  // Return estimated values
  return opt_results;
}





//' EM algorithm with L1-penalty
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
  theta.fill(1) ; // initialize thetas to 1/M
  arma::mat probMat = X ; // define prob matrix
  probMat.fill(1) ;


  arma::mat theta_trace(M,maxiter);
  int iter = 0; // iterator
  theta_trace.col(iter) = theta.t();

  // initialize Xb

  arma::mat Xb = X ;
  if(d){
    Xb = X * diagmat(exp(b));
  }

  arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(),0);

  // normalized Xb
  // for(int i = 0 ; i < M ; i++){
  //   Xb.row(i) = Xb.row(i) / summed_by_row_Xb(i);
  // }
  // arma:: mat Xb_tr = Xb.t();
  arma::mat Xb_tr = Xb.t() * diagmat(summed_by_row_Xb);
  Xb = Xb_tr.t();

  //Xb = X


  arma::mat adjP_Xb = Xb; // adjusted p-values = sampling probability X theta X b, with bias parameter include
  arma::rowvec summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0); // will be used in logLikelihood estimation as well as probMat update

  double deltaTheta = 1;
  while(deltaTheta > conv && iter < (maxiter-1)){

    //update iterator
    iter++;


    adjP_Xb =  (Xb.t() * diagmat(theta)).t();


    //
    summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0);

    probMat = adjP_Xb * diagmat(Y / summed_by_col_adjP_Xb);
    probMat.replace(arma::datum::nan, 0);

    theta = arma::sum(probMat.t(),0);

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




//' EM algorithm with L1-penalty
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



  int iter = 0; // iterator
  est_trace.col(iter) = est.t();
  List theta_out(2); // create a empty list of size 5
  arma::rowvec theta = est.head(M);



  List b_out(5); // create a empty list of size 5
  arma::rowvec b = est.tail(J);




  double deltaEst = 1;

  // algorithm
  while(deltaEst > conv && iter < (maxiter-1)){

    //update iterator
    iter++;

    theta = est.head(M) ;
    b = est.tail(J) ;


    // at each estimation, run EM one time
    theta_out = em_theta(X, Y, lambda,b, d, maxiter, conv) ;
    theta = Rcpp::as<arma::rowvec>(theta_out[0]) ;




    est.head(M) = theta ;
    if(d){
      b_out = optim_b_rcpp(theta,X, Y, lambda,b) ;
      b = Rcpp::as<arma::rowvec>(b_out[0]) ;
      // to shink b so that it won't change theta estimates significantly
      if(b.size()>2){
        b = (b - mean(b))/max(abs(b)); // in case very different bs will produce very unstable estimates, normalized to N(0,0.1)
      }
      if(b.size()<=2){
        b = (b - mean(b))/max(abs(b));
      }
      b.elem(find(b <= 0)).zeros();
      est.tail(J) = b ;
    }


    est_trace.col(iter) = est.t();
    deltaEst = dot(est_trace.col(iter) - est_trace.col(iter-1), est_trace.col(iter) - est_trace.col(iter-1));

  }
  Rcout << "iter : " << iter << "\n" ;
  Rcout << "theta : " << theta << "\n" ;
  Rcout << "b : " << b << "\n" ;


  // returns
  List ret ;
  ret["theta"] = est.head(M);
  ret["b"] = est.tail(J) ;
  ret["nzindex"] = find(est.tail(J) > 0) ;
  return(ret) ;


}


//' Create likelihood objective function to estimate bias parameters
//'
//' @param par A row vector of bias parameters for each read class
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param nzindex estimates for isoform expression
double obj_fun_rcpp(const arma::rowvec par,
                    const arma::mat X,
                    const arma::rowvec Y,
                    const double lambda,
                    const arma::rowvec nzindex
){
  int J = X.n_cols; //number of equivalent read class
  int M = X.n_rows; //number of isoforms
  int p = nzindex.size() ; // number of non-zero b estimates
  // predefine adjusted sampling matrix to be the sampling matrix
  arma::mat a_adj_mat = X;

  // of the parameters inputed, the first M are theta estimates
  // the rest are bias estimates
  arma::rowvec theta = par.head(M);
  arma::rowvec b = par.tail(p);

  // for read class with non-zero bias estimates, adjust sampling matrix by bias parameter
  for(int i = 0 ; i < p ; i++){
    a_adj_mat.col(nzindex(i)) = X.col(nzindex(i)) * exp(b(i));
  }

  // obtain the rowsum of the adjusted sampling matrix
  arma::rowvec a_adj_mat_rowSum = arma::sum(a_adj_mat.t(),0);

  // for each row, divide the adjusted sampling matrix by rowsum and then multiply by theta estimates
  for(int i = 0 ; i < M ; i++){
    a_adj_mat.row(i) = a_adj_mat.row(i) / a_adj_mat_rowSum(i) * theta(i);
  }

  // calculate the column sum of the result obtained above
  arma::rowvec column_agg = arma::sum(a_adj_mat,0);

  // define log of column_agg
  arma::rowvec log_column_agg(X.n_cols);

  // for each obtain the log value of the column_agg
  for(int i = 0 ; i < J ; i++){
    log_column_agg(i) = std::log(column_agg(i));
  }


  //log_column_agg.replace(arma::datum::nan, 0); // replace each NaN with 0

  double obj_val = -(sum(Y%log_column_agg-column_agg)); //define objective value

  //Rcout << "obj_fun : 7 \n";
  // Return a single value
  return obj_val;
}

//' Optimization function for bias parameters
//'
//' @param par A row vector of bias parameters for each read class
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param theta estimates for isoform expression
// [[Rcpp::export]]
List optim_thetaandb_rcpp(const arma::rowvec est,
                          const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
                          const arma::rowvec Y, // observed number of reads for each read class j
                          const double lambda,  // the tuning parameter for bias estimation, take as
                          const arma::rowvec nzindex){

  // Extract Rs optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  // arma::uvec lowervalue(est.size());
  // lowervalue.fill(10^(-8));
  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = est,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                 Rcpp::_["method"] = "SANN", // objective function is non-differentiable
                                 Rcpp::_["hessian"] = "FALSE",
                                 Rcpp::_["X"] = X,
                                 Rcpp::_["Y"] = Y,
                                 Rcpp::_["lambda"] = lambda,
                                 Rcpp::_["nzindex"] = nzindex);


  // Return estimated values
  return opt_results;
}

//' EM algorithm with L1-penalty
//'
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param conv convergence threshold of likelihoods
//' @export
// [[Rcpp::export]]
List emWithoutL1 (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
               const arma::rowvec Y, // observed number of reads for each read class j
               const double lambda,  // the tuning parameter for bias estimation, take as
               const arma::rowvec nzindex, // nzindex of non-zero read classes
               const double conv // , const int nThr = 1
){

  // initialization
  int J = X.n_cols; //number of equivalent read class
  int M = X.n_rows; //number of isoforms
  int p = nzindex.size() ; // number of non-zero b estimates

  arma::rowvec est(M+J - p); // bias parameter
  est.fill(1);

  arma::rowvec estPre = est - 1; // define theta


  int iter = 0; // iterator
  List est_out(5); // create a empty list of size 5

  // algorithm
  while(!all(abs(est - estPre) < conv)){

    //update iterator
    iter++;

    //store log likelihood
    estPre = est;

    //Rcout << "iter : " << iter << "\n";
    est_out = optim_thetaandb_rcpp(est,X, Y, lambda,nzindex);
    est = Rcpp::as<arma::rowvec>(est_out[0]);
    //Rcout << "est : " << est << "\n";


  }

  arma::rowvec b(X.n_cols);
  b.zeros();
  for(int i ; i < p ; i++){
    b(nzindex(i)) = est(X.n_rows+nzindex(i));
  }
   // returns
  List ret ;
  ret["theta"] = est.head(M) ;
  ret["b"] = b ;
  return(ret) ;


}





//' run_em
//'
//' @param XList sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param YList observed number of reads for each read class j
//' @param lambdaList, the tuning parameter for bias estimation
//' @param conv convergence threshold of likelihoods
//' @export
// [[Rcpp::export]]
List run_em (const List XList, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
             const List YList, // observed number of reads for each read class j
             const List lambdaList,
             const bool d,
             const int maxiter,
             const double conv, // , const int nThr = 1
             const int nthr,
             bool display_progress=true
             ){

  // initialization
  int ngene = XList.size() ; // number of genes

  Progress p(ngene, display_progress);

  List estOutput(ngene);


  omp_set_dynamic(0);
  omp_set_num_threads( nthr);


#pragma omp parallel
{
#pragma omp for
  for(int i = 0; i < ngene ; i++){
    if ( ! Progress::check_abort() ) {
      p.increment(); // update progress
      //Rcout << "i : " << i << "\n";
      //Rcout << "XList[i] : " << XList[i] << "\n";

      estOutput[i] = emWithL1(Rcpp::as<arma::mat>(XList[i]),
                              Rcpp::as<arma::rowvec>(YList[i]),
                              lambdaList[i], d, maxiter, conv);
    }
  }

}


  return(estOutput) ;


}



