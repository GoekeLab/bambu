# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;

//' Create likelihood objective function to estimate bias parameters
//'
//' @param par A row vector of bias parameters for each read class
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param theta estimates for isoform expression
double obj_fun_rcpp(const arma::rowvec par,
                    const arma::mat X,
                    const arma::rowvec Y,
                    const double lambda,
                    const arma::rowvec theta
){
  int M = X.n_rows ; // number of transcripts
  int N = X.n_cols ; // number of read classes

  //Rcout << "obj_fun : 1 \n";
  arma::mat a_adj_mat(M,N); // define tmp variables

  //Rcout << "obj_fun : 2 \n";
  for(int i = 0; i < N; i++){
    a_adj_mat.col(i) = X.col(i) * exp(par(i));
  }

  //Rcout << "obj_fun : 3 \n";
  arma::rowvec a_adj_mat_rowSum = arma::sum(a_adj_mat.t(),0);
  for(int i = 0; i < M; i++){
    a_adj_mat.row(i) = a_adj_mat.row(i) / a_adj_mat_rowSum(i) * theta(i);
  }

  //Rcout << "obj_fun : 4 \n";
  arma::rowvec column_agg = arma::sum(a_adj_mat,0);

  //Rcout << "colum_agg : " << column_agg  << "\n";
  //Rcout << "obj_fun : 5 \n";
  arma::rowvec log_column_agg(N);
  for(int i = 0; i < N; i++){
    //Rcout << "colum_agg : " << i  << "\n";
    log_column_agg(i) = std::log(column_agg(i));
    //Rcout << "colum_agg : " << i  << "\n";
  }
  //  = log(column_agg);

  log_column_agg.replace(arma::datum::nan, 0); // replace each NaN with 0
  // Compute objective value
  double obj_val = -(sum(Y%log_column_agg-column_agg)-lambda*sum(abs(par)));

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
List optim_rcpp(const arma::rowvec par,
                const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
                const arma::rowvec Y, // observed number of reads for each read class j
                const double lambda,  // the tuning parameter for bias estimation, take as
                const arma::rowvec theta){

  // Extract Rs optim function
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  // Call the optim function from R in C++
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = par,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                 Rcpp::_["method"] = "Nelder-Mead",
                                 Rcpp::_["hessian"] = "FALSE",
                                 Rcpp::_["X"] = X,
                                 Rcpp::_["Y"] = Y,
                                 Rcpp::_["lambda"] = lambda,
                                 Rcpp::_["theta"] = theta);


  // Return estimated values
  return opt_results;
}


//' EM algorithm without L1-penalty
//'
//' @param par A row vector of parameters for isoform expression
//' @param X sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
//' @param Y observed number of reads for each read class j
//' @param lambda, the tuning parameter for bias estimation
//' @param conv convergence threshold of likelihoods
//' @export
// [[Rcpp::export]]
List emWithoutL1 (const arma::mat X, // sampling probability matrix, (i,j) = 1 if read class j is potentially from transcript i, otherwise 0
                  const arma::rowvec Y, // observed number of reads for each read class j
                  const arma::rowvec par, // bias parameter treated as fixed while estimating theta
                  const double lambda,  // the tuning parameter for bias estimation, take as
                  const double conv = 0.001 // , const int nThr = 1

) {
  // inputs
  int M = X.n_rows ; // number of transcripts
  int N = X.n_cols ; // number of read classes
  double logLDiff = 1; // initialize

  int iter = 0; // iterator
  // omp_set_num_threads(nThr) ; // using multiple threads

  // containers
  arma::rowvec theta(M) ; // define theta
  theta.fill(1) ; // initialize thetas to 1/M
  arma::mat probMat(M, N) ; // define prob matrix
  probMat.fill(1) ;



  // initialize Xb
  arma::mat Xb(M,N) ;

  for (int i = 0 ; i < N ; i++) {
    Xb.col(i) = X.col(i)*std::exp(par(i));
  }
  arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(),0); // to be used in theta estimation

  // normalized Xb
  for(int i = 0; i < M; i++){
    Xb.row(i) = Xb.row(i) / summed_by_row_Xb(i);
  }


  arma::mat adjP_Xb(M,N) ; // adjusted p-values = sampling probability X theta X b, with bias parameter included

  for (int i = 0 ; i < M ; i++) {
    adjP_Xb.row(i) = Xb.row(i)*theta(i);
  }



  arma::rowvec summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0); // will be used in logLikelihood estimation as well as probMat update
  // normalize adjusted p-values within each read
  // the complete version of probMat, it will be equivalent to
  for (int i = 0; i < N; i++) {
    probMat.col(i) = adjP_Xb.col(i) * Y(i) / summed_by_col_adjP_Xb(i);
  }

  // initialize logLikelihoods
  arma::mat probMat_int(M,N);
  for (int i = 0; i < N; i++) {
    probMat_int.col(i) = probMat.col(i) * std::log(summed_by_col_adjP_Xb(i));
  }

  arma::mat logLMat = probMat_int - adjP_Xb;
  logLMat.replace(arma::datum::nan, 0); // replace each NaN with 0
  float logLNext = accu(logLMat) ; // logLikelihood to be maximized over
  float logLPre;


  // algorithm
  while(!((logLDiff <= conv)&(logLDiff >= 0))){

    //update iterator
    iter++;

    //store log likelihood
    logLPre = logLNext;

    // update theta
    // probMat.replace(arma::datum::nan, 0);
    //  );
    for(int i = 0 ; i < M ; i++){
      theta(i) =  sum(probMat.row(i)) ;

    }


    // update adjust_Xb as theta is updated
    for (int i = 0 ; i < M ; i++) {
      adjP_Xb.row(i) = Xb.row(i)*theta(i);
    }

    // adjP.replace(arma::datum::nan, 0);
    summed_by_col_adjP_Xb = arma::sum(adjP_Xb,0);
    for (int i = 0; i < N; i++) {
      probMat.col(i) = adjP_Xb.col(i) * Y(i) / summed_by_col_adjP_Xb(i);
    }

    //update log likelihood
    for (int i = 0; i < N; i++) {
      probMat_int.col(i) = probMat.col(i) * std::log(summed_by_col_adjP_Xb(i));
    }
    logLMat =  probMat_int - adjP_Xb;
    logLMat.replace(arma::datum::nan, 0); // replace each NaN with 0
    logLNext = accu(logLMat) ;



    //update log likelihood difference
    logLDiff = logLNext - logLPre;

    // printing value of vector

  }


  // Rcout << "The increment in log likihood : " << logLNext << "\n";
  //  Rcout << "The increment in log likihood : " << logLPre << "\n";
  //  Rcout << "The increment in log likihood : " << logLDiff << "\n";


  // returns
  List ret ;
  ret["theta"] = theta ;
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
               const double conv = 0.001 // , const int nThr = 1
){

  // initialization
  double logLDiff = 1; // logLikelihood difference
  int M = X.n_rows ; // number of transcripts
  int N = X.n_cols ; // number of read classes

  arma::rowvec par(N); // bias parameter
  par.fill(0);
  arma::rowvec theta(M); // theta estimate
  theta.fill(1);

  double logLPre;
  //Rcout << "begin : 1 \n";
  double logLNext = obj_fun_rcpp(par, X, Y, lambda, theta);
  //Rcout << "logLNext : 2 \n";

  int iter = 0; // iterator
  List theta_out(2); // create a empty list of size 5
  List b_out(5);
  // algorithm
  while(!((logLDiff <= conv)&(logLDiff >= 0))){

    //update iterator
    iter++;

    //store log likelihood
    logLPre = logLNext;

    //Rcout << "theta_est : " << iter << "\n";
    theta_out = emWithoutL1(X,Y,par,lambda, conv);
    theta = Rcpp::as<arma::rowvec>(theta_out[0]);

    //Rcout << "b_est : " << iter << "\n";
    b_out = optim_rcpp(par,X, Y, lambda,theta);
    par = Rcpp::as<arma::rowvec>(b_out[0]);

    //Rcout << "log_lik : " << iter << "\n";

    logLNext = obj_fun_rcpp(par, X, Y, lambda, theta);



    //update log likelihood difference
    logLDiff = logLPre - logLNext;

    // printing value of vector

  }


  // Rcout << "The increment in log likihood : " << logLNext << "\n";
  //  Rcout << "The increment in log likihood : " << logLPre << "\n";
  //  Rcout << "The increment in log likihood : " << logLDiff << "\n";


  // returns
  List ret ;
  ret["theta"] = theta ;
  ret["b"] = par ;
  ret["iter"] = iter;
  return(ret) ;


}
