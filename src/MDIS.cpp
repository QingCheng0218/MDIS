#include "RcppArmadillo.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void Leastsquare(arma::mat X, arma::vec y, vec& bh, vec& residual){
  //contail the intercept term.
  int p = X.n_cols;
  int n = X.n_rows;
  mat Xtmp = join_rows(ones(n, 1), X);
  
  // bh = solve(X.t() * X, X.t()*y);
  bh = solve(Xtmp, y);
  vec fittedvalue = zeros(p, 1);
  fittedvalue = Xtmp*bh;
  residual = y - fittedvalue;
}

// [[Rcpp::export]]
List statfun1(mat X, vec y){
  
  // ---------------------------------------------------------------------------------------------
  // _Yres1: E(Y-Xbeta)(X1-rhoX2)(X2-rhoX1);residual-based screening procedure 
  // _Yres2: EY(X1-rhoX2)(X2-rhoX1); (Y center)The centralized response-based screening procedure
  // ---------------------------------------------------------------------------------------------
  
  mat XC = join_rows(X, X.col(0) % X.col(1));
  int n = X.n_rows; int p = X.n_cols;
  
  double RhoEst; 
  
  RhoEst = mean(X.col(0) % X.col(1));
  
  mat RESX = zeros(n, p);
  RESX.col(0) = X.col(0) - RhoEst*X.col(1);
  RESX.col(1) = X.col(1) - RhoEst*X.col(0);
  
  
  vec RES1 = zeros(n, 1); vec bh1 = zeros(p + 2, 1);
  //vec bh = zeros(p, 1), residual = zeros(n, 1);
  Leastsquare(XC, y, bh1, RES1);
  // cout <<" bh1:" << bh1.t() << endl;
  // cout << "RES1:" << RES1.t()<< endl;
  
  double sigma1, sigma2;
  sigma1 = stddev(RESX.col(0) % RESX.col(1))*stddev(RES1);
  sigma2 = stddev(bh1[1] * X.col(0) % RESX.col(0) % RESX.col(1) + bh1[2] * X.col(1) % RESX.col(0) % RESX.col(1));
  
  // cout << "sigma1:" << sigma1<< endl;
  // cout << "sigma2:" << sigma2<< endl;
  double ASYMSTD_Yres1, ASYMSTD_Yres2, STATIC_Yres1, STATIC_Yres2, zscore_Yres1, zscore_Yres2;
  
  // residual-based method
  vec RES_Yres1 = zeros(n, 1), RES_Yres2 = zeros(n, 1), bh = zeros(p + 1, 1);
  Leastsquare(X, y, bh, RES_Yres1);
  ASYMSTD_Yres1 = sigma1 / sqrt(n);
  STATIC_Yres1 = mean(RESX.col(0) % RESX.col(1) % RES_Yres1);
  zscore_Yres1 = STATIC_Yres1 / ASYMSTD_Yres1;
  
  // centralized response-based method.
  RES_Yres2 = y;
  ASYMSTD_Yres2 = sqrt(pow(sigma1, 2) + pow(sigma2, 2)) / sqrt(n);
  STATIC_Yres2 = mean(RESX.col(0) % RESX.col(1) % RES_Yres2);
  zscore_Yres2 = STATIC_Yres2 / ASYMSTD_Yres2;
  
  // cout << "STATIC_Yres2:" << STATIC_Yres2 << endl;
  // cout << "ASYMSTD_Yres2:" << ASYMSTD_Yres2 << endl;
  List output = List::create(
    Rcpp::Named("zscore_Yres1") = zscore_Yres1,
    Rcpp::Named("zscore_Yres2") = zscore_Yres2
  );
  return output;
}


double statfun(mat X, mat XC, vec y, int Yres){
  
  int n = X.n_rows; int p = X.n_cols;
  
  double RhoEst; double STATIC; double ASYMSTD;
  
  RhoEst = mean(X.col(0) % X.col(1));
  mat RESX = zeros(n, p);
  RESX.col(0) = X.col(0) - RhoEst*X.col(1);
  RESX.col(1) = X.col(1) - RhoEst*X.col(0);
  
  vec RES = zeros(n, 1); vec bh = zeros(p, 1);
  vec RES1 = zeros(n, 1); vec bh1 = zeros(p + 1, 1);
  //vec bh = zeros(p, 1), residual = zeros(n, 1);
  Leastsquare(XC, y, bh1, RES1);
  
  double sigma1, sigma2;
  sigma1 = stddev(RESX.col(0) % RESX.col(1))*stddev(RES1);
  sigma2 = stddev(bh1[0] * X.col(0) % RESX.col(0) % RESX.col(1) + bh1[1] * X.col(1) % RESX.col(0) % RESX.col(1));
  
  if (Yres == 1) {
    Leastsquare(X, y, bh, RES);
    ASYMSTD = sigma1 / sqrt(n);
  }
  else {
    RES = y;
    ASYMSTD = sqrt(pow(sigma1, 2) + pow(sigma2, 2)) / sqrt(n);
  }
  
  
  STATIC = mean(RESX.col(0) % RESX.col(1) % RES);
  
  double zscore = STATIC / ASYMSTD;
  
  
  /*return Rcpp::List::create(
      Rcpp::Named(" sigma1 ") = sigma1,
  Rcpp::Named(" sigma2 ") = sigma2,
  Rcpp::Named(" STATIC ") = STATIC,
  Rcpp::Named(" ASYMSTD ") = ASYMSTD,
  Rcpp::Named(" zscore ") = zscore
  );*/
  return zscore;
}


mat statloop(int n, int p, double Rho, vec Kappa, int IterMax, double b1, double b2, int Yres){
  
  arma::dmat S0 = abs(repmat(linspace(1, p, p), 1, p) - repmat(linspace(1, p, p), 1, p).t());
  arma::dmat S = zeros(p, p);
  for (int i = 0; i < p; i = i + 1){
    for (int j = 0; j < p; j = j + 1){
      S(i, j) = pow(Rho, S0(i, j));
    }
  }
  vec Aone = ones(n, 1);
  
  mat ZScore = zeros(IterMax, Kappa.size());
  //mat POWER = zeros(Kappa.size(), pVALUE.size());
  vec B = zeros(3, 1);
  
  //double RhoEst; double STATIC; double ASYMSTD;
  int kMax = Kappa.size();
  for (int k = 0; k < kMax; k = k + 1){
    B[0] = b1; B[1] = b2; B[2] = Kappa(k);
    vec temp = zeros(IterMax, 1);
    for (int iter = 0; iter < IterMax; iter = iter + 1){
      
      mat X0 = zeros(n, p); mat X = zeros(n, p); mat XC = zeros(n, p + 1);
      vec Err = zeros(n, 1);
      vec y = zeros(n, 1);
      
      
      X0 = mvnrnd(zeros(p, 1), S, n).t();
      X = (X0 - repmat(mean(X0, 0), n, 1)) / repmat(stddev(X0, 0), n, 1);
      XC = join_rows(X, X.col(0) % X.col(1));
      Err = randn(n);
      //y = X.col(0)*b1 + X.col(1)*b2 + X.col(0) % X.col(1)*Kappa(k) + Err;
      y = XC*B + Err;
      y = y - mean(y);
      temp[iter] = statfun(X, XC, y, Yres);
    }
    ZScore.col(k) = temp;
  }
  
  
  return ZScore;
}

// [[Rcpp::export]]
imat comb(int p)
{
  //2- combinations from a given set S of n elements
  vector<int> myVector;
  std::string bitmask(2, 1); // K leading 1's
  bitmask.resize(p, 0); // N-K trailing 0's
  // cout << bitmask << endl;
  // print integers and permute bitmask
  do {
    for (int i = 0; i < p; ++i) // [0..N-1] integers
    {
      if (bitmask[i]){
        myVector.push_back(i);
      }
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  
  imat at = zeros<imat>(p*(p - 1) / 2, 2);
  for (int j = 0; j < p*(p - 1) / 2; j++){
    at(j, 0) = as_scalar(myVector[2 * j]);
    at(j, 1) = as_scalar(myVector[2 * j + 1]);
  }
  
  return at;
}

// [[Rcpp::export]]
List MDIS(mat X, vec y){
  
  int p = X.n_cols;
  int n = X.n_rows;
  imat at = comb(p);
  
  y = y - mean(y);
  uvec Trank_Yres1 = zeros<uvec>(at.n_rows, 1), Trank_Yres2 = zeros<uvec>(at.n_rows, 1);
  vec Zscore_Yres1 = zeros(at.n_rows, 1), Zscore_Yres2 = zeros(at.n_rows, 1);
  for (int j = 0; j < (int)at.n_rows; j++){
    mat X_tmp = zeros(n, 2), XC = zeros(n, 3);
    X_tmp = join_rows(X.col(at(j, 0)), X.col(at(j, 1)));
    // XC = join_rows(X, X.col(at(j, 0)) % X.col(at(j, 1)));
    List tmp = statfun1(X_tmp, y);
    double tmp1, tmp2;
    tmp1 = tmp["zscore_Yres1"]; tmp2 = tmp["zscore_Yres2"];
    Zscore_Yres1[j] = abs(tmp1);
    Zscore_Yres2[j] = abs(tmp2);
  }
  //return the index descend direction. the first one is the lagest index location.
  Trank_Yres1 = sort_index(Zscore_Yres1, "descend");
  Trank_Yres2 = sort_index(Zscore_Yres2, "descend");
  
  // return the location from 1 but not 0.
  Trank_Yres1 = Trank_Yres1 + 1;
  Trank_Yres2 = Trank_Yres2 + 1;
  
  List output = List::create(
    Rcpp::Named("Zscore_Yres1") = Zscore_Yres1,
    Rcpp::Named("Zscore_Yres2") = Zscore_Yres2,
    Rcpp::Named("Trank_Yres1") = Trank_Yres1,
    Rcpp::Named("Trank_Yres2") = Trank_Yres2
  );
  return output;
}
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
double jcisfun(mat X, vec y){
  int n = X.n_rows;
  double term1, term2, term3, term4;
  double res;
  vec x1 = X.col(0);
  vec x2 = X.col(1);
  term1 = abs(sum((x1 - mean(x1))%(x2 - mean(x2))%(y - mean(y))));
  term2 = sum(pow((x1 - mean(x1)), 2));
  term3 = sum(pow((x2 - mean(x2)), 2));
  term4 = sum(pow((y - mean(y)), 2));
  res = sqrt(n)*term1/sqrt(term2*term3*term4);
  return res;
}
// [[Rcpp::export]]
List JCIS(imat at, mat X, vec y){
  int n = X.n_rows;
  int p = X.n_cols;
  uvec rank_i = zeros<uvec>(p*(p - 1) / 2, 1);
  vec Zscore = zeros(at.n_rows, 1);
  for (int j = 0; j < p*(p - 1) / 2; j++){
    mat X_tmp = zeros(n, 2);
    X_tmp = join_rows(X.col(at(j, 0)), X.col(at(j, 1)));
    Zscore[j] = jcisfun(X_tmp, y);
  }
  //return the index descend direction. the first one is the lagest index location.
  rank_i = sort_index(Zscore, "descend");
  // return the location from 1 but not 0.
  rank_i = rank_i + 1;
  
  List output = List::create(
    Rcpp::Named("Zscore") = Zscore,
    Rcpp::Named("Rank") = rank_i
  );
  return output;
}
// [[Rcpp::export]]
List DIS(imat at, mat X, vec y){
  // direct interaction screening(DIS)
  int n = X.n_rows;
  int p = X.n_cols;
  uvec rank_i = zeros<uvec>(p*(p - 1) / 2, 1);
  vec Zscore = zeros(at.n_rows, 1);
  for (int j = 0; j < p*(p - 1) / 2; j++){
    mat X_tmp = zeros(n, 2);
    vec x1 = X.col(at(j, 0));
    vec x2 = X.col(at(j, 1));
    Zscore[j] = abs(sum(x1%x2%y));
  }
  //return the index descend direction. the first one is the lagest index location.
  rank_i = sort_index(Zscore, "descend");
  
  // return the location from 1 but not 0.
  rank_i = rank_i + 1;
  
  List output = List::create(
    Rcpp::Named("Zscore") = Zscore,
    Rcpp::Named("Rank") = rank_i
  );
  return output;
}
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
uvec proposal(imat at, mat X, vec y, int Yres){
  int n = X.n_rows;
  int p = X.n_cols;
  uvec rank_i = zeros<uvec>(p*(p - 1) / 2, 1);
  vec Zscore = zeros(at.n_rows, 1);
  for (int j = 0; j < p*(p - 1) / 2; j++){
    mat X_tmp = zeros(n, 2), XC = zeros(n, 3);
    X_tmp = join_rows(X.col(at(j, 0)), X.col(at(j, 1)));
    XC = join_rows(X, X.col(at(j, 0)) % X.col(at(j, 1)));
    Zscore[j] = abs(statfun(X_tmp, XC, y, Yres));
  }
  //return the index descend direction. the first one is the lagest index location.
  rank_i = sort_index(Zscore, "descend");
  return rank_i;
}


struct ResultStructure
{
  vec y;
  mat X;
};

ResultStructure simfun(mat S, mat B, vec beta, int n){
  ResultStructure res;
  int p = S.n_rows;
  mat X = zeros(n, p);
  vec y = zeros(n, 1);
  
  X = mvnrnd(zeros(p, 1), S, n).t();
  vec Err = randn(n, 1);
  vec INX = zeros(n, 1);
  for (int i = 0; i < n; i++){
    INX[i] = as_scalar(X.row(i)*B*X.row(i).t());
  }
  y = X*beta + INX + Err;
  res.X = X;
  res.y = y;
  return res;
}


// [[Rcpp::export]]
List proploop(mat S, mat B, vec beta, int n, uword IterMax){
  int p = S.n_rows;
  uword introw = p*(p - 1) / 2;
  umat rank_pr1 = zeros<umat>(introw, IterMax);
  umat rank_pr2 = zeros<umat>(introw, IterMax);
  imat at = comb(p);
  for (int iter = 0; iter < (int)(IterMax); iter = iter + 1){
    ResultStructure data;
    data = simfun(S, B, beta, n);
    mat X;
    vec y;
    X = data.X;
    y = data.y;
    rank_pr1.col(iter) = proposal(at, X, y, 1);
    rank_pr2.col(iter) = proposal(at, X, y, 2);
  }
  
  List output = List::create(
    Rcpp::Named("rank_pr1") = rank_pr1,
    Rcpp::Named("rank_pr2") = rank_pr2
  );
  return output;
}
