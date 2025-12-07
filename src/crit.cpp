
#include "common.h"
#include "Approx_E_WPoly_IID.h"

//[[Rcpp::export]]
Rcpp::List deisgnCritCpp(arma::rowvec x, arma::uword nSupp, arma::uword dSupp, arma::rowvec poly, double u, double v, double a, double b) {
  
  // Get Design
  arma::mat DESIGN;	arma::rowvec WT;
  v2d(DESIGN, WT, x, nSupp, dSupp);
  arma::mat m;
  double val = Approx_E_WPoly_IID(poly, u, v, a, b, DESIGN, WT, m);
  return List::create(Named("crit") = wrap(val),
                      Named("fim") = wrap(m));
}


//[[Rcpp::export]]
Rcpp::List equivCpp(arma::rowvec grid, arma::rowvec x, arma::uword nSupp, arma::uword dSupp, arma::uword &ctype, arma::rowvec poly, double u, double v, double a, double b) {
  
  // Get Design
  arma::mat DESIGN;	arma::rowvec WT;
  v2d(DESIGN, WT, x, nSupp, dSupp);
  
  Rcpp::List tmpL = deisgnCritCpp(x, nSupp, dSupp, poly, u, v, a, b);
  double critVal = tmpL["crit"]; 
  arma::mat FIM = tmpL["fim"]; 
  
  arma::uword dpoly = poly.n_elem;
  
  arma::mat inv_FIM = FIM;
  if (ctype != 2) {
    inv_FIM = inv_sympd(FIM);
  }
  
  arma::uword nGrid = grid.n_elem;
  
  arma::mat ETA_Dev(dpoly, nGrid, fill::zeros);
  for (arma::uword i = 0; i < dpoly; i++) {
    ETA_Dev.row(i) = Approx_E_WPoly_IID_fDev(poly(i), u, v, a, b, grid.t()); 
  }
  
  arma::mat DD(1, nGrid, fill::zeros);
  arma::rowvec inequ(2);
  for (arma::uword i = 0; i < nGrid; i++) {
    if (ctype == 2) {
      inequ = directionalDerivative(ETA_Dev.col(i).t(), FIM, critVal, ctype);
    } else {
      inequ = directionalDerivative(ETA_Dev.col(i).t(), inv_FIM, critVal, ctype);
    }
    DD(0, i) = inequ(0) - inequ(1);
  } 
  
  return List::create(Named("grid") = wrap(grid),
                      Named("DD") = wrap(DD));
}

