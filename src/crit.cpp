
#include "common.h"
#include "Approx_HeterPoly.h"

//[[Rcpp::export]]
Rcpp::List deisgnCritCpp(arma::rowvec x, arma::uword nSupp, arma::uword dSupp, arma::uword &ctype,
                         arma::rowvec poly, double u, double v, double a, double b) {
  // Get Design
  arma::mat DESIGN;	arma::rowvec WT;
  v2d(DESIGN, WT, x, nSupp, dSupp);
  // Set a empty information matrix 
  arma::mat m;
  // Input the design requirements and the empty matrix
  // to calculate the criterion value
  double val = Approx_HeterPoly(poly, u, v, a, b, DESIGN, WT, ctype, m);
  // Return results
  return List::create(Named("crit") = wrap(val), // Criterion value
                      Named("fim") = wrap(m));   // the corresponding information matrix
}


//[[Rcpp::export]]
Rcpp::List equivCpp(arma::rowvec &grid, arma::rowvec &x, arma::uword &nSupp, arma::uword &dSupp, 
                    arma::uword &ctype, 
                    arma::rowvec &poly, double &u, double &v, double &a, double &b) {
  
  // Get Design
  arma::mat DESIGN;	arma::rowvec WT;
  v2d(DESIGN, WT, x, nSupp, dSupp);
  // Compute Design Criterion
  Rcpp::List tmpL = deisgnCritCpp(x, nSupp, dSupp, ctype, poly, u, v, a, b);
  double critVal = tmpL["crit"]; 
  arma::mat FIM = tmpL["fim"]; 
  // Prepare for the Computation of Directional Derivatives
  arma::uword dpoly = poly.n_elem;
  arma::mat inv_FIM = FIM;
  if (ctype != 2) {
    inv_FIM = inv_sympd(FIM);
  }
  // Compute Directional Derivatives over the given grid
  arma::uword nGrid = grid.n_elem;
  arma::mat ETA_Dev(dpoly, nGrid, fill::zeros);
  for (arma::uword i = 0; i < dpoly; i++) {
    ETA_Dev.row(i) = Approx_HeterPoly_fDev(poly(i), u, v, a, b, grid.t()); 
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
  // Return results
  return List::create(Named("grid") = wrap(grid), // user-specified grids
                      Named("DD") = wrap(DD));    // Directional Derivatives over the given grid
}

