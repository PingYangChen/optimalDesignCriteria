// Rcpp Header File
#include <cmath>
#include <RcppArmadillo.h>
#include <R.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// DECLARE FUNCTIONS
void matrixPrintf(const mat &m);
void rvecPrintf(const rowvec &v);
void v2d(arma::mat &DESIGN, arma::rowvec &WT, const arma::rowvec &x, const arma::uword &nSupp, const arma::uword &dSupp);
double matrixMean(const arma::mat &m, const arma::uword &ctype = 0, const arma::rowvec &linear = arma::rowvec());
arma::rowvec directionalDerivative(const arma::rowvec &fDev, const arma::mat &inv_FIM, const double &criVal, const arma::uword &ctype = 0, const arma::rowvec &linear = arma::rowvec());

// BODY
void matrixPrintf(const arma::mat &m) 
{
  for (arma::uword i = 0; i < m.n_rows; i++) {
    for (arma::uword j = 0; j < m.n_cols; j++) Rprintf("%4.4f\t", m(i,j));
    Rprintf("\n");
  }
  Rprintf("\n\n");
}

void rvecPrintf(const arma::rowvec &v) 
{
  for (arma::uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i)); 	
  Rprintf("\n\n");
}

// Approximation Design Format Change
void v2d(arma::mat &DESIGN, arma::rowvec &WT, const arma::rowvec &x, const arma::uword &nSupp, const arma::uword &dSupp) {
  DESIGN.set_size(nSupp, dSupp); DESIGN.zeros();
  WT.set_size(nSupp);
  for (arma::uword i = 0; i < nSupp; i++) { 
    DESIGN.row(i) = x.subvec(i*dSupp, (i+1)*dSupp - 1); 
  }
  // Get Design Weight
  arma::rowvec wcumsin(nSupp, fill::zeros), wcos(nSupp, fill::zeros);
  arma::rowvec ang = x.subvec(nSupp*dSupp, nSupp*dSupp + nSupp - 2);
  arma::rowvec wsin = arma::sin(ang);
  wcumsin(0) = 1.0;
  for (arma::uword i = 1; i < nSupp; i++) { 
    wcumsin(i) = wcumsin(i-1)*wsin(i-1); 
  }
  wcos(nSupp - 1) = 1.0; 
  wcos.subvec(0, nSupp - 2) = arma::cos(ang);
  arma::rowvec w_root = wcumsin % wcos;
  WT = w_root % w_root;
}


double matrixMean(const arma::mat &m, const arma::uword &ctype, const arma::rowvec &linear)
{
	double val = 1e20;
	if (!m.has_nan()) {
		switch (ctype) {
			case 0: { // D
				//if (std::abs(arma::det(m)) > 0.0) {
					//val = (-1.0)*std::log(arma::det(m)); 
					double tmp_v, sign_v; 
					arma::log_det(tmp_v, sign_v, m); 
					if (sign_v > 0) { val = -1.0*tmp_v; }
				//}
				break;  
			}
			case 1: { // A
				if (rcond(m) > (datum::eps*((double)m.n_cols))) val = arma::trace(m.i()); 
				break;						 
			}
			case 2: { // E
				//if (rcond(m) > (datum::eps*((double)m.n_cols))) {
				//vec eigval; eig_sym(eigval, m.i()); val = eigval.max(); 	
				arma::vec eigval; 
			  arma::eig_sym(eigval, m);
			  val = -1.0*eigval.min(); 	
				//}
				break;
			}
			case 3: { // c
				if (rcond(m) > (datum::eps*((double)m.n_cols))) {
					val = as_scalar(linear * m.i() * linear.t()); 
				}
				break; 
			}
		}	
	}
	if (std::isnan(val)) val = 1e20;
	return val;
}


arma::rowvec directionalDerivative(const arma::rowvec &fDev, const arma::mat &inv_FIM, const double &criVal, const arma::uword &ctype, const arma::rowvec &linear) 
{
	arma::rowvec DD(2); 
	switch (ctype) {
		case 0: { // D
	    DD << arma::as_scalar(fDev*inv_FIM*fDev.t()) << (double)inv_FIM.n_cols << endr; 
	    break; 
	  }
	  case 1: { // A
		  DD << arma::as_scalar(fDev*inv_FIM*inv_FIM*fDev.t()) << criVal << endr; 
		  break; 
		}
		case 2: { // E
			arma::mat FIM = inv_FIM;
			arma::vec eigVal; 
			arma::mat eigVec; 
			arma::eig_sym(eigVal, eigVec, FIM); 
			arma::uword min_eig_val = index_min(eigVal);
			arma::vec min_eig_vec = eigVec.col(min_eig_val);
			double tmp = arma::as_scalar(fDev*(min_eig_vec*min_eig_vec.t())*fDev.t());
			DD << tmp << eigVal.min() << endr;
			break;
		}
		case 3: {
			double tmp = arma::as_scalar(fDev*inv_FIM*linear.t()); 
			DD << tmp*tmp << criVal << endr;
			break;
		}
	}
	return DD;
}
