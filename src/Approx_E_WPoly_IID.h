// DECLARE FUNCTIONS
rowvec Approx_E_WPoly_IID_fDev(const double &ipoly, const double &u, const double &v, const double &a, const double &b, const vec &x);

// BODY
double Approx_E_WPoly_IID(const arma::rowvec &poly, const double &u, const double &v, const double &a, const double &b,
                          const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &m)
{
  arma::uword nSupp = DESIGN.n_rows;
  arma::uword dSupp = DESIGN.n_cols;
	arma::uword dpoly = poly.n_elem;
	m.set_size(dpoly, dpoly); m.zeros();
	
	arma::mat ETA_Dev(dpoly, nSupp, fill::zeros);
	for (arma::uword i = 0; i < dpoly; i++) {
	  ETA_Dev.row(i) = Approx_E_WPoly_IID_fDev(poly(i), u, v, a, b, DESIGN.col(0)); 
	}
	for (arma::uword iSupp = 0; iSupp < nSupp; iSupp++) {
		m += WT(iSupp) * (ETA_Dev.col(iSupp) * ETA_Dev.col(iSupp).t());
	}
	double val = matrixMean(m, 2);
	return val;
}

// SUBFUNCTIONS
arma::rowvec Approx_E_WPoly_IID_fDev(const double &ipoly, const double &u, const double &v, const double &a, const double &b, const vec &x)
{
	arma::vec f(x.n_elem, fill::ones);
	arma::vec lambda = arma::sqrt(arma::pow(x - a, u) % arma::pow(b - x, v));
	if (!std::isfinite(ipoly)) {
		f = arma::log(x) % lambda;
	} else {
		f = arma::pow(x, ipoly) % lambda;
	}
	return f.t();
}



