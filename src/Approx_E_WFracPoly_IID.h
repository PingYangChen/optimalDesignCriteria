// DECLARE FUNCTIONS

// BODY
double Approx_E_WFracPoly_IID(const arma::irowvec &poly, const double &u, const double &v, 
                              const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &m)
{
  arma::uword nSupp = DESIGN.n_rows;
  arma::uword dSupp = DESIGN.n_cols;
  
	m.set_size(OBJ.dmPara, OBJ.dmPara); m.zeros();

	double u = p(p.n_elem - 2); 
	double v = p(p.n_elem - 1); 
		
	for (arma::uword i = 0; i < nSupp; i++) {
		arma::rowvec supp = DESIGN.row(i);
		double lambda = std::pow(1.0 + supp(0), u)*std::pow(1.0 - supp(0), v);
		rowvec f(OBJ.dmPara, fill::ones);
		for (arma::uword j = 0; j < (OBJ.dmPara - 1); j++) {
			f(j+1) = std::pow(supp(0), -1.0*((double)(j+1)));
		}
		m += WT(i) * lambda * (f.t() * f);
	}
	double val = matrixMean(m, 2, 0); 
	return val;
}

