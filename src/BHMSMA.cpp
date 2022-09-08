//######## BHMSMAfMRI package C++ functions #########//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
 

//[[Rcpp::export(rng = false)]]
List glmcoef_sub(uword grid, arma::cube D_sub, arma::mat X)
{
	int p = X.n_cols;
	arma::cube GLM_coef_nonst(grid,grid,p), GLM_coef_st(grid,grid,p), GLM_coef_se(grid,grid,p);
	for(uword i=0; i<grid; i++)
	for(uword j=0; j<grid; j++)
	{
		arma::vec y = D_sub.subcube(i,j,0,  i,j,D_sub.n_slices-1);
		if(sum(y)!=0)
		{
			arma::vec coef = solve(X,y);
			GLM_coef_nonst.subcube(i,j,0, i,j,p-1) = coef;

			arma::vec res = y - X*coef;
			double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(y.n_elem-p);
			arma::vec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

			GLM_coef_se.subcube(i,j,0, i,j,p-1) = std_err;  // sqrt(coef_var);

			GLM_coef_st.subcube(i,j,0, i,j,p-1) = coef/std_err; 
		}
	}
	return Rcpp::List::create(
    _["GLM_coef_st"]=GLM_coef_st, 
    _["GLM_coef_se"]=GLM_coef_se 
  );
}


//[[Rcpp::export(rng = false)]]
double minus_ll(double C0, double C1, double C2, double C3, double C4, double C5, arma::uvec subs, uword grid, arma::mat waveletcoefmat )
{
	double aux = 0.0;
	double aux1;

	for(uword i=0; i<subs.n_elem; i++)
	for(uword l=0; l<log2(grid); l++)
	{
		arma::vec auxind = 2*regspace(0,l-1);
		uword sum1 = round(sum( exp(log(2)*auxind)*3 ));

		for(uword j=0; j<pow(2,2*l)*3; j++)
		{
			double d;
			if(l==0) d = waveletcoefmat(subs(i)-1,j); else d = waveletcoefmat(subs(i)-1, sum1 + j);
			double a_l = C0 * pow(2,-C1*l);
			double b_l = C2 * pow(2,-C3*l);
			double expect_w = std::min( 1.0, a_l/(a_l+b_l) );
			double c_l = C4 * pow(2,-C5*l);
			if (abs(d)<=39.6) 
			{
			aux1 = log(  exp( log(expect_w * pow(1+c_l,-.5)) - 0.5*pow(d,2)/(1+c_l) ) + exp( log(1-expect_w) - 0.5*pow(d,2) )  );  
			} else 
			{
			aux1 = log(  exp( log(expect_w * pow(1+c_l,-.5)) - 0.5*pow(39.6,2)/(1+c_l) ) + exp( log(1-expect_w) - 0.5*pow(39.6,2) )  );
			}
			aux = aux + aux1;
		}
	}
	return(-aux);
}


//[[Rcpp::export(rng = false)]]
double ll(double C0, double C1, double C2, double C3, double C4, double C5, arma::uvec subs, uword grid, arma::mat waveletcoefmat)
{
	return(- minus_ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) );
}


//[[Rcpp::export(rng = false)]]
arma::mat var_mle(double C0, double C1, double C2, double C3, double C4, double C5, arma::uvec subs, uword grid, arma::mat waveletcoefmat, arma::vec h)
{
	arma::mat IM = mat(6,6);

	IM(0,0) = -( ll(C0+h(0),C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0-h(0),C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) ) / pow(h(0),2);
	IM(0,1) = -( ll(C0+h(0),C1+h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0+h(0),C1-h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0-h(0),C1+h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0-h(0),C1-h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) ) / (4*h(0)*h(1));
	IM(0,2) = -( ll(C0+h(0),C1,C2+h(2),C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0+h(0),C1,C2-h(2),C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0-h(0),C1,C2+h(2),C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0-h(0),C1,C2-h(2),C3,C4,C5,subs,grid,waveletcoefmat) ) / (4*h(0)*h(2));
	IM(0,3) = -( ll(C0+h(0),C1,C2,C3+h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0+h(0),C1,C2,C3-h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0-h(0),C1,C2,C3+h(3),C4,C5,subs,grid,waveletcoefmat) + ll(C0-h(0),C1,C2,C3-h(3),C4,C5,subs,grid,waveletcoefmat) ) / (4*h(0)*h(3));
	IM(0,4) = -( ll(C0+h(0),C1,C2,C3,C4+h(4),C5,subs,grid,waveletcoefmat) - ll(C0+h(0),C1,C2,C3,C4-h(4),C5,subs,grid,waveletcoefmat) - ll(C0-h(0),C1,C2,C3,C4+h(4),C5,subs,grid,waveletcoefmat) + ll(C0-h(0),C1,C2,C3,C4-h(4),C5,subs,grid,waveletcoefmat) ) / (4*h(0)*h(4));
	IM(0,5) = -( ll(C0+h(0),C1,C2,C3,C4,C5+h(5),subs,grid,waveletcoefmat) - ll(C0+h(0),C1,C2,C3,C4,C5-h(5),subs,grid,waveletcoefmat) - ll(C0-h(0),C1,C2,C3,C4,C5+h(5),subs,grid,waveletcoefmat) + ll(C0-h(0),C1,C2,C3,C4,C5-h(5),subs,grid,waveletcoefmat) ) / (4*h(0)*h(5));
	IM(1,1) = -( ll(C0,C1+h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1-h(1),C2,C3,C4,C5,subs,grid,waveletcoefmat) ) / pow(h(1),2);
	IM(1,2) = -( ll(C0,C1+h(1),C2+h(2),C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1+h(1),C2-h(2),C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1-h(1),C2+h(2),C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1-h(1),C2-h(2),C3,C4,C5,subs,grid,waveletcoefmat) ) / (4*h(1)*h(2));
	IM(1,3) = -( ll(C0,C1+h(1),C2,C3+h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1+h(1),C2,C3-h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1-h(1),C2,C3+h(3),C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1-h(1),C2,C3-h(3),C4,C5,subs,grid,waveletcoefmat) ) / (4*h(1)*h(3));
	IM(1,4) = -( ll(C0,C1+h(1),C2,C3,C4+h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1+h(1),C2,C3,C4-h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1-h(1),C2,C3,C4+h(4),C5,subs,grid,waveletcoefmat) + ll(C0,C1-h(1),C2,C3,C4-h(4),C5,subs,grid,waveletcoefmat) ) / (4*h(1)*h(4));
	IM(1,5) = -( ll(C0,C1+h(1),C2,C3,C4,C5+h(5),subs,grid,waveletcoefmat) - ll(C0,C1+h(1),C2,C3,C4,C5-h(5),subs,grid,waveletcoefmat) - ll(C0,C1-h(1),C2,C3,C4,C5+h(5),subs,grid,waveletcoefmat) + ll(C0,C1-h(1),C2,C3,C4,C5-h(5),subs,grid,waveletcoefmat) ) / (4*h(1)*h(5));
	IM(2,2) = -( ll(C0,C1,C2+h(2),C3,C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2-h(2),C3,C4,C5,subs,grid,waveletcoefmat) ) / pow(h(2),2);
	IM(2,3) = -( ll(C0,C1,C2+h(2),C3+h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2+h(2),C3-h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2-h(2),C3+h(3),C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2-h(2),C3-h(3),C4,C5,subs,grid,waveletcoefmat) ) / (4*h(2)*h(3));
	IM(2,4) = -( ll(C0,C1,C2+h(2),C3,C4+h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2+h(2),C3,C4-h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2-h(2),C3,C4+h(4),C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2-h(2),C3,C4-h(4),C5,subs,grid,waveletcoefmat) ) / (4*h(2)*h(4));
	IM(2,5) = -( ll(C0,C1,C2+h(2),C3,C4,C5+h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2+h(2),C3,C4,C5-h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2-h(2),C3,C4,C5+h(5),subs,grid,waveletcoefmat) + ll(C0,C1,C2-h(2),C3,C4,C5-h(5),subs,grid,waveletcoefmat) ) / (4*h(2)*h(5));
	IM(3,3) = -( ll(C0,C1,C2,C3+h(3),C4,C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3-h(3),C4,C5,subs,grid,waveletcoefmat) ) / pow(h(3),2);
	IM(3,4) = -( ll(C0,C1,C2,C3+h(3),C4+h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3+h(3),C4-h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3-h(3),C4+h(4),C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3-h(3),C4-h(4),C5,subs,grid,waveletcoefmat) ) / (4*h(3)*h(4));
	IM(3,5) = -( ll(C0,C1,C2,C3+h(3),C4,C5+h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3+h(3),C4,C5-h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3-h(3),C4,C5+h(5),subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3-h(3),C4,C5-h(5),subs,grid,waveletcoefmat) ) / (4*h(3)*h(5));
	IM(4,4) = -( ll(C0,C1,C2,C3,C4+h(4),C5,subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3,C4-h(4),C5,subs,grid,waveletcoefmat) ) / pow(h(4),2);
	IM(4,5) = -( ll(C0,C1,C2,C3,C4+h(4),C5+h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4+h(4),C5-h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4-h(4),C5+h(5),subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3,C4-h(4),C5-h(5),subs,grid,waveletcoefmat) ) / (4*h(4)*h(5));
	IM(5,5) = -( ll(C0,C1,C2,C3,C4,C5+h(5),subs,grid,waveletcoefmat) - ll(C0,C1,C2,C3,C4,C5,subs,grid,waveletcoefmat) + ll(C0,C1,C2,C3,C4,C5-h(5),subs,grid,waveletcoefmat) ) / pow(h(5),2);

	IM(1,0) = IM(0,1);
	IM(2,0) = IM(0,2);
	IM(2,1) = IM(1,2);
	IM(3,0) = IM(0,3);
	IM(3,1) = IM(1,3);
	IM(3,2) = IM(2,3);
	IM(4,0) = IM(0,4);
	IM(4,1) = IM(1,4);
	IM(4,2) = IM(2,4);
	IM(4,3) = IM(3,4);
	IM(5,0) = IM(0,5);
	IM(5,1) = IM(1,5);
	IM(5,2) = IM(2,5);
	IM(5,3) = IM(3,5);
	IM(5,4) = IM(4,5);

	arma::mat out = inv(IM);
	return(out);
}

 
/// Writing the loglikelihood ln(p(w_lj/y)) as a function of w_lj
//[[Rcpp::export(rng = false)]]
arma::vec LL(arma::vec w, uword l, uword j, uword n, arma::mat waveletcoefmat, double C0, double C1, double C2, double C3, double C4, double C5)
{
	double d;
	arma::vec aux(w.n_elem), aux1;
	arma::vec auxind = 2*regspace(0,l-1);
	uword sum1 = round(sum( exp(log(2)*auxind)*3 ));

	double c_l = C4 * pow(2,-C5*l);
		
	for(uword i=0; i<n; i++)
	{
		if(l==0) d = waveletcoefmat(i,j-1); else d = waveletcoefmat(i, sum1+j-1);
		
		if (abs(d) <= 39.6)
	    {
	      aux1 = log(  exp( log(w * pow(1+c_l,-.5)) - 0.5*pow(d,2)/(1+c_l) ) + exp( log(1-w) - 0.5*pow(d,2) )  );
	    } else
	    {
	      aux1 = log(  exp( log(w * pow(1+c_l,-.5)) - 0.5*pow(39.6,2)/(1+c_l) ) + exp( log(1-w) - 0.5*pow(39.6,2) )  );
	    }
    	aux =  aux + aux1 + (C0 * pow(2,-C1*l) - 1)*log(w) + (C2 * pow(2,-C3*l) - 1) * log(1-w);
	}
	return(aux);
}


/// Writing piklj as a function of w_lj 
//[[Rcpp::export(rng = false)]]
arma::vec pklj(arma::vec w, uword l, uword j, uword i, arma::mat waveletcoefmat, double C4, double C5)
{
	double d,aux;
	arma::vec O,p;
	arma::vec auxind = 2*regspace(0,l-1);
	uword sum1 = round(sum( exp(log(2)*auxind)*3 ));
	if(l==0) d = waveletcoefmat(i-1,j-1); else d = waveletcoefmat(i-1, sum1+j-1);
	double c_l = C4 * pow(2,-C5*l);

	if(0.5 * c_l/(1+c_l) * pow(d,2) <= 700) aux = exp(0.5 * c_l/(1+c_l) * pow(d,2)); else aux = exp(700);  // Note: Adjustment
	O = pow(1+c_l,-0.5) * w/(1-w) * aux;
	p = O/(1+O);
	return(p);
}

    
		
/// Using Trapezoidal rule to evaluate p_klj bar 
//[[Rcpp::export(rng = false)]]
arma::mat pklj_bar(uword grid, uword n, arma::mat waveletcoefmat, double C0, double C1, double C2, double C3, double C4, double C5)
{
	uword N = 1000;
	uword a = 0;
	uword b = 1;
	arma::mat pklj_bar(n, pow(grid,2)-1);

	arma::vec w_grid = a + regspace(0,N) * (b-a)/N;

	for(uword l=0; l<log2(grid); l++)
	for(uword j=0; j<pow(2,2*l)*3; j++)
	{
		arma::vec w_density_unnorm = exp( LL(w_grid,l,j+1,n,waveletcoefmat,C0,C1,C2,C3,C4,C5) );
		for(uword id=0; id<w_grid.n_elem; id++) 
		{
		  if(!is_finite(w_density_unnorm(id))) w_density_unnorm(id) = 1;
		}	
		double sum_w = sum(w_density_unnorm.elem(regspace<uvec>(1,N-1))) + 0.5 * (w_density_unnorm(0)+w_density_unnorm(N));
		arma::vec  w_density_norm = w_density_unnorm / sum_w;

		for(uword i=0; i<n; i++)
		{
			arma::vec pklj_func_w = pklj(w_grid,l,j+1,i+1,waveletcoefmat,C4,C5);

			pklj_func_w(N) = 1;

			arma::vec prod_p_wd = w_density_norm % pklj_func_w;
			double integral = sum(prod_p_wd.elem(regspace<uvec>(1,N-1)))  + 0.5*(prod_p_wd(0)+prod_p_wd(N));

			arma::vec auxind = 2*regspace(0,l-1);
			uword sum1 = round(sum( exp(log(2)*auxind)*3 ));
			if(l==0)  pklj_bar(i,j) = integral; else pklj_bar(i,sum1+j) = integral;
		}
	}	
	return(pklj_bar);
}


/// Posterior mean/median of wavelet coefficients
//[[Rcpp::export(rng = false)]]
List post_wavelet_coef(uword grid, uword n, arma::mat waveletcoefmat, arma::mat pkljbar, double C4, double C5)
{
	arma::mat PostMeanWaveletCoeff(n,pow(grid,2)-1);
	arma::mat PostMedianWaveletCoeff(n,pow(grid,2)-1);
	for(uword i=0; i<n; i++)
	{
		arma::vec PostMean(pow(grid,2)-1), PostMedian(pow(grid,2)-1);
		double d, p;
		int counter = 0;
		for(uword l=0; l<log2(grid); l++)
	    {
	        arma::vec auxind = 2*regspace(0,l-1);
	      	uword sum1 = round(sum( exp(log(2)*auxind)*3 ));
    		for(uword j=0; j<pow(2,2*l)*3; j++)
    		{
    			if(l==0) 
    			{
    				d = waveletcoefmat(i,j);
    				p = pkljbar(i,j);
    			} else
    			{
    				d = waveletcoefmat(i, sum1+j);
    				p = pkljbar(i, sum1+j);
    			}
    			
    			double c_l = C4 * pow(2,-C5 * l);
    			double mean = p * c_l/(1+c_l) * d;
    			PostMean(counter) = mean;
    			
    			NumericVector p1 = wrap(p);
    			NumericVector q_norm = qnorm(0.5/p1);
    			double q_norm1 = q_norm(0);
          
          		double median;
    			if(p <= 0.5) median = 0; else median = c_l/(1+c_l) * d - sign(d) * pow(c_l/(1+c_l),.5) * q_norm1;
    			PostMedian(counter) = median;
    			counter = counter + 1;
    		}
		}
		PostMeanWaveletCoeff.row(i) = PostMean.t();
		PostMedianWaveletCoeff.row(i) = PostMedian.t();
	}
	return Rcpp::List::create(
    _["PostMeanWaveletCoeff"]=PostMeanWaveletCoeff, 
    _["PostMedianWaveletCoeff"]=PostMedianWaveletCoeff 
	);
}


/// Set seed, supplied from R, inside Rcpp
// [[Rcpp::export]]
void set_seed(unsigned int seed) 
{
	Rcpp::Environment base_env("package:base");
	Rcpp::Function set_seed_r = base_env["set.seed"];
	set_seed_r(seed);  
}


/// Posterior sample generation
// [[Rcpp::export]]
arma::cube post_samp(uword nsample, uword grid, uword n, arma::mat waveletcoefmat, arma::mat pkljbar, double C4, double C5, uword seed)
{
	set_seed(seed);
	arma::cube dsamp_cube(n,pow(grid,2)-1,nsample);
	for(uword i=0; i<n; i++)
	{
		arma::mat dsamp_mat(pow(grid,2)-1,nsample);
		for(uword g=0; g<nsample; g++)
		{
			arma::vec dsamp(pow(grid,2)-1);
			for(uword l=0; l<log2(grid); l++)
		    {
		        arma::vec auxind = 2*regspace(0,l-1);
		      	uword sum1 = round(sum( exp(log(2)*auxind)*3 ));
				for(uword j=0; j<pow(2,2*l)*3; j++)
				{
					double d, p;
					if(l==0) 
					{
						d = waveletcoefmat(i,j);
						p = pkljbar(i,j);
					} else
					{
						d = waveletcoefmat(i, sum1+j);
						p = pkljbar(i, sum1+j);
					}
					double c_l = C4 * pow(2,-C5 * l);

					double prob = std::min(1.0,p);
					NumericVector U = rbinom(1,1,prob);
					double U1 = U(0);
					if(l==0) 
				    {
				        if(U1==1) dsamp(j) = R::rnorm( c_l/(1+c_l) * d, sqrt(c_l/(1+c_l) )); else dsamp(j) = 0; 
				    } else 
				    {
				        if(U1==1) dsamp(sum1+j) = R::rnorm( c_l/(1+c_l) * d, sqrt(c_l/(1+c_l) ));  else dsamp(sum1+j) = 0;
				    }			
				}
			}
			dsamp_mat.col(g) = dsamp;
		}	
		dsamp_cube.row(i) = dsamp_mat;	
	}
	return(dsamp_cube);
}



















