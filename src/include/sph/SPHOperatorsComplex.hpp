/*
 * SPHOperatorsComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: martin
 */

#ifndef SPHOPERATORS_COMPLEX_HPP_
#define SPHOPERATORS_COMPLEX_HPP_

#include <sph/SPHData.hpp>
#include <sph/SPHIdentities.hpp>

class SPHOperatorsComplex	:
		public SPHIdentities
{
	friend SPHConfig;

public:
	void setup(
			SPHConfig *i_sphConfig
	)
	{
	}

	~SPHOperatorsComplex()
	{
	}


#if 0
public:
	/**
	 * Compute the coriolis effect on the given velocity
	 *
	 * multiply (n=1, m=0) with \Omega sqrt(8/3)
	 */
	SPHDataComplex coriolis(
			const SPHDataComplex &i_sph_data,
			double i_coriolis_frequency
	)	const
	{
		// THIS IS NOT WORKING
		assert(false);
		SPHDataComplex out_sph_data(i_sph_data);

		out_sph_data.request_data_spectral();

		std::size_t idx = LM(i_sph_data.sphConfig->shtns, 1, 0);
		out_sph_data.data_spat[idx] *= i_coriolis_frequency*std::sqrt(8.0/3.0);

		return out_sph_data;
	}
#endif

public:
	/**
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
	SPHDataComplex diff_lon(
			const SPHDataComplex &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data(i_sph_data.sphConfig);

		// compute d/dlambda in spectral space
#pragma omp parallel for
		for (int n = 0; n <= i_sph_data.sphConfig->spec_n_max; n++)
		{
			int idx = i_sph_data.sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out_sph_data.data_spec[idx] = i_sph_data.data_spec[idx]*std::complex<double>(0, m);
				idx++;
			}
		}
		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d mu f(lambda,mu)
	 *
	 * sqrt(1-mu*mu)*d/dmu P_n^m = ...
	 */
	SPHDataComplex diff_lat_mu(
			const SPHDataComplex &i_sph_data
	)	const
	{

		SPHDataComplex out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, std::complex<double> &o_data)
				{
					o_data /= (1.0-mu*mu);
				}
		);

		return out_sph_data;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d phi f(lambda,mu)
	 */
	SPHDataComplex diff_lat_phi(
			const SPHDataComplex &i_sph_data
	)	const
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SPHDataComplex grad_lon(
			const SPHDataComplex &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data = diff_lon(i_sph_data);

		out_sph_data.request_data_spatial();

		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, std::complex<double> &o_data)
				{
					double cos_phi = sqrt(1.0-mu*mu);
					o_data /= cos_phi;
				}
		);

		return out_sph_data;
	}


#if 0
	static double fac_double(double v)
	{
		double retval = 1;
	    for (double i = 1; i <= v; ++i)
	    	retval *= i;
	    return retval;
	}
#endif


	SPHDataComplex spec_one_minus_mu_squared_diff_lat_mu(
			const SPHDataComplex &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();
		SPHConfig *sphConfig = i_sph_data.sphConfig;

		SPHDataComplex out_sph_data = SPHDataComplex(sphConfig);


#pragma omp parallel for
		for (int n = 0; n <= i_sph_data.sphConfig->spec_n_max; n++)
		{
			int idx = i_sph_data.sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out_sph_data.data_spec[idx] =
						((-n+1.0)*R(n-1,m))*i_sph_data.spec_get(n-1, m) +
						((n+2.0)*S(n+1,m))*i_sph_data.spec_get(n+1, m);

				idx++;
			}
		}

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SPHDataComplex mu(
			const SPHDataComplex &i_sph_data
	)	const
	{
		SPHConfig *sphConfig = i_sph_data.sphConfig;
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data = SPHDataComplex(sphConfig);

#pragma omp parallel for
		for (int n = 0; n <= i_sph_data.sphConfig->spec_n_max; n++)
		{
			int idx = i_sph_data.sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out_sph_data.data_spec[idx] =
					R(n-1,m)*i_sph_data.spec_get(n-1, m)
					+ S(n+1,m)*i_sph_data.spec_get(n+1, m);

				idx++;
			}
		}

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}

	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SPHDataComplex mu2(
			const SPHDataComplex &i_sph_data
	)	const
	{
		SPHConfig *sphConfig = i_sph_data.sphConfig;
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data = SPHDataComplex(sphConfig);

#pragma omp parallel for
		for (int n = 0; n <= i_sph_data.sphConfig->spec_n_max; n++)
		{
			int idx = i_sph_data.sphConfig->getArrayIndexByModes_Complex(n, -n);
			for (int m = -n; m <= n; m++)
			{
				out_sph_data.data_spec[idx] =
						+A(n-2,m)*i_sph_data.spec_get(n-2, m)
						+B(n+0,m)*i_sph_data.spec_get(n+0, m)
						+C(n+2,m)*i_sph_data.spec_get(n+2, m)
						;
				idx++;
			}
		}

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	/**
	 * Compute gradient component along latitude
	 */
	SPHDataComplex grad_lat(
			const SPHDataComplex &i_sph_data
	)	const
	{
		/*
		 * compute sin(theta)*d/d theta
		 * theta is the colatitude
		 *
		 * Hence, we have to
		 * 	first divide by sin(M_PI*0.5-phi) and
		 * 	second multiply by sqrt(1-mu*mu)
		 */
//		SPHDataComplex out_sph_data = spec_sinD(i_sph_data);
		SPHDataComplex out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.request_data_spatial();
		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, std::complex<double> &o_data)
				{
					double phi = asin(mu);

					//o_data /= sin(M_PI*0.5-phi);
					//o_data /= ::cos(phi);
					o_data /= sqrt(1.0-mu*mu);
				}
		);

#if 0
		/**
		 * WARNING: Leave this code here
		 * We can see that the following operations would cancel out.
		 * Therefore this was commented.
		 */
		// undo the sin(theta) and multiply with sqrt(1-mu*mu)
		out_sph_data.request_data_spatial();
		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					double phi = asin(mu);

					//o_data /= sin(M_PI*0.5-phi);
					o_data /= ::cos(phi);

					double cos_phi = std::sqrt((double)(1.0-mu*mu));
					o_data *= cos_phi;
				}
			);
#endif

		return out_sph_data;
	}



	/**
	 * Divergence Operator along longitude
	 *
	 * Identical to gradient operator along longitude
	 */
	SPHDataComplex div_lon(
			const SPHDataComplex &i_sph_data
	)	const
	{
		return grad_lon(i_sph_data);
	}




	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	SPHDataComplex div_lat(
			const SPHDataComplex &i_sph_data
	)	const
	{
		SPHDataComplex out_sph_data(i_sph_data);

		out_sph_data.spat_update_lambda_cogaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					//o_data *= cos(phi);
					o_data *= mu;
				}
			);

		// grad_lat = diff_lat_phi
		out_sph_data = grad_lat(out_sph_data);

		// undo the sin(theta) which is cos(phi)
#if 0
		out_sph_data.spat_update_lambda_cogaussian_grid(
				[](double lambda, double mu, std::complex<double> &o_data)
				{
					o_data /= mu;
					//o_data /= cos(phi);
				}
			);
#else
		out_sph_data.spat_update_lambda(
				[](double lambda, double phi, std::complex<double> &o_data)
				{
					//o_data /= mu;
					o_data /= cos(phi);
				}
			);
#endif
		return out_sph_data;
	}



	/**
	 * Laplace operator
	 */
	SPHDataComplex laplace(
			const SPHDataComplex &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data(i_sph_data);

		out_sph_data.spec_update_lambda(
				[](int n, int m, std::complex<double> &o_data)
				{
					o_data *= -(double)n*((double)n+1.0);
				}
			);

		return out_sph_data;
	}


public:
	/**
	 * Compute vorticity
	 *
	 * \eta = grad_lat(V_lon) - grad_lon(V_lat)
	 */
	SPHDataComplex vort(
			const SPHDataComplex &i_lon,
			const SPHDataComplex &i_lat
	)	const
	{
		return grad_lon(i_lat) - grad_lat(i_lon);
	}



public:
	/**
	 * Compute divergence
	 *
	 * \delta = div_lon(i_lon) + div_lan(i_lan)
	 */
	SPHDataComplex div(
			const SPHDataComplex &i_lon,
			const SPHDataComplex &i_lat
	)	const
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}
};






#endif /* SPHOPERATORS_HPP_ */
