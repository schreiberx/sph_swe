/*
 * SPHOperators.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: martin
 */

#ifndef SPHOPERATORS_HPP_
#define SPHOPERATORS_HPP_

#include <sph/SPHData.hpp>
#include <sph/SPHIdentities.hpp>

class SPHOperators	:
		public SPHIdentities
{
	double *sin_mx = nullptr;
	double *cos_mx = nullptr;

	friend SPHConfig;

public:
	void setup(
			SPHConfig *i_sphConfig
	)
	{
		// prepare matrix which represents applying sin(theta)*d/dtheta()
		sin_mx = (double*)fftw_malloc(i_sphConfig->spec_num_elems * 2 * sizeof(double));
		st_dt_matrix(i_sphConfig->shtns, sin_mx);

		/*
		 * SHTNS uses the co-latitude.
		 * Hence, we have to flip the sign.
		 */
		for (std::size_t i = 0; i < i_sphConfig->spec_num_elems * 2; i++)
			sin_mx[i] *= -1.0;

		cos_mx = (double*)fftw_malloc(i_sphConfig->spec_num_elems * 2 * sizeof(double));
		mul_ct_matrix(i_sphConfig->shtns, cos_mx);
		for (std::size_t i = 0; i < i_sphConfig->spec_num_elems * 2; i++)
			cos_mx[i] *= -1.0;
	}

	~SPHOperators()
	{
		fftw_free(sin_mx);
		fftw_free(cos_mx);
	}


#if 0
public:
	/**
	 * Compute the coriolis effect on the given velocity
	 *
	 * multiply (n=1, m=0) with \Omega sqrt(8/3)
	 */
	SPHData coriolis(
			const SPHData &i_sph_data,
			double i_coriolis_frequency
	)	const
	{
		// THIS IS NOT WORKING
		assert(false);
		SPHData out_sph_data(i_sph_data);

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
	SPHData diff_lon(
			const SPHData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHData out_sph_data(i_sph_data.sphConfig);

		// compute d/dlambda in spectral space
		int idx = 0;
		for (int m = 0; m <= i_sph_data.sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= i_sph_data.sphConfig->spec_n_max; n++)
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
	SPHData diff_lat_mu(
			const SPHData &i_sph_data
	)	const
	{

#if 1
		SPHData out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					o_data /= (1.0-mu*mu);
				}
		);
#else
		SPHData out_sph_data = grad_lat(i_sph_data);

		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					o_data /= sqrt(1.0-mu*mu);
				}
		);
#endif
		return out_sph_data;
	}



	/**
	 * Compute differential along latitude
	 *
	 * Compute d/d phi f(lambda,mu)
	 */
	SPHData diff_lat_phi(
			const SPHData &i_sph_data
	)	const
	{
		return grad_lat(i_sph_data);
	}



	/**
	 * Compute gradient component along longitude (lambda)
	 */
	SPHData grad_lon(
			const SPHData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHData out_sph_data = diff_lon(i_sph_data);

		out_sph_data.request_data_spatial();

		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
				{
					double cos_phi = sqrt(1.0-mu*mu);
					o_data /= cos_phi;
				}
		);

		return out_sph_data;
	}



	static double fac_double(double v)
	{
		double retval = 1;
	    for (double i = 1; i <= v; ++i)
	    	retval *= i;
	    return retval;
	}



	SPHData spec_one_minus_mu_squared_diff_lat_mu(
			const SPHData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();
		SPHConfig *sphConfig = i_sph_data.sphConfig;

		SPHData out_sph_data = SPHData(sphConfig);

#pragma omp parallel for
		for (std::size_t idx = 0; idx < sphConfig->spec_num_elems; idx++)
		{
			// modes
			int ni = sphConfig->shtns->li[idx];
			int mi = sphConfig->shtns->mi[idx];

			out_sph_data.data_spec[idx] =	((-ni+1.0)*R(ni-1,mi))*i_sph_data.spec_get(ni-1, mi) +
											((ni+2.0)*S(ni+1,mi))*i_sph_data.spec_get(ni+1, mi);
		}

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	/**
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
	SPHData mu(
			const SPHData &i_sph_data
	)	const
	{
		SPHConfig *sphConfig = i_sph_data.sphConfig;
		i_sph_data.request_data_spectral();

		SPHData out_sph_data = SPHData(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
		{
			// modes
			int ni = sphConfig->shtns->li[idx];
			int mi = sphConfig->shtns->mi[idx];

			std::complex<double> P0, P2;
			i_sph_data.spec_getElement_im_in(ni-1, mi, P0);
			i_sph_data.spec_getElement_im_in(ni+1, mi, P2);

			std::complex<double> P = R(ni-1,mi)*P0 + S(ni+1,mi)*P2;

			out_sph_data.data_spec[idx] = P;
		}

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	/**
	 * Compute gradient component along latitude
	 */
	SPHData grad_lat(
			const SPHData &i_sph_data
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
//		SPHData out_sph_data = spec_sinD(i_sph_data);
		SPHData out_sph_data = spec_one_minus_mu_squared_diff_lat_mu(i_sph_data);

		out_sph_data.request_data_spatial();
		out_sph_data.spat_update_lambda_gaussian_grid(
				[this](double lambda, double mu, double &o_data)
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
	SPHData div_lon(
			const SPHData &i_sph_data
	)	const
	{
		return grad_lon(i_sph_data);
	}


#if 0
	SPHData spec_sinD(
			const SPHData &i_sph_data
	)	const
	{
		assert(sin_mx != nullptr);

		i_sph_data.request_data_spectral();

		SPHData out_sph_data(i_sph_data.sphConfig);
		SH_mul_mx(out_sph_data.sphConfig->shtns, sin_mx, i_sph_data.data_spec, out_sph_data.data_spec);

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	SPHData spec_cos(
			const SPHData &i_sph_data
	)	const
	{
		assert(sin_mx != nullptr);

		i_sph_data.request_data_spectral();

		SPHData out_sph_data(i_sph_data.sphConfig);
		SH_mul_mx(out_sph_data.sphConfig->shtns, cos_mx, i_sph_data.data_spec, out_sph_data.data_spec);

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}
#endif



	/**
	 * Divergence Operator along latitude
	 *
	 * d(sqrt(1-mu*mu)*F)/dmu
	 */
	SPHData div_lat(
			const SPHData &i_sph_data
	)	const
	{
		assert(sin_mx != nullptr);

		SPHData out_sph_data(i_sph_data);

		out_sph_data.request_data_spatial();
		out_sph_data.spat_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					//o_data *= cos(phi);
					o_data *= mu;
				}
			);

		// grad_lat = diff_lat_phi
		out_sph_data = grad_lat(out_sph_data);

		// undo the sin(theta) which is cos(phi
		out_sph_data.request_data_spatial();
		out_sph_data.spat_update_lambda_cogaussian_grid(
				[](double lambda, double mu, double &o_data)
				{
					o_data /= mu;
					//o_data /= cos(phi);
				}
			);

		return out_sph_data;
	}



	/**
	 * Laplace operator
	 */
	SPHData laplace(
			const SPHData &i_sph_data
	)	const
	{
		i_sph_data.request_data_spectral();

		SPHData out_sph_data(i_sph_data);

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
	SPHData vort(
			const SPHData &i_lon,
			const SPHData &i_lat
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
	SPHData div(
			const SPHData &i_lon,
			const SPHData &i_lat
	)	const
	{
		return div_lon(i_lon) + div_lat(i_lat);
	}
};






#endif /* SPHOPERATORS_HPP_ */
