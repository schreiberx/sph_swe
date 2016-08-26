/*
 * SPHDataComplex.hpp
 *
 *  Created on: 9 Aug 2016
 *      Author: martin
 */

#ifndef SPHDATACOMPLEX_HPP_
#define SPHDATACOMPLEX_HPP_

#include <complex.h>
#include <fftw3.h>
#include <functional>
#include <array>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "SPHConfig.hpp"
#include "SPHData.hpp"



class SPHDataComplex
{
public:
	SPHConfig *sphConfig;

public:
	std::complex<double> *data_spat;
	std::complex<double> *data_spec;

	bool data_spec_valid;
	bool data_spat_valid;

public:
	SPHDataComplex(
			SPHConfig *i_sphConfig
	)	:
		sphConfig(i_sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		assert(i_sphConfig != 0);

		setup();
	}


public:
	SPHDataComplex(
			const SPHDataComplex &i_sph_data
	)	:
		sphConfig(i_sph_data.sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		setup();

		operator=(i_sph_data);
	}

public:
	SPHDataComplex(
			const SPHData &i_sph_data
	)	:
		sphConfig(i_sph_data.sphConfig),
		data_spat(nullptr),
		data_spec(nullptr)
	{
		setup();

		operator=(i_sph_data);
	}


public:
	SPHDataComplex& operator=(
			const SPHDataComplex &i_sph_data
	)
	{
		data_spec_valid = i_sph_data.data_spec_valid;
		data_spat_valid = i_sph_data.data_spat_valid;

		if (data_spat_valid)
			memcpy(data_spat, i_sph_data.data_spat, sizeof(cplx)*sphConfig->spat_num_elems);

		if (data_spec_valid)
			memcpy(data_spec, i_sph_data.data_spec, sizeof(cplx)*sphConfig->cplx_spec_num_elems);

		return *this;
	}


public:
	SPHDataComplex& operator=(
			const SPHData &i_sph_data
	)
	{
		if (i_sph_data.data_spat_valid)
		{
			for (std::size_t i = 0; i < sphConfig->spat_num_elems; i++)
				data_spat[i] = {i_sph_data.data_spat[i], 0};
		}

		if (i_sph_data.data_spec_valid)
			memcpy(data_spec, i_sph_data.data_spec, sizeof(cplx)*sphConfig->cplx_spec_num_elems);

		data_spec_valid = i_sph_data.data_spec_valid;
		data_spat_valid = i_sph_data.data_spat_valid;

		return *this;
	}


	void request_data_spectral()	const
	{
		if (data_spec_valid)
			return;

		assert(data_spat_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		spat_cplx_to_SH(sphConfig->shtns, data_spat, data_spec);

		SPHDataComplex *this_var = (SPHDataComplex*)this;

		this_var->data_spat_valid = false;
		this_var->data_spec_valid = true;
	}


	void request_data_spatial()	const
	{
		if (data_spat_valid)
			return;

		assert(data_spec_valid);

		/**
		 * Warning: This is an in-situ operation.
		 * Therefore, the data in the source array will be destroyed.
		 */
		SH_to_spat_cplx(sphConfig->shtns, data_spec, data_spat);

		SPHDataComplex *this_var = (SPHDataComplex*)this;

		this_var->data_spec_valid = false;
		this_var->data_spat_valid = true;
	}


	SPHDataComplex operator+(
			const SPHDataComplex &i_sph_data
	)
	{
		request_data_spectral();
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] + i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SPHDataComplex& operator+=(
			const SPHDataComplex &i_sph_data
	)
	{
		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] += i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}


	SPHDataComplex& operator-=(
			const SPHDataComplex &i_sph_data
	)
	{
		request_data_spectral();
		i_sph_data.request_data_spectral();

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			data_spec[idx] -= i_sph_data.data_spec[idx];

		data_spec_valid = true;
		data_spat_valid = false;

		return *this;
	}



	SPHDataComplex operator-(
			const SPHDataComplex &i_sph_data
	)
	{
		request_data_spectral();
		i_sph_data.request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx] - i_sph_data.data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SPHDataComplex operator-()
	{
		request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = -data_spec[idx];

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


	SPHDataComplex operator*(
			const SPHDataComplex &i_sph_data
	)	const
	{
		request_data_spatial();
		i_sph_data.request_data_spatial();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_elems; i++)
			out_sph_data.data_spat[i] = i_sph_data.data_spat[i]*data_spat[i];

		out_sph_data.data_spec_valid = false;
		out_sph_data.data_spat_valid = true;

		return out_sph_data;
	}


	SPHDataComplex operator*(
			const double i_value
	)	const
	{
		request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]*i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	SPHDataComplex operator+(
			const double i_value
	)	const
	{
		SPHDataComplex out_sph_data(*this);

		out_sph_data.request_data_spectral();

		/**
		 * TODO: Check normalization
		 */
		out_sph_data.data_spec[0] += i_value;

		return out_sph_data;
	}



	SPHDataComplex operator+(
			const std::complex<double> &i_value
	)	const
	{
		SPHDataComplex out_sph_data(*this);

		out_sph_data.request_data_spectral();
		out_sph_data.data_spec[0] += i_value*std::sqrt(4.0*M_PI);

		return out_sph_data;
	}



	SPHDataComplex operator*(
			const std::complex<double> &i_value
	)	const
	{
		request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]*i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}



	SPHDataComplex operator/(
			double i_value
	)	const
	{
		request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);

#pragma omp parallel for
		for (int idx = 0; idx < sphConfig->spec_num_elems; idx++)
			out_sph_data.data_spec[idx] = data_spec[idx]/i_value;

		out_sph_data.data_spec_valid = true;
		out_sph_data.data_spat_valid = false;

		return out_sph_data;
	}


private:
	void setup()
	{
		data_spec_valid = false;
		data_spat_valid = false;

		data_spat = (std::complex<double> *) fftw_malloc(sphConfig->spat_num_elems * sizeof(cplx));
		data_spec = (std::complex<double> *) fftw_malloc(sphConfig->cplx_spec_num_elems * sizeof(cplx));
	}

#if 0
public:
	/**
	 * Truncate modes which are not representable in spectral space
	 */
	SPHDataComplex spat_truncate()
	{
		request_data_spatial();

		SPHDataComplex out_sph_data(sphConfig);
		spat_cplx_to_SH(sphConfig->shtns, data_spat, out_sph_data.data_spec);
		SH_to_spat_cplx(sphConfig->shtns, out_sph_data.data_spec, out_sph_data.data_spat);

		out_sph_data.data_spat_valid = true;
		out_sph_data.data_spec_valid = false;

		return out_sph_data;
	}
#endif

	/**
	 * Truncate modes which are not representable in spectral space
	 */
public:
	SPHDataComplex spec_truncate()
	{
		request_data_spectral();

		SPHDataComplex out_sph_data(sphConfig);
		SH_to_spat_cplx(sphConfig->shtns, data_spec, out_sph_data.data_spat);
		spat_cplx_to_SH(sphConfig->shtns, out_sph_data.data_spat, out_sph_data.data_spec);

		out_sph_data.data_spat_valid = false;
		out_sph_data.data_spec_valid = true;

		return out_sph_data;
	}

#if 0
	void spec_set_lambda(
			std::function<void(int,int,cplx&)> i_lambda
	)
	{
		if (data_spat_valid)
			request_data_spectral();

#pragma omp parallel for
		for (int m = 0; m <= spec_m_max; m++)
		{
			for (int n = m; n <= spec_n_max; n++)
			{
				/*
				 * TODO: Iteration range might be different!!!!!
				 * TODO: Check this before uncommenting this part
				 */
				assert(LM(sphConfig->shtns, n, m) < sphConfig->spec_num_elems);
				i_lambda(m, n, data_spec[LM(sphConfig->shtns, n, m)]);
			}
		}

		data_spat_valid = false;
		data_spec_valid = true;
	}
#endif



	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-M_PI/2;M_PI/2])
	 */

public:
	void spat_update_lambda(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_spatial();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				//double colatitude = acos(shtns->ct[j]);

				/*
				 * Colatitude is 0 at the north pole and 180 at the south pole
				 *
				 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
				 */
				//double lat_degree = M_PI*0.5 - colatitude;
				double lat_degree = sphConfig->lat[j];

				i_lambda(lon_degree, lat_degree, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */
	void spat_update_lambda_gaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_spatial();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				double mu = sphConfig->lat_gaussian[j];

				i_lambda(lon_degree, mu, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}




	/*
	 * Set values for all latitude and longitude degrees
	 *
	 * lambda function parameters: (longitude \in [0;2*pi], Gaussian latitude \in [-1;1])
	 */

	void spat_update_lambda_cogaussian_grid(
			std::function<void(double,double,cplx&)> i_lambda	///< lambda function to return value for lat/mu
	)
	{
		if (data_spec_valid)
			request_data_spatial();

#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
		{
			double lon_degree = (((double)i)/(double)sphConfig->spat_num_lon)*2.0*M_PI;

			for (int j = 0; j < sphConfig->spat_num_lat; j++)
			{
				double comu = sphConfig->shtns->st[j];
				/*
				 * IDENTITAL FORMULATION
				double mu = shtns->ct[j];
				double comu = sqrt(1.0-mu*mu);
				*/

				i_lambda(lon_degree, comu, data_spat[i*sphConfig->spat_num_lat + j]);
			}
		}

		data_spat_valid = true;
		data_spec_valid = false;
	}



	/*
	 * Set all values to zero
	 */
	void spat_set_zero()
	{
#pragma omp parallel for
		for (int i = 0; i < sphConfig->spat_num_lon; i++)
			for (int j = 0; j < sphConfig->spat_num_lat; j++)
				data_spat[i*sphConfig->spat_num_lat + j] = 0;

		data_spat_valid = true;
		data_spec_valid = false;
	}


	/**
	 * Return the maximum error norm
	 */
	double spat_reduce_error_max(
			const SPHDataComplex &i_sph_data
	)
	{
#if 0
		SPHDataComplex tmp = *this-i_sph_data;

		double error = -1;

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(
						error,
						std::abs(tmp.data_spat[j])
						);
		}
		return error;

#else

		request_data_spatial();
		i_sph_data.request_data_spatial();

		double error = -1;

		for (int j = 0; j < sphConfig->spat_num_elems; j++)
		{
			error = std::max(
						error,
						std::abs(
								data_spat[j].real() - i_sph_data.data_spat[j].real()
							)
						);

			error = std::max(
						error,
						std::abs(
								data_spat[j].imag() - i_sph_data.data_spat[j].imag()
							)
						);
		}
		return error;
#endif
	}


public:
	~SPHDataComplex()
	{
		fftw_free(data_spat);

//		std::size_t spec_num_elems = sphConfig->spec_n_max*(sphConfig->spec_n_max+1)+sphConfig->spec_m_max;
		fftw_free(data_spec);
	}



	void spat_write_file(
			const char *i_filename,
			const char *i_title = "",
			int i_precision = 8
	)
	{
		request_data_spatial();

		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (std::size_t i = 0; i < sphConfig->spat_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			lon_degree = lon_degree/M_PI*180.0;

			file << lon_degree;
			if (i < sphConfig->spat_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphConfig->spat_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphConfig->spat_num_lon; i++)
        		{
        			file << data_spat[i*sphConfig->spat_num_lat+j];
        			if (i < sphConfig->spat_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}

	void spat_write_file_lon_pi_shifted(
			const char *i_filename,
			std::string i_title = "",
			int i_precision = 8
	)
	{
		request_data_spatial();

		std::ofstream file(i_filename, std::ios_base::trunc);

		file << std::setprecision(i_precision);
		file << "#TI " << i_title << std::endl;
		file << "#TX Longitude" << std::endl;
		file << "#TY Latitude" << std::endl;

		//file << "lat\\lon\t";
		// Use 0 to make it processable by python
		file << "0\t";

		for (std::size_t i = 0; i < sphConfig->spat_num_lon; i++)
		{
//			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			double lon_degree = ((double)i/(double)sphConfig->spat_num_lon)*2.0*M_PI;
			lon_degree = (lon_degree-M_PI)/M_PI*180.0;

			file << lon_degree;
			if (i < sphConfig->spat_num_lon-1)
				file << "\t";
		}
		file << std::endl;

        for (int j = sphConfig->spat_num_lat-1; j >= 0; j--)
        {
//        		double lat_degree =  M_PI*0.5 - acos(shtns->ct[j]);
        		double lat_degree = sphConfig->lat[j];
        		lat_degree = lat_degree/M_PI*180.0;

        		file << lat_degree << "\t";

        		for (int i = 0; i < sphConfig->spat_num_lon; i++)
        		{
        			int ia = i+sphConfig->spat_num_lon/2;
        			if (ia >= sphConfig->spat_num_lon)
        				ia -= sphConfig->spat_num_lon;

        			file << data_spat[ia*sphConfig->spat_num_lat+j];
        			if (i < sphConfig->spat_num_lon-1)
        				file << "\t";
        		}
        		file << std::endl;
        }
        file.close();
	}
};




/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SPHDataComplex operator*(
		const double i_value,
		const SPHDataComplex &i_array_data
)
{
	return ((SPHDataComplex&)i_array_data)*i_value;
}


/**
 * operator to support operations such as:
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
SPHDataComplex operator*(
		const std::complex<double> &i_value,
		const SPHDataComplex &i_array_data
)
{
	return i_array_data*i_value;
}


/**
 * operator to support operations such as:
 *
 * 1.5 - arrayData;
 *
 * Otherwise, we'd have to write it as arrayData-1.5
 *
 */
#if 0
inline
static
SPHDataComplex operator-(
		const double i_value,
		const SPHDataComplex &i_array_data
)
{
	return ((SPHDataComplex&)i_array_data).valueMinusThis(i_value);
//	return -(((SPHDataComplex&)i_array_data).operator-(i_value));
}
#endif
/**
 * operator to support operations such as:
 *
 * 1.5 + arrayData;
 *
 * Otherwise, we'd have to write it as arrayData+1.5
 *
 */

inline
static
SPHDataComplex operator+(
		const double i_value,
		const SPHDataComplex &i_array_data
)
{
	return ((SPHDataComplex&)i_array_data)+i_value;
}

inline
static
SPHDataComplex operator+(
		const std::complex<double> &i_value,
		const SPHDataComplex &i_array_data
)
{
	return i_array_data+i_value;
}

#endif /* SPHDATA_HPP_ */
