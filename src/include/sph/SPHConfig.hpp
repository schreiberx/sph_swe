/*
 * SPHSetup.hpp
 *
 *  Created on: 12 Aug 2016
 *      Author: martin
 */

#ifndef SPHSETUP_HPP_
#define SPHSETUP_HPP_


// undef for C++11 support
#undef _COMPLEX_H
#include <shtns.h>
#include <cmath>
#include <fftw3.h>
#include <iostream>


class SPHConfig
{
	friend class SPHOperators;
	friend class SPHData;
	friend class SPHDataComplex;

private:
	shtns_cfg shtns;

	/**
	 * Number of longitudes
	 */
public:
	int spat_num_lon;

	/**
	 * Number of latitudes
	 */
public:
	int spat_num_lat;

	/**
	 * Number of total longitudes and latitudes
	 */
public:
	int spat_num_elems;


	/**
	 * Number of modes
	 */
public:
	int spec_n_max;
	int spec_m_max;

	/**
	 * Number of total mode variations
	 */
	int spec_num_elems;

	/**
	 * Number of elements for complex-valued spatial data
	 */
	int cplx_spec_num_elems;


	/**
	 * Array with latitude phi angles
	 *
	 * WARNING: Phi is not the phi from SHTNS
	 */
public:
	double *lat;

	/**
	 * Array with mu = cos(phi) values
	 */
public:
	double *lat_gaussian;

public:
	SPHConfig()	:
		shtns(nullptr),
		spat_num_lat(-1),
		spat_num_lon(-1),
		spat_num_elems(-1),

		spec_n_max(-1),
		spec_m_max(-1),
		spec_num_elems(-1),
		cplx_spec_num_elems(-1),

		lat(nullptr),
		lat_gaussian(nullptr)
	{
	}


	std::size_t getPIndexByModes(int l, int im)
	{
//		return (spec_n_max-im)*im + ((im+1)*im)/2+l;
		return (im*(2*spec_n_max-im+1)>>1)+l;
	}


private:
	void setup_data()
	{
		spat_num_lat = shtns->nlat;
		spat_num_lon = shtns->nphi;
		spat_num_elems = shtns->nspat;

		spec_n_max = shtns->lmax;
		spec_m_max = shtns->mmax;
		spec_num_elems = shtns->nlm;
		cplx_spec_num_elems = (spec_n_max+1)*(spec_n_max+1);

		/**
		 * Some safety checks to make sure that we really get what we've asked for
		 */

		/**
		 * TEST: iteration over the modes n,m
		 */
		{
			int idx = 0;
			for (int m = 0; m <= spec_m_max; m++)
			{
				for (int n = m; n <= spec_n_max; n++)
				{
					int test_idx = getPIndexByModes(n,m);

					if (test_idx != idx)
					{
						std::cerr << "IDX TEST NOT SUCCESSFUL" << std::endl;
						std::cout << "n=" << n << ", m=" << m << "     " << idx << ", " << test_idx << std::endl;
						exit(1);
					}
					idx++;
				}
			}
		}

		lat = (double*)fftw_malloc(sizeof(double)*shtns->nlat);

		/*
		 * Colatitude is 0 at the north pole and 180 at the south pole
		 *
		 * WARNING: The latitude degrees are not equidistant spaced in the angles!!!! We have to use the shtns->ct lookup table
		 */
		for (int i = 0; i < shtns->nlat; i++)
			lat[i] = M_PI*0.5 - ::acos(shtns->ct[i]);

		lat_gaussian = (double*)fftw_malloc(sizeof(double)*shtns->nlat);
		for (int i = 0; i < shtns->nlat; i++)
			lat_gaussian[i] = shtns->ct[i];;
	}

public:
	void setup(
			int mmax,
			int nmax,
			int nphi,
			int nlat
	)
	{
		shtns_verbose(1);			// displays informations during initialization.
		shtns_use_threads(0);		// enable multi-threaded transforms (if supported).

		shtns = shtns_create(
				nmax,
				mmax,
				1,
				(shtns_norm)((int)sht_orthonormal /*| SHT_NO_CS_PHASE*/)
			);

		shtns_set_grid(
				shtns,
				// TODO: replace this with sht_gauss
				(shtns_type)(sht_quick_init | SHT_THETA_CONTIGUOUS),
				//sht_gauss | SHT_THETA_CONTIGUOUS,	// use gaussian grid
				0,
				nlat,		// number of latitude grid points
				nphi		// number of longitude grid points
			);

		setup_data();
	}


	/**
	 * Setup with given modes.
	 * Spatial resolution will be determined automatically
	 */
	void setup(
			int mmax,
			int nmax,
			int *nphi,
			int *nlat
	)
	{
		shtns_verbose(1);			// displays informations during initialization.
		shtns_use_threads(0);		// enable multi-threaded transforms (if supported).


		shtns = shtns_create(
				nmax,
				mmax,
				1,
				sht_orthonormal
			);

		*nphi = 0;
		*nlat = 0;

		shtns_set_grid_auto(
				shtns,
				// TODO: replace this with sht_gauss
				(shtns_type)(sht_quick_init | SHT_THETA_CONTIGUOUS),
				//sht_gauss | SHT_THETA_CONTIGUOUS,	// use gaussian grid
				0,
				2,		// use order 2
				nlat,
				nphi
			);

		setup_data();
	}


	void shutdown()
	{
		if (shtns != nullptr)
		{
			shtns_destroy(shtns);
			shtns = nullptr;
		}

		fftw_free(lat);
		lat = nullptr;

		fftw_free(lat_gaussian);
		lat_gaussian = nullptr;
	}

	~SPHConfig()
	{
		shutdown();
	}
};



#endif /* SPHSETUP_HPP_ */
