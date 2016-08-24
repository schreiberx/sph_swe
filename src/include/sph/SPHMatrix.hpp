/*
 * SPHMatrix.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SPH_SPHMATRIX_HPP_
#define SRC_INCLUDE_SPH_SPHMATRIX_HPP_

#include <sph/SPHConfig.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <shtns.h>


/**
 * Matrix to store coefficients related to spherical harmonics
 */
template<typename T = double>
class SPHMatrix
{
	/**
	 * Data stores the diagonal and off-diagonal components for a matrix which is to be inverted.
	 * This matrix is partitioned by independent chunks for each mode m.
	 *
	 * The off-diagonal components connect the l-modes only!
	 *
	 * EXAMPLE:
	 *
	 *   Multiplying the data matrix with a vector P would require P to be in the following format:
	 *   P is given by P_n^m
	 *
	 * The the Vector P is given by
	 *   P = (P_0^0, P_1^0, P_2^0, P_3^0, P_1^1, P_2^1, P_3^1, P_2^2, P_3^2, P_3^3)^T
	 */

	T *data;
	int halosize_off_diagonal;
	int number_diagonals;

	SPHConfig *sphConfig;


public:
	SPHMatrix()	:
		data(nullptr),
		halosize_off_diagonal(-1),
		number_diagonals(-1),
		sphConfig(nullptr)
	{
	}



	/**
	 * Zero all matrix coefficients
	 */
	void zero()
	{
		for (std::size_t i = 0; i < sphConfig->spec_num_elems*number_diagonals; i++)
			data[i] = (T)0.0;
	}


	/**
	 * Setup data storage
	 */
	void setup(
			SPHConfig *i_sphConfig,				///< Handler to sphConfig
			int i_halosize_offdiagonal = 0		///< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		assert(data == nullptr);

		sphConfig = i_sphConfig;

		halosize_off_diagonal = i_halosize_offdiagonal;
		number_diagonals = 2*halosize_off_diagonal+1;

		data = MemBlockAlloc::alloc<T>(sizeof(T)*sphConfig->spec_num_elems*number_diagonals);
	}


	/**
	 * Return matrix row  which is related to the specified modes
	 */
	T *getMatrixRow(int in, int im)
	{

		sphConfig->shtns
	}

	void shutdown()
	{
		if (data != nullptr)
		{
			MemBlockAlloc::free(data, sizeof(T)*sphConfig->spec_num_elems*number_diagonals);
			data = nullptr;
		}
	}

	~SPHMatrix()
	{
		shutdown();
	}

};


#endif /* SRC_INCLUDE_SPH_SPHMATRIX_HPP_ */
