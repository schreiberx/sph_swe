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
	void zeroAll()
	{
		for (std::size_t i = 0; i < sphConfig->spec_num_elems*number_diagonals; i++)
			data[i] = T(0);
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

		zeroAll();
	}



	/**
	 * Return matrix row  which is related to the specified modes
	 */
	T *getMatrixRow(
			int n,		///< row related to P Legendre mode n
			int m		///< row related to P Fourier mode n
	)
	{
		std::size_t idx = sphConfig->getPIndexByModes(n, m);
		return data+idx*number_diagonals;
	}



	/**
	 * Set an element in the row to the specified value
	 */
	void setRowElement(
			T *io_row,		///< pointer to current row
			int i_row_n,	///< row related to P Legendre mode n
			int i_row_m,	///< row related to P Fourier mode n
			int rel_n,		///< Relative Legendre mode n (e.g. -1 or +2)
			T &i_value		///< Value to set in the matrix element
	)
	{
		assert(i_row_n > i_row_m);
		assert(i_row_m >= 0);
		assert(i_row_m <= sphConfig->spec_m_max);

		int n = i_row_n+rel_n;

		if (n < 0 || n < i_row_m || n > sphConfig->spec_n_max)
			return;

		int idx = rel_n + halosize_off_diagonal+1;

		assert(idx >= 0 || idx <= halosize_off_diagonal*2+1);

		io_row[idx] = i_value;
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
