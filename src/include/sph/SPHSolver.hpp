/*
 * SPHSolver.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SPH_SPHSOLVER_HPP_
#define SRC_INCLUDE_SPH_SPHSOLVER_HPP_

#include <sph/SPHMatrix.hpp>


template <typename T>
class SPHSolver
{
	SPHMatrix<T> lhs;
	SPHMatrix<T> inv_lhs;

	bool inv_lhs_valid;

	SPHConfig *sphConfig;

	/**
	 * Setup the SPH solver
	 */
public:
	void setup(
			SPHConfig *i_sphConfig,				///< Handler to sphConfig
			int i_halosize_offdiagonal = -1		///< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		sphConfig = i_sphConfig;

		if (i_halosize_offdiagonal == -1)
			i_halosize_offdiagonal = 4;

		lhs.setup(sphConfig, i_halosize_offdiagonal);
		inv_lhs.setup(sphConfig, i_halosize_offdiagonal);

		inv_lhs_valid = false;
	}



	/**
	 * Add a solver component for multiplication with a constant
	 */
	void solver_component_const(
			T i_value
	)
	{
		T set_value = i_value;

		int idx = 0;
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.setRowElement(row, n, m, 0, set_value);
			}
		}
	}


	void setupInverse()
	{
		/*
		 * TODO: compute inverse
		 */
		inv_lhs_valid = true;
	}



	SPHData solve(SPHData &i_sph_data)
	{
		assert(inv_lhs_valid);

		SPHData out(sphConfig);


		int idx = 0;
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			/*
			 * TODO: matrix multiplication
			 */
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				for (int i = 0; i < inv_lhs.matrix_width; i++)
				{
					//
				}
				idx++;
			}
		}

		return out;
	}
};


#endif /* SRC_INCLUDE_SPH_SPHSOLVER_HPP_ */
