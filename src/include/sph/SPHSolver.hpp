/*
 * SPHSolver.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SPH_SPHSOLVER_HPP_
#define SRC_INCLUDE_SPH_SPHSOLVER_HPP_

#include <libmath/BandedMatrixSolver.hpp>
#include <sph/SPHMatrix.hpp>
#include <sph/SPHIdentities.hpp>



/**
 *
 * phi(lambda,mu) denotes the solution
 */
template <typename T>
class SPHSolver	:
		SPHIdentities
{
public:
	/**
	 * Matrix on left-hand side
	 */
	SPHMatrix<T> lhs;

	/**
	 * SPH configuration
	 */
	SPHConfig *sphConfig;

	/**
	 * Solver for banded matrix
	 */
	BandedMatrixSolver< std::complex<double> > bandedMatrixSolver;


	/**
	 * Setup the SPH solver
	 */
public:
	void setup(
			SPHConfig *i_sphConfig,		///< Handler to sphConfig
			int i_halosize_offdiagonal	///< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		sphConfig = i_sphConfig;

		lhs.setup(sphConfig, i_halosize_offdiagonal);

		bandedMatrixSolver.setup(i_sphConfig->spec_n_max+1, i_halosize_offdiagonal);
	}



	/**
	 * Solver for
	 * 	a*phi(lambda,mu)
	 */
	void solver_component_scalar_phi(
			const std::complex<double> &i_value
	)
	{
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, 0, i_value);
			}
		}
	}



	/**
	 * Solver for
	 * 	mu*phi(lambda,mu)
	 */
	void solver_component_mu_phi()
	{
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, R(n-1,m));
				lhs.rowElement_add(row, n, m, +1, S(n+1,m));
			}
		}
	}




	/**
	 * Solver for
	 * 	mu*mu*phi(lambda,mu)
	 */
	void solver_component_mu_mu_phi()
	{
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -2, R(n-1,m)*R(n-2,m));
				lhs.rowElement_add(row, n, m,  0, R(n-1,m)*S(n,m) + S(n+1,m)*R(n,m));
				lhs.rowElement_add(row, n, m, +2, S(n+1,m)*S(n+2,m));
			}
		}
	}



	/**
	 * Solver for
	 * 	(1-mu*mu)*d/dmu phi(lambda,mu)
	 */
	void solver_component_one_minus_mu_mu_diff_mu_phi()
	{
		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, (-(double)n+1.0)*R(n-1,m));
				lhs.rowElement_add(row, n, m, +1, ((double)n+2.0)*S(n+1,m));
			}
		}
//		lhs.print();
	}


	/**
	 * Apply the solver matrix.
	 * This function is intended to be used for debugging.
	 * WARNING: This only multiplies the i_x values with the matrix.
	 * Use solve(...) to solve for the matrix
	 */
	SPHData apply(
			const SPHData &i_x	///< solution to be searched
	)
	{
		SPHData out(sphConfig);

		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			std::size_t idx = sphConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphConfig->spec_n_max; n++)
			{
				out.data_spec[idx] = 0;

				std::complex<double> *row = lhs.getMatrixRow(n, m);
				for (int i = 0; i < lhs.num_diagonals; i++)
				{
					int delta = i-lhs.halosize_off_diagonal;
					out.data_spec[idx] += lhs.rowElement_getRef(row, n, m, delta)*i_x.spec_get(n+delta, m);
				}

				idx++;
			}
		}

		out.data_spat_valid = false;
		out.data_spec_valid = true;

		return out;
	}


	SPHData solve(
			const SPHData &i_rhs
	)
	{
		i_rhs.request_data_spectral();

		SPHData out(sphConfig);

		for (int m = 0; m <= sphConfig->spec_m_max; m++)
		{
			int idx = sphConfig->getArrayIndexByModes(m,m);

			bandedMatrixSolver.solve_diagBandedInverse_Carray(
							&lhs.data[idx*lhs.num_diagonals],
							&i_rhs.data_spec[idx],
							&out.data_spec[idx],
							sphConfig->spec_n_max+1-m	// size of block
					);
		}

		out.data_spec_valid = true;
		out.data_spat_valid = false;

		return out;
	}
};


#endif /* SRC_INCLUDE_SPH_SPHSOLVER_HPP_ */
