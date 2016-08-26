/*
 * DiagBandedMatrix.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_
#define SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_

#include <complex>
#include <string.h>
#include <stdlib.h>


/*
 * We can't store the inverse of the general band matrix since
 * the inverse would result in a dense matrix.
 *
 * http://www.netlib.org/lapack/explore-html/d9/dbb/group__complex16_g_bsolve_ga908abc0aad64131b9a32edb08510eb00.html#ga908abc0aad64131b9a32edb08510eb00
 *
 * ZGBSV computes the solution to system of linear equations A * X = B for GB matrices (simple driver)
 *
subroutine zgbsv
	(
		integer  	N,
		integer  	KL,
		integer  	KU,
		integer  	NRHS,
		complex*16, dimension( ldab, * )  	AB,
		integer  	LDAB,
		integer, dimension( * )  	IPIV,
		complex*16, dimension( ldb, * )  	B,
		integer  	LDB,
		integer  	INFO
	)
*/

/*
 * We create the interface to this function here
 */
extern "C"
{
	void zgbsv_(
			int &N,
			int &KL,
			int &KU,
			int &NRHS,
			std::complex<double> *AB,
			int &LDAB,
			int *IPIV,
			std::complex<double> *B,
			int &LDB,
			int &INFO
	);

	void zlapmt_(
			bool &forward,
			int &M,	// rows
			int &N,	// cols
			std::complex<double> *X,
			int &LDX,
			int *K
	);
}


template <typename T>
class MatrixSolversCommon
{
public:
	int N;
	int num_diagonals;
	int num_off_diagonals;

	int LDAB;

	std::complex<double>* AB;
	int *IPIV;

	MatrixSolversCommon()	:
		AB(nullptr)
	{
	}



	~MatrixSolversCommon()
	{
		shutdown();
	}

	void setup(
			int i_N,			///< size of the matrix
			int i_num_diagonals	///< number of block diagonals
	)
	{
		N = i_N;
		num_diagonals = i_num_diagonals;
		num_off_diagonals = num_diagonals>>1;
		assert(2*num_off_diagonals+1 == num_diagonals);

		LDAB = 2*num_off_diagonals + num_off_diagonals + 1;

		AB = (std::complex<double>*)malloc(LDAB*i_N);
		IPIV = (int*)malloc(sizeof(int)*i_N);
	}

	void shutdown()
	{
		free(IPIV);
		free(AB);
	}
};

template <typename T>
class BandedMatrixSolver	:
		public MatrixSolversCommon<T>
{
public:
	/**
	 * Solve a diagonal banded matrix
	 *
	 * A*X = B
	 */
	void solve_diagBandedInverse(
		const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
		const std::complex<double>* i_b,		///< RHS of equation and output of solution X
		std::complex<double>* o_x,
		int i_num_diagonals,
		int i_size
	);
};


template <>
class BandedMatrixSolver<std::complex<double>>	:
		public MatrixSolversCommon<std::complex<double>>
{
public:
	void solve_diagBandedInverse(
		const std::complex<double>* i_A,
		const std::complex<double>* i_b,
		std::complex<double>* o_x
	)
	{
		assert(num_diagonals & 1 == 1);
		assert(AB != nullptr);

#if 1
		/*
		 * Make a copy of the array data since this is a destructive function
		 */
		memcpy((void*)AB, (const void*)i_A, sizeof(std::complex<double>)*num_diagonals*N);
		memcpy((void*)o_x, (const void*)i_b, sizeof(std::complex<double>)*N);


		int one = 1;
		int info;
		zgbsv_(
				N,					// number of linear equations
				num_off_diagonals,	// number of subdiagonals
				num_off_diagonals,	// number of superdiagonals
				one,				// number of columns of matrix B
				AB,					// array with matrix A to solve for
				LDAB,				// leading dimension of matrix A
				IPIV,				// integer array for pivoting
				o_x,				// output array
				N,					// leading dimension of array o_x
				info
			);
#if 0
		bool bvalue = true;
		int one = 1;
		zlapmt_(
				bvalue,	// true = forward permutation
				i_size,	// rows
				one,		// cols
				o_x,
				one,		// leading dimension
				IPIV
			);
#endif

#elif 1
		int mid = i_num_diagonals >> 1;

		// zero everything
		for (std::size_t j = 0; j < i_size; j++)
			for (int i = 0; i < i_num_diagonals; i++)
				o_out[j*i_num_diagonals+i] = std::complex<double>(0);

		for (std::size_t idx = 0; idx < i_size; idx++)
		{
			std::size_t i = idx*i_num_diagonals+mid;
			o_out[i] = 1.0/i_in[i];
		}

#else

#endif

	}

};



#endif /* SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_ */
