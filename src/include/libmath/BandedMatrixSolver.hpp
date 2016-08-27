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
 * LAPACK description:
 *
 * For an array
 * 	REAL A( LDA, * )
 *
 * LDA is called the leading dimension of the array.
 * In Fortran, the values in this dimension
 * are consecutively stored in memory.
 *
 * aij = A(i,j)
 *
 *    a11 a12 a13 a14   |
 *    a21 a22 a23 a24   |
 *    a31 a32 a33 a34   |LDA
 *    a41 a42 a43 a44   v
 *
 */
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
	int max_N;
	int num_diagonals;
	int num_halo_size_diagonals;

	int LDAB;

	std::complex<double>* AB;
	int *IPIV;

	MatrixSolversCommon()	:
		AB(nullptr),
		IPIV(nullptr)
	{
	}



	~MatrixSolversCommon()
	{
		shutdown();
	}

	void setup(
			int i_max_N,			///< size of the matrix
			int i_num_off_diagonals		///< number of block diagonals
	)
	{
		max_N = i_max_N;
		num_diagonals = 2*i_num_off_diagonals+1;
		num_halo_size_diagonals = i_num_off_diagonals;

		assert(2*num_halo_size_diagonals+1 == num_diagonals);

		LDAB = 2*num_halo_size_diagonals + num_halo_size_diagonals + 1;

		AB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*LDAB*i_max_N);
		IPIV = (int*)malloc(sizeof(int)*i_max_N);
	}

	void shutdown()
	{
		if (IPIV != nullptr)
		{
			free(IPIV);
			IPIV = nullptr;
		}

		if (AB != nullptr)
		{
			free(AB);
			AB = nullptr;
		}
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
	void solve_diagBandedInverse_C(
		const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
		const std::complex<double>* i_b,		///< RHS of equation and output of solution X
		std::complex<double>* o_x,
		int i_size
	);


	void solve_diagBandedInverse_Fortran(
		const std::complex<double>* i_A,		///< Matrix for input and in-place transformations
		const std::complex<double>* i_b,		///< RHS of equation and output of solution X
		std::complex<double>* o_x,
		int i_size
	);
};


template <>
class BandedMatrixSolver<std::complex<double>>	:
		public MatrixSolversCommon<std::complex<double>>
{
	/**
	 * Solve for input matrix
	 *
	 * i_A: width: num_diagonals
	 *      height: i_size
	 *
	 * i_b: RHS of equation
	 *
	 * o_x: Solution
	 */
public:
	void solve_diagBandedInverse_C(
		const std::complex<double>* i_A,
		const std::complex<double>* i_b,
		std::complex<double>* o_x,
		int i_size
	)
	{
		assert(max_N >= i_size);

		// zero everything
		for (std::size_t j = 0; j < i_size; j++)
			for (int i = 0; i < num_diagonals; i++)
				o_x[j*num_diagonals+i] = std::complex<double>(0);

#if 0
		for (std::size_t idx = 0; idx < i_size; idx++)
		{
			o_x[idx] = i_b[idx]/i_A[idx*num_diagonals+num_halo_size_diagonals];
		}
#else
		/*
		 * Convert to Fortran storage array
		 */
#ifndef NDEBUG
		//std::cout << "num_diagonals: " << num_diagonals << std::endl;
		for (int i = 0; i < LDAB*i_size; i++)
			AB[i] = 666.0;
#endif
		assert(max_N >= i_size);

		// c columns / fortran rows
		for (int i = 0; i < num_diagonals; i++)
		{
			// c rows / fortran columns
			for (int j = 0; j < i_size; j++)
			{
				// AB is LDAB large!
				assert(LDAB*max_N > i*i_size+j);

				// note, that the array is upside-down
				//AB[j*LDAB + (LDAB-i-1)] = i_A[j*num_diagonals + i];
				AB[j*LDAB + (LDAB-(num_diagonals-i-1)-1)] = i_A[j*num_diagonals + i];
			}
		}

#if 1
		std::cout << "*** A MATRIX IN FORTRAN FORMAT ***" << std::endl;
		for (int i = 0; i < LDAB; i++)
		{
			for (int j = 0; j < i_size; j++)
			{
				std::cout << AB[j*LDAB + i] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
#endif
		solve_diagBandedInverse_Fortran_largeA(AB, i_b, o_x, i_size);

#endif
	}


public:
	void solve_diagBandedInverse_Fortran_largeA(
		const std::complex<double>* i_A,	///< A of max size
		const std::complex<double>* i_b,
		std::complex<double>* o_x,
		int i_size
	)
	{
		assert(num_diagonals & 1 == 1);
		assert(AB != nullptr);

		/*
		 * Make a copy of the array data since this is a destructive function
		 */
		if (AB != i_A)
			memcpy((void*)AB, (const void*)i_A, sizeof(std::complex<double>)*num_diagonals*LDAB);

		memcpy((void*)o_x, (const void*)i_b, sizeof(std::complex<double>)*i_size);

		int one = 1;
		int info;

#if 0
		std::cout << "************************************" << std::endl;
		std::cout << "i_size: " << i_size << std::endl;
		std::cout << "num_halo_size_diagonals: " << num_halo_size_diagonals << std::endl;
		std::cout << "LDAB: " << LDAB << std::endl;
#endif

		zgbsv_(
				i_size,				// number of linear equations
				num_halo_size_diagonals,	// number of subdiagonals
				num_halo_size_diagonals,	// number of superdiagonals
				one,				// number of columns of matrix B
				AB,					// array with matrix A to solve for
				LDAB,				// leading dimension of matrix A
				IPIV,				// integer array for pivoting
				o_x,				// output array
				i_size,				// leading dimension of array o_x
				info
			);

		if (info != 0)
		{
			std::cerr << "zgbsv returned INFO != 0: " << info << std::endl;
			assert(false);
			exit(1);
		}

#if 1
		bool bvalue = false;
		zlapmt_(
				bvalue,	// true = forward permutation
				i_size,	// rows
				one,		// cols
				o_x,
				i_size,		// leading dimension
				IPIV
			);
#endif
	}

};



#endif /* SRC_INCLUDE_LIBMATH_BANDEDMATRIXSOLVER_HPP_ */
