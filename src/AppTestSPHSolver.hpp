/*
 * AppTestSPHSolver.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTSPHSOLVERS_HPP_
#define SRC_TESTSPHSOLVERS_HPP_

#include <benchmarks/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks/SphereTestSolutions_SPH.hpp>
#include <SimVars.hpp>
#include <sph/SPHData.hpp>
#include <sph/SPHDataComplex.hpp>
#include <sph/SPHOperators.hpp>
#include <sph/SPHConfig.hpp>
#include <sph/SPHSolver.hpp>



class AppTestSPHSolvers
{
public:
	SimVars simVars;
	SPHOperators op;

	SPHConfig *sphConfig;

	/**
	 * Run with
	 *
	 * 	$ ./build/sh_example T32 P2
	 */
	void run(
			SPHConfig *i_sphConfig
	)
	{
		sphConfig = i_sphConfig;

		{
			SphereTestSolutions_Gaussian testSolutions;

			if (true)
			{
				//assert(sphConfig->spec_n_max >= 32);

//				std::cout << sphConfig->spec_n_max << std::endl;
				SPHSolver<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				/*
				 * Solve
				 * a * x = b
				 *
				 * With b the manufactured solution
				 */
				double scalar_a = 2.0;
				sphSolver.solver_component_const(scalar_a);

#if 0
				sphSolver.lhs.print();

				sphSolver.lhs.convertToFortranArray();
				sphSolver.lhs.printFortran();
//				exit(1);
#endif


#if 0
				std::cout << std::endl;
				std::cout << "LHS" << std::endl;
				sphSolver.lhs.print();
				std::cout << std::endl;
#endif

				// d/d phi
				SPHData b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &io_data){
							testSolutions.test_function__grid_gaussian(a,b,io_data);
						}
				);

				SPHData x_numerical = sphSolver.solve(b);

				/*
				 * Compute expected result
				 */
				SPHData x_result(i_sphConfig);
				x_result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &io_data){
							testSolutions.test_function__grid_gaussian(a,b,io_data);
							io_data /= scalar_a;
						}
				);
				//result = result.spat_truncate();
				//x_result.spat_write_file("O_diff_phi_correct_result.csv");

//				x_result.spat_truncate();
//				x_numerical.spat_truncate();

#if 0
				std::cout << "****************************************" << std::endl;
				std::cout << "x_result" << std::endl;
				std::cout << "****************************************" << std::endl;
				x_result.print();
				std::cout << std::endl;

				std::cout << "****************************************" << std::endl;
				std::cout << "x_numerical" << std::endl;
				std::cout << "****************************************" << std::endl;
				x_numerical.print();
				std::cout << std::endl;
#endif
				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << "TEST a*x=b - Nmax(x,x') = " << error_max << std::endl;
			}
		}
	}
};


#endif /* SRC_TESTSPHSOLVERS_HPP_ */
