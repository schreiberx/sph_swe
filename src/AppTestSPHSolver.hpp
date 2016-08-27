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

			if (sphConfig->spec_n_max < 32)
			{
				std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
			}

			/**
			 * Use test function as expected result
			 */
			SPHData x_result(i_sphConfig);
			x_result.spat_update_lambda_gaussian_grid(
					[&](double lat, double mu, double &io_data){
						testSolutions.test_function__grid_gaussian(lat,mu,io_data);
					}
			);

			std::complex<double> scalar_a = 3.0;

			/*
			 * Test Z1 = c*Phi(mu)
			 */
			if (true)
			{
				std::cout << "Test Z1 = c*Phi(mu)";

				SPHSolver<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				sphSolver.solver_component_scalar_phi(scalar_a);

				/*
				 * Setup RHS = scalar_a * phi(lambda,mu)
				 */
				SPHData b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data){
							testSolutions.test_function__grid_gaussian(lat,mu,io_data);
							io_data *= scalar_a.real();
						}
				);

				SPHData x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
			}

			/*
			 * Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)";

				SPHSolver<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				sphSolver.solver_component_scalar_phi(scalar_a);
				sphSolver.solver_component_mu_phi();
//				sphSolver.lhs.print();

				SPHData b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data){
							double phi;
							testSolutions.test_function__grid_gaussian(lat,mu,phi);

							io_data = mu*phi+scalar_a.real()*phi;
						}
				);

				SPHData x_numerical = sphSolver.solve(b);

//				x_numerical.spat_write_file("O_numerical.csv");
//				x_result.spat_write_file("O_result.csv");
//				(x_result-x_numerical).spat_write_file("O_diff.csv");

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
			}


			/*
			 * Test Zx = mu*mu*Phi(lam,mu) + a*Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = mu*mu*Phi(lam,mu) + a*Phi(lam,mu)";

				SPHSolver<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				sphSolver.solver_component_scalar_phi(scalar_a);
				sphSolver.solver_component_mu_mu_phi();
//				sphSolver.lhs.print();

				SPHData b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data){
							testSolutions.test_function__grid_gaussian(lat,mu,io_data);
							io_data = (mu*mu+scalar_a.real())*io_data;
						}
				);

				SPHData x_numerical = sphSolver.solve(b);

//				x_numerical.spat_write_file("O_numerical.csv");
//				x_result.spat_write_file("O_result.csv");

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
			}


			/*
			 * Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)";

				SPHSolver<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 1);

				sphSolver.solver_component_scalar_phi(scalar_a);
				sphSolver.solver_component_one_minus_mu_mu_diff_mu_phi();
//				sphSolver.lhs.print();

				SPHData b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data)
						{
							double phi;
							testSolutions.test_function__grid_gaussian(lat,mu,phi);

							double om_dphi;
							testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(lat,mu,om_dphi);

							io_data = om_dphi + scalar_a.real()*phi;
						}
				);

				SPHData x_numerical = sphSolver.solve(b);


				x_numerical.spat_write_file("O_numerical.csv");
				x_result.spat_write_file("O_result.csv");
				(x_result-x_numerical).spat_write_file("O_diff.csv");

				std::cout << std::endl;
				std::cout << "SPECTRAL DIFFERENCE" << std::endl;
				(x_result-x_numerical).spec_print();
				std::cout << std::endl;

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
			}
		}
	}
};


#endif /* SRC_TESTSPHSOLVERS_HPP_ */
