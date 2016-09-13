/*
 * AppTestSPHSolverComplex.hpp
 *
 *  Created on: 31 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTSPHSOLVERS_COMPLEX_HPP_
#define SRC_TESTSPHSOLVERS_COMPLEX_HPP_

#include <benchmarks/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks/SphereTestSolutions_SPH.hpp>
#include <SimVars.hpp>
#include <sph/SPHDataComplex.hpp>
#include <sph/SPHOperatorsComplex.hpp>
#include <sph/SPHConfig.hpp>
#include <sph/SPHSolverComplex.hpp>
#include <rexi/REXI.hpp>


class AppTestSPHSolversComplex
{
public:
	SimVars simVars;
	SPHOperatorsComplex opComplex;

	SPHConfig *sphConfig;

	double threshold = 1e-10;

	AppTestSPHSolversComplex()
	{
	}



	void test_threshold(double i_value)
	{
		if (std::abs(i_value) > threshold)
		{
			std::cout << "!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!! threshold exceeded !!!!!!!!!!!!!!!!" << std::endl;
			exit(1);
		}
	}



	/**
	 * Run with
	 *
	 * 	$ ./build/sh_example T32 P2
	 */
	void run(
			SPHConfig *i_sphConfig
	)
	{
		std::cout << std::setprecision(10);

		sphConfig = i_sphConfig;

		SphereTestSolutions_Gaussian testSolutions;

		if (sphConfig->spec_n_max < 32)
		{
			std::cerr << "WARNING: AT LEAST 32 MODES REQUIRED for proper accuracy!!!" << std::endl;
		}

		/**
		 * Use test function as expected result
		 */
		SPHDataComplex x_result(i_sphConfig);
		x_result.spat_update_lambda_gaussian_grid(
				[&](double lat, double mu, std::complex<double> &io_data){
					double tmp;
					testSolutions.test_function__grid_gaussian(lat,mu,tmp);
					io_data = tmp;
				}
		);

		double r = simVars.earth_radius;
		double two_omega = 2.0*simVars.coriolis_omega;


		REXI rexi;
		rexi.setup(0.2, 64);


		for (int i = 0; i < rexi.alpha.size(); i++)
		{
			std::complex<double> &alpha = rexi.alpha[i];
			std::cout << std::endl;
			std::cout << "***************************************************" << std::endl;
			std::cout << "Testing for alpha " << alpha << std::endl;
			std::cout << "***************************************************" << std::endl;


			/*
			 * Test Zx = c*Phi(mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = c*Phi(mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				sphSolver.solver_component_scalar_phi(alpha);

				/*
				 * Setup RHS = scalar_a * phi(lambda,mu)
				 */
				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat,mu,tmp);
							io_data = tmp*alpha;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}

			/*
			 * Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = mu*Phi(lam,mu) + a*Phi(lam,mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				sphSolver.solver_component_scalar_phi(alpha);
				sphSolver.solver_component_mu_phi();
//				sphSolver.lhs.print();

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double fun;
							testSolutions.test_function__grid_gaussian(lat,mu,fun);

							io_data = mu*fun+alpha*fun;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

//				x_numerical.spat_write_file("O_numerical.csv");
//				x_result.spat_write_file("O_result.csv");
//				(x_result-x_numerical).spat_write_file("O_diff.csv");

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}

			/*
			 * Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Zx = (1-mu*mu)*d/dmu Phi(lam,mu) + a*Phi(lam,mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> a(alpha);
				sphSolver.solver_component_scalar_phi(a);
				sphSolver.solver_component_one_minus_mu_mu_diff_mu_phi();

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double fun;
							testSolutions.test_function__grid_gaussian(lat, mu, fun);

							double one_minus_m2_dfun;
							testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(lat, mu, one_minus_m2_dfun);

							io_data = one_minus_m2_dfun + alpha*fun;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}

			std::cout << "************************************************************" << std::endl;
			std::cout << "*** BASIC TESTS FINISHED ***********************************" << std::endl;
			std::cout << "************************************************************" << std::endl;
			std::cout << "*** TESTING REXI COMPONENTS ********************************" << std::endl;
			std::cout << "************************************************************" << std::endl;

			/*
			 * Test F1 = alpha^4 * Phi(mu)
			 */
			if (true)
			{
				std::cout << "Test Z1 = alpha^4*Phi(mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = (alpha*alpha)*(alpha*alpha);
				sphSolver.solver_component_rexi_z1(scalar, r);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat,mu,tmp);
							io_data = tmp*scalar;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}

			/*
			 * Test Z2 = mu^2 * Phi(mu)
			 */
			if (true)
			{
				std::cout << "Test Z2 = mu^2*Phi(lam,mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = alpha*alpha*two_omega*two_omega;
				//sphSolver.solver_component_scalar_phi(scalar_a);
				sphSolver.solver_component_rexi_z2(scalar, r);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data){
							double tmp;
							testSolutions.test_function__grid_gaussian(lat,mu,tmp);
							io_data = scalar*(mu*mu)*tmp;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test Z3 = mu^4 * Phi(mu)
			 */
			if (true)
			{
				std::cout << "Test Z3 = mu^4*Phi(lam,mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 4);

				std::complex<double> scalar = two_omega*two_omega*two_omega*two_omega;
				//sphSolver.solver_component_scalar_phi(scalar_a);
				sphSolver.solver_component_rexi_z3(scalar, r);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat,mu,tmp);
							io_data = scalar*(mu*mu)*(mu*mu)*tmp;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)
			 */
			if (true)
			{
				std::cout << "Test Z4 = grad_j(mu) grad_i(Phi(lam,mu)) = d/dlambda Phi(lam,mu)";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = -alpha*alpha*two_omega;
				sphSolver.solver_component_rexi_z4(scalar, r);

				// ADD OFFSET FOR NON-SINGULAR SOLUTION
				sphSolver.solver_component_scalar_phi(alpha);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double fun;
							testSolutions.test_function__grid_gaussian(lat, mu, fun);

							double dlambda_fun;
							testSolutions.correct_result_diff_lambda__grid_gaussian(lat, mu, dlambda_fun);

							io_data = scalar/(r*r)*dlambda_fun + fun*alpha;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))
			 */
			if (true)
			{
				std::cout << "Test Z5 = grad_j(mu) mu^2 grad_i(Phi(lam,mu))";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = two_omega*two_omega*two_omega;
				sphSolver.solver_component_rexi_z5(scalar, r);

				// ADD OFFSET FOR NON-SINGULAR SOLUTION
				sphSolver.solver_component_scalar_phi(alpha);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double fun;
							testSolutions.test_function__grid_gaussian(lat, mu, fun);

							double dgrad_lambda_fun;
							testSolutions.correct_result_grad_lambda__grid_gaussian(lat, mu, dgrad_lambda_fun);

							io_data = scalar/(r*r)*std::sqrt(1.0-mu*mu)*mu*mu*dgrad_lambda_fun + fun*alpha;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))
			 */
			if (true)
			{
				std::cout << "Test Z6 = grad_j(mu) mu grad_j(Phi(lam,mu))";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = 2.0*alpha*two_omega*two_omega;
				sphSolver.solver_component_rexi_z6(scalar, r);

				// ADD OFFSET FOR NON-SINGULAR SOLUTION
				sphSolver.solver_component_scalar_phi(alpha);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double fun;
							testSolutions.test_function__grid_gaussian(lat, mu, fun);

							double dgrad_mu_fun;
							testSolutions.correct_result_grad_phi__grid_gaussian(lat, mu, dgrad_mu_fun);

							io_data = scalar/(r*r)*std::sqrt(1.0-mu*mu)*mu*dgrad_mu_fun + fun*alpha;
						}
				);

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test Z7 = laplace(Phi(lam,mu))
			 */
			if (true)
			{
				std::cout << "Test Z7 = laplace(Phi(lam,mu))";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = 1.0;
				sphSolver.solver_component_rexi_z7(scalar, r);

				// ADD OFFSET FOR NON-SINGULAR SOLUTION
				sphSolver.solver_component_scalar_phi(alpha);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat, mu, tmp);
							io_data = tmp;
						}
				);

				b = opComplex.laplace(b)*scalar/(r*r) + b*alpha;

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}



			/*
			 * Test Z8 = mu*mu*laplace(Phi(lam,mu))
			 */
			if (true)
			{
				std::cout << "Test Z8 = mu*mu*laplace(Phi(lam,mu))";

				SPHSolverComplex<std::complex<double>> sphSolver;
				sphSolver.setup(sphConfig, 2);

				std::complex<double> scalar = 1.0;
				sphSolver.solver_component_rexi_z8(scalar, r);

				// ADD OFFSET FOR NON-SINGULAR SOLUTION
				sphSolver.solver_component_scalar_phi(alpha);

				SPHDataComplex b(i_sphConfig);
				b.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat, mu, tmp);
							io_data = tmp;
						}
				);

				b = opComplex.mu2(opComplex.laplace(b))*scalar/(r*r) + b*alpha;

				SPHDataComplex x_numerical = sphSolver.solve(b);

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}


			/*
			 * Test direct Helmholtz problem solver
			 * (a + b D^2) x = rhs
			 */
			if (true)
			{
				std::cout << "Test Zx = a + b*laplace";

				SPHDataComplex testb(i_sphConfig);
				testb.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, std::complex<double> &io_data)
						{
							double tmp;
							testSolutions.test_function__grid_gaussian(lat, mu, tmp);
							io_data = tmp;
						}
				);

				std::complex<double> a = alpha;
				std::complex<double> b = alpha*2.0-1.0;
				double r = simVars.earth_radius;

				testb = a*testb + (b/(r*r))*opComplex.laplace(testb);

				SPHDataComplex x_numerical = testb.spec_solve_helmholtz(a, b, r);

				//testb.spat_print();
				//x_numerical.spat_print();

				double error_max = x_numerical.spat_reduce_error_max(x_result);
				std::cout << " ||| Error: " << error_max << std::endl;
				test_threshold(error_max);
			}

		}
	}
};


#endif /* SRC_TESTSPHSOLVERS_HPP_ */
