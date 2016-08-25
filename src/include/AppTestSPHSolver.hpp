/*
 * AppTestSPHSolver.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTSPHSOLVERS_HPP_
#define SRC_TESTSPHSOLVERS_HPP_

#include <SimVars.hpp>
#include <sph/SPHData.hpp>
#include <sph/SPHDataComplex.hpp>
#include <sph/SPHOperators.hpp>
#include <sph/SPHConfig.hpp>
#include <sph/SPHTestSolutions_Gaussian.hpp>
#include <sph/SPHTestSolutions_SPH.hpp>



class AppTestSPHSolvers
{
public:
	SimVars simVars;
	SPHOperators op;


	void run(
			SPHConfig *i_sphConfig
	)
	{
		op.setup(i_sphConfig);

		double epsilon = 1e-10;

		if (true)
		{
			SPHTestSolutions_SPH testSolutionsSph(2,1);

			int tn = 1;
			int tm = 0;
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutionsSph.test_function__grid_gaussian(a,b,c);}
			);
			h.request_data_spectral();


			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(
					[&](double a, double b, double &c){testSolutionsSph.correct_result_diff_mu__grid_gaussian(a,b,c);}
			);
			result.request_data_spectral();

			h.spat_write_file("O_SPHbasis_test_function.csv");


			// one_minus_mu_squared_diff_lat
			h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
			h.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_sph_result.csv");

			result.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_correct_result.csv");


			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;
		}


		{
			SPHTestSolutions_Gaussian testSolutions;


			if (true)
			{
				SPHData data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, double &io_data)
						{
							io_data = y;
						}
				);

				data = data + 100.0;
				data.spec_truncate();

				SPHData data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, double &io_data)
						{
							io_data = y+100.0;
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator+(double): " << error << std::endl;
			}

			if (true)
			{
				cplx offset(200.0, 100.0);

				SPHDataComplex data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, cplx &io_data)
						{
							io_data.real(x*y);
							io_data.imag(x+y);
						}
				);

				data = data+offset;
				data.spec_truncate();


				SPHDataComplex data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, cplx &io_data)
						{
							io_data.real(x*y+offset.real());
							io_data.imag(x+y+offset.imag());
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator+(cplx): " << error << std::endl;
			}



			if (true)
			{
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h.spat_write_file("O_test_function.csv");

				SPHData hphi(i_sphConfig);
				hphi.spat_update_lambda(
						[&](double a, double b, double &c){testSolutions.test_function_phi__grid_phi(a,b,c);}
				);
				hphi.spat_write_file("O_test_function_phi.csv");

				double error_max = h.spat_reduce_error_max(hphi);
				std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;
			}

			if (true)
			{
				// identity
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				//h = h.spat_truncate();

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();

				double error_max = h.spat_reduce_error_max(result);

				std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d lambda
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				//h = h.spat_truncate();

				h = op.diff_lon(h);

				h.spat_write_file("O_diff_lambda_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_diff_lambda__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d phi
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.diff_lat_phi(h);
				h.spat_write_file("O_diff_phi_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_diff_phi__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_phi_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// d/d mu
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.diff_lat_mu(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_diff_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_diff_mu__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_diff_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// grad lambda
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.grad_lon(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_grad_lambda_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_grad_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// mu*F(\lambda,\mu)
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.mu(h);
				h.spat_write_file("O_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_mu__grid_gaussian(a,b,c);}
				);
				result.spat_write_file("O_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST mu*() - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// one_minus_mu_squared_diff_lat
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
				h.spat_write_file("O_one_minus_mu_squared_diff_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian(a,b,c);}
				);

				result.spat_write_file("O_one_minus_mu_squared_diff_mu_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// grad mu
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.grad_lat(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_grad_phi_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_grad_phi__grid_gaussian(a,b,c);}
				);

				result.spat_write_file("O_grad_phi_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// div lambda
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.div_lon(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_div_lambda_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_div_lambda__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();
				result.spat_write_file("O_div_lambda_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// div mu
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.div_lat(h);
				//h = h.spat_truncate();
				h.spat_write_file("O_div_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_div_mu__grid_gaussian(a,b,c);}
				);
				//result = result.spat_truncate();

				result.spat_write_file("O_div_mu_correct_result.csv");

				(h-result).spat_write_file("O_div_mu_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// Laplace
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);

				h = op.div_lon(op.grad_lon(h))
					+ op.div_lat(op.grad_lat(h));

				h.spat_write_file("O_laplace_div_grad_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				result = op.laplace(result);

				result.spat_write_file("O_laplace_laplace_result.csv");

				(h-result).spat_write_file("O_laplace_laplace_z_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST LAPLACE  - max error: " << error_max << std::endl;
			}
		}
	}
};


#endif /* SRC_TESTSPHSOLVERS_HPP_ */
