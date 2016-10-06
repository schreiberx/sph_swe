/*
 * AppTestSPHOperators.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTOPERATORS_HPP_
#define SRC_TESTOPERATORS_HPP_

#include <benchmarks/SphereTestSolutions_Gaussian.hpp>
#include <benchmarks/SphereTestSolutions_SPH.hpp>
#include <SimVars.hpp>
#include <sph/SPHData.hpp>
#include <sph/SPHOperators.hpp>
#include <sph/SPHConfig.hpp>



class AppTestSPHOperators
{
public:
	SimVars simVars;
	SPHOperators op;


	void run(
			SPHConfig *i_sphConfig
	)
	{
		double epsilon = 1e-10;

		if (true)
		{
			SphereTestSolutions_SPH testSolutionsSph(2,1);

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
			SphereTestSolutions_Gaussian testSolutions;


			if (true)
			{
				SPHData data(i_sphConfig);
				data.spat_update_lambda(
						[&](double x, double y, double &io_data)
						{
							io_data = y;
						}
				);

				data = data*123.0;
				data.spec_truncate();

				SPHData data2(i_sphConfig);
				data2.spat_update_lambda(
						[&](double x, double y, double &io_data)
						{
							io_data = y*123.0;
						}
				);
				data2.spec_truncate();

				double error = data.spat_reduce_error_max(data2);
				std::cout << "ERROR operator*123.0: " << error << std::endl;
			}

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
				// mu*mu*F(\lambda,\mu)
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.mu2(h);
				h.spat_write_file("O_mu2_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &i_data){
							testSolutions.test_function__grid_gaussian(lat, mu, i_data);
							i_data *= mu*mu;
						}
				);
				result.spat_write_file("O_mu2_correct_result.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST mu*mu*() - max error: " << error_max << std::endl;
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
				h.spat_write_file("O_div_lambda_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_div_lambda__grid_gaussian(a,b,c);}
				);
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
				h.spat_write_file("O_div_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.correct_result_div_mu__grid_gaussian(a,b,c);}
				);

				result.spat_write_file("O_div_mu_correct_result.csv");
				(h-result).spat_write_file("O_div_mu_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV LAT  - max error: " << error_max << std::endl;
			}
#if 0

			if (true)
			{
				// div mu TEST
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.div_lat_TEST(h);
				h.spat_write_file("O_div_mu_TEST_sph_result.csv");
#if 0
				h.spat_update_lambda_cogaussian_grid(
						[&](double a, double mu, double &c){
							c *= mu;
						}
				);
#endif
				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double mu, double &c){
							testSolutions.correct_result_div_mu__grid_gaussian(a,mu,c);
						}
				);
				result.spat_update_lambda_cogaussian_grid(
						[&](double a, double mu, double &c){
							c *= mu;
						}
				);

				result.spat_write_file("O_div_mu_TEST_correct_result.csv");
				(h-result).spat_write_file("O_div_mu_TEST_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIV TEST LAT  - max error: " << error_max << std::endl;
			}
#endif

			if (true)
			{
				// divergence
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.div(h, h);
//				h.spat_update_lambda_gaussian_grid([&](double a, double mu, double &c){c *= mu*mu*mu*mu;});
				h.spat_write_file("O_divergence_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){

							double data1;
							testSolutions.correct_result_div_mu__grid_gaussian(a,b,data1);

							double data2;
							testSolutions.correct_result_div_lambda__grid_gaussian(a,b,data2);

							c = data1 + data2;
						}
				);
//				result.spat_update_lambda_gaussian_grid([&](double a, double mu, double &c){c *= mu*mu*mu*mu;});
				//result = result.spat_truncate();

				result.spat_write_file("O_divergence_correct_result.csv");

				(h-result).spat_write_file("O_divergence_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST DIVERGENCE - max error: " << error_max << std::endl;
			}


			if (true)
			{
				// vorticity
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){testSolutions.test_function__grid_gaussian(a,b,c);}
				);
				h = op.vort(h, h);
				//h = h.spat_truncate();
				h.spat_write_file("O_vort_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double a, double b, double &c){

							double data1;
							testSolutions.correct_result_grad_lambda__grid_gaussian(a,b,data1);

							double data2;
							testSolutions.correct_result_grad_phi__grid_gaussian(a,b,data2);

							c = data1 - data2;
						}
				);
				//result = result.spat_truncate();

				result.spat_write_file("O_vort_correct_result.csv");

				(h-result).spat_write_file("O_vort_correct_diff.csv");

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST VORTICITY - max error: " << error_max << std::endl;
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
				std::cout << "TEST LAPLACE (div.grad vs. sph laplace) - max error: " << error_max << std::endl;
			}
		}
	}
};


#endif /* SRC_TESTOPERATORS_HPP_ */
