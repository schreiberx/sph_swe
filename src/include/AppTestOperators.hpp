/*
 * test_operators.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTOPERATORS_HPP_
#define SRC_TESTOPERATORS_HPP_

#include <SimVars.hpp>
#include <sph/SPHData.hpp>
#include <sph/SPHDataComplex.hpp>
#include <sph/SPHOperators.hpp>
#include <sph/SPHConfig.hpp>



class AppTestOperators
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



#if 1
		{
			omp_set_num_threads(1);

			int tn = 1;
			int tm = 0;
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(
					[&](double lon, double mu, double &io_data)
					{
							/**
							 * List of normalized SPH basis functions
							 */

							double m;

							// n=0, m=0  OK
							if (tn == 1 && tm == 0)
							{
								io_data = sqrt(2.0)/2.0;
							}

							// n=1, m=0  OK
							if (tn == 1 && tm == 0)
							{
								io_data = sqrt(6.0)*mu/2.0;
							}

							// n=2, m=0  OK
							if (tn == 2 && tm == 0)
							{
								io_data = sqrt(10.0)/4.0*(3.0*mu*mu-1.0);
							}

							// n=2, m=1
							if (tn == 2 && tm == 1)
							{
								io_data = -mu/2.0*sqrt(-15.0*mu*mu + 15.0);
							}

							// n=2, m=2
							if (tn == 2 && tm == 2)
							{
								io_data = sqrt(15.0)/4.0*(-mu*mu + 1.0);
							}

							// n=3, m=0  OK
							if (tn == 3 && tm == 0)
							{
								io_data = sqrt(14.0)/4.0*mu*(5*mu*mu-3.0);
							}

							// n=3, m=1
							if (tn == 3 && tm == 1)
							{
								io_data = sqrt(-42.0*mu*mu+42.0)*(-5.0*mu*mu/8.0 + 1.0/8.0);
							}

							// n=3, m=2  OK
							if (tn == 3 && tm == 2)
							{
								io_data = sqrt(105.0)*mu/4.0*(-mu*mu+1.0);
							}

							// n=3, m=3  OK
							if (tn == 3 && tm == 3)
							{
								io_data = -sqrt(70.0)/8.0*pow(-mu*mu+1.0, 3.0/2.0);
							}

							io_data *= std::cos(lon*tm);

							io_data *= 1.0/std::sqrt(2.0*M_PI);
					}
			);
			h.request_data_spectral();


			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(
					[&](double lon, double mu, double &io_data)
					{
						/*
						 * List of results
						 */
						double m;

						// n=0, m=0  OK
						if (tn == 1 && tm == 0)
						{
							io_data = 0;
						}

						// n=1, m=0  OK
						if (tn == 1 && tm == 0)
						{
							io_data = sqrt(6.0)/2.0;
						}

						// n=2, m=0  OK
						if (tn == 2 && tm == 0)
						{
							io_data = sqrt(10.0)/4.0*6.0*mu;
						}

						// n=2, m=1
						if (tn == 2 && tm == 1)
						{
//							io_data = -mu/2.0*sqrt(-15.0*mu*mu + 15.0);
							io_data = (15.0*mu*mu-7.5)/(sqrt(15.0)*sqrt(1.0-mu*mu));
						}


#if 0						// n=2, m=2
						if (tn == 2 && tm == 2)
						{
							io_data = sqrt(15.0)/4.0*(-mu*mu + 1.0);
						}

						// n=3, m=0  OK
						if (tn == 3 && tm == 0)
						{
							io_data = sqrt(14.0)/4.0*mu*(5*mu*mu-3.0);
						}

						// n=3, m=1
						if (tn == 3 && tm == 1)
						{
							io_data = sqrt(-42.0*mu*mu+42.0)*(-5.0*mu*mu/8.0 + 1.0/8.0);
						}

						// n=3, m=2  OK
						if (tn == 3 && tm == 2)
						{
							io_data = sqrt(105.0)*mu/4.0*(-mu*mu+1.0);
						}

						// n=3, m=3  OK
						if (tn == 3 && tm == 3)
						{
							io_data = -sqrt(70.0)/8.0*pow(-mu*mu+1.0, 3.0/2.0);
						}
#endif
						io_data *= std::cos(lon*tm);

						io_data *= (1.0-mu*mu);

						io_data *= 1.0/std::sqrt(2.0*M_PI);
					}
			);
			result.request_data_spectral();


#if 0
			std::cout << "TEST FUNCTION" << std::endl;
			h.spec_update_lambda(
					[&](int n, int m, cplx &io_data)
					{
						double r = io_data.real();
						double i = io_data.imag();
						if (std::abs(r) < 1e-10)
							r = 0;
						if (std::abs(i) < 1e-10)
							i = 0;

						std::cout << "n=" << n << ", m=" << m << ": (" << r << ", " << i << ")" << std::endl;
					}
			);
#endif
			h.spat_write_file("O_SPHbasis_test_function.csv");


			// one_minus_mu_squared_diff_lat
			h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
			h.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_sph_result.csv");

			result.spat_write_file("O_SPHbasis_one_minus_mu_squared_diff_mu_correct_result.csv");

#if 0
			std::cout << "SPH RESULT" << std::endl;
			h.spec_update_lambda(
					[&](int n, int m, cplx &io_data)
					{
						double r = io_data.real();
						double i = io_data.imag();
						if (std::abs(r) < 1e-10)
							r = 0;
						if (std::abs(i) < 1e-10)
							i = 0;
						std::cout << "n=" << n << ", m=" << m << ": (" << r << ", " << i << ")" <<  std::endl;
					}
			);
#endif

#if 0
			std::cout << "ANALYTICAL RESULT" << std::endl;
			result.spec_update_lambda(
					[&](int n, int m, cplx &io_data)
					{
						double r = io_data.real();
						double i = io_data.imag();
						if (std::abs(r) < 1e-10)
							r = 0;
						if (std::abs(i) < 1e-10)
							i = 0;
						std::cout << "n=" << n << ", m=" << m << ": (" << r << ", " << i << ")" << std::endl;
					}
			);
#endif

			if (false)
			{
				// mu*F(\lambda,\mu)
				SPHData h(i_sphConfig);
				h.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data)
						{
							//io_data = sqrt(2.0)/2.0;
							//io_data = sqrt(6.0)*mu/2.0;
							io_data = sqrt(10.0)/4.0*(3.0*mu*mu-1.0);
						}
				);
				h = op.mu(h);
				//h.spat_write_file("O_mu_sph_result.csv");

				SPHData result(i_sphConfig);
				result.spat_update_lambda_gaussian_grid(
						[&](double lat, double mu, double &io_data)
						{
							//io_data = sqrt(2.0)/2.0*mu;
							//io_data = sqrt(6.0)*mu/2.0*mu;
							io_data = sqrt(10.0)/4.0*(3.0*mu*mu-1.0)*mu;
						}
				);
				//result.spat_write_file("O_mu_correct_result.csv");

#if 0
				if (false)
				{
					std::cout << "mu P analytical result" << std::endl;
					result.spec_update_lambda(
							[&](int n, int m, cplx &io_data)
							{
								double r = io_data.real();
								double i = io_data.imag();
								if (std::abs(r) < 1e-10)
									r = 0;
								if (std::abs(i) < 1e-10)
									i = 0;
								std::cout << "n=" << n << ", m=" << m << ": (" << r << ", " << i << ")" << std::endl;
							}
					);
				}
#endif

				double error_max = h.spat_reduce_error_max(result);
				std::cout << "TEST mu*() - max error: " << error_max << std::endl;
			}


			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST SPHbasis (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;

			//exit(1);
		}
#endif


#if 0
		{
			SPHData data(i_sphConfig);
			data.spat_update_lambda(
					[&](double x, double y, double &io_data)
					{
						io_data = std::sqrt(2.0)/2.0/sqrt(2.0*M_PI);
					}
			);
			data.request_data_spectral();
			std::cout << data.data_spec[0] << std::endl;
			exit(1);
		}
#endif


		{
			omp_set_num_threads(1);

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

			omp_set_num_threads(0);
		}

		{
			omp_set_num_threads(1);

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

	#if 1
		double exp_fac = 10.0;
	#if 0
		double center_lon = M_PI/3;
		double center_lat = M_PI/2;
	#else
		double center_lon = M_PI/5;
		double center_lat = M_PI/4;
	#endif


		// EXPONENTIAL
		auto test_function__grid_gaussian = [&](double lon, double mu, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp(-d*d*exp_fac);
		};

		auto test_function_phi__grid_phi = [&](double lon, double phi1, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp(-d*d*exp_fac);
		};
		auto correct_result_diff_lambda__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			// http://www.wolframalpha.com/input/?i=diff(exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+c)
			o_data = -
					(exp_fac*2.0*cos(phi1)*cos(phi2)*sin(lambda1-lambda2))
					*d
					*exp(-d*d*exp_fac)
					/sqrt(1-pow(cos(phi1)*cos(phi2)*cos(lambda1-lambda2) + sin(phi1)*sin(phi2), 2.0))
				;
		};
		auto correct_result_diff_phi__grid_gaussian = [&](double lon, double mu, double &o_data)
		{
			// OK
			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			// http://www.wolframalpha.com/input/?i=diff(exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+a)
			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp_fac*2.0
					*d
					*(cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda1-lambda2))
					*exp(-d*d*exp_fac)
					/sqrt(1.0-pow(cos(phi1)*cos(phi2)*cos(lambda1-lambda2) + sin(phi1)*sin(phi2), 2.0))
				;
		};

		auto correct_result_diff_mu__grid_gaussian =
				[&](double lon, double mu, double &o_data)
		{
			correct_result_diff_phi__grid_gaussian(lon, mu, o_data);
			o_data /= sqrt(1.0-mu*mu);
		};

		auto correct_result_grad_lambda__grid_gaussian =
				[&](double lon, double mu, double &o_data)
		{
			correct_result_diff_lambda__grid_gaussian(lon, mu, o_data);
			o_data /= sqrt(1.0-mu*mu);
		};

		auto correct_result_mu__grid_gaussian =
				[&](double lon, double mu, double &o_data)
		{
			test_function__grid_gaussian(lon, mu, o_data);
			o_data *= mu;
		};
		auto correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian =
				[&](double lon, double mu, double &o_data)
		{
			correct_result_diff_mu__grid_gaussian(lon, mu, o_data);
			o_data *= (1.0-mu*mu);
		};
		auto correct_result_grad_phi__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			correct_result_diff_phi__grid_gaussian(lon, mu, o_data);
		};

		auto correct_result_div_lambda__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			correct_result_grad_lambda__grid_gaussian(lon, mu, o_data);
		};

		auto correct_result_div_mu__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			// http://www.wolframalpha.com/input/?i=diff(cos(a)*exp(-pow(arccos(sin(a)*sin(b)+%2B+cos(a)*cos(b)*cos(c-d)),2)*A),+a)
			double da = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2);
			double d = acos(da);

			o_data = exp(-exp_fac*d*d)*
					(
						2.0*exp_fac*cos(phi1)
							*d
							*(cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda1-lambda2))
							/sqrt(1.0-da*da)

						-sin(phi1)
					)
				;
			o_data /= cos(phi1);
		};
	#endif


	#if 0
		// ALP
		auto test_function__grid_gaussian = [&](double lon, double mu, double &o_data)
		{
			// P_0^1
			o_data = sqrt(6.0)*mu*0.5;
		};
		auto correct_result_diff_lambda__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			o_data = 0;
		};
		auto correct_result_grad_lambda__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			o_data = 0;
		};

		auto correct_result_diff_mu__grid_gaussian = [&](double lon, double mu, double &o_data)
		{
			o_data = sqrt(6.0)*0.5;
		};
		auto correct_result_grad_phi__grid_gaussian = [=](double lon, double mu, double &o_data)
		{
			correct_result_diff_mu__grid_gaussian(lon, mu, o_data);
			o_data *= sqrt(1.0-mu*mu);
		};

	#endif

	#if 0
		// f = cos(lambda)*cos(phi)
		auto test_function__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = cos(lon)*cos_phi;
		};
		auto correct_result_diff_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = -sin(lon)*cos_phi;
		};
		auto correct_result_diff_mu__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = -cos(lon)*mu/cos_phi;
		};
		auto correct_result_grad_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = -sin(lon)*cos_phi;
			o_data /= cos_phi;
		};
	#endif


	#if 0
		// f = sin(phi)
		auto test_function__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = sin_phi;
		};

		auto correct_result_diff_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = 0;
		};
		auto correct_result_grad_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = 0;
			o_data /= cos_phi;
		};
	#endif


	#if 0
		// f = sin(lambda)
		auto test_function__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = sin_lambda;
		};

		auto correct_result_diff_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = sin_lambda;
		};

		auto correct_result_grad_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			double sin_phi = mu;
			double cos_phi = sqrt(1.0-mu*mu);
			double sin_lambda = sin(lon);
			double cos_lambda = cos(lon);

			o_data = sin_lambda;
			o_data /= cos_phi;
		};
	#endif

	#if 0
		// f = sin(lambda)
		auto test_function__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			o_data = (double)random()/(double)RAND_MAX;
		};

		auto correct_result_diff_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			o_data = 0;
		};

		auto correct_result_grad_lambda__grid_gaussian = [](double lon, double mu, double &o_data)
		{
			o_data = 0;
		};
	#endif

#if 1
		{
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h.spat_write_file("O_test_function.csv");

			SPHData hphi(i_sphConfig);
			hphi.spat_update_lambda(test_function_phi__grid_phi);
			hphi.spat_write_file("O_test_function_phi.csv");

			double error_max = h.spat_reduce_error_max(hphi);
			std::cout << "TEST PHI vs. MU: max error: " << error_max << std::endl;
		}

		if (true)
		{
			// identity
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			//h = h.spat_truncate();

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			//result = result.spat_truncate();

			double error_max = h.spat_reduce_error_max(result);

			std::cout << "TEST IDENTITY: max error: " << error_max << std::endl;
		}


		if (true)
		{
			// d/d lambda
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			//h = h.spat_truncate();

			h = op.diff_lon(h);

			h.spat_write_file("O_diff_lambda_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_diff_lambda__grid_gaussian);
			//result = result.spat_truncate();
			result.spat_write_file("O_diff_lambda_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST DIFF LON - max error: " << error_max << std::endl;
		}


		if (true)
		{
			// d/d phi
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.diff_lat_phi(h);
			h.spat_write_file("O_diff_phi_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_diff_phi__grid_gaussian);
			//result = result.spat_truncate();
			result.spat_write_file("O_diff_phi_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST DIFF PHI - max error: " << error_max << std::endl;
		}

		if (true)
		{
			// d/d mu
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.diff_lat_mu(h);
			//h = h.spat_truncate();
			h.spat_write_file("O_diff_mu_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_diff_mu__grid_gaussian);
			//result = result.spat_truncate();
			result.spat_write_file("O_diff_mu_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST DIFF LAT MU - max error: " << error_max << std::endl;
		}


		if (true)
		{
			// grad lambda
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.grad_lon(h);
			//h = h.spat_truncate();
			h.spat_write_file("O_grad_lambda_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_grad_lambda__grid_gaussian);
			//result = result.spat_truncate();
			result.spat_write_file("O_grad_lambda_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST GRAD LON - max error: " << error_max << std::endl;
		}


		if (true)
		{
			// mu*F(\lambda,\mu)
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.mu(h);
			h.spat_write_file("O_mu_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_mu__grid_gaussian);
			result.spat_write_file("O_mu_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST mu*() - max error: " << error_max << std::endl;
		}

		if (true)
		{
			// one_minus_mu_squared_diff_lat
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.spec_one_minus_mu_squared_diff_lat_mu(h);
			h.spat_write_file("O_one_minus_mu_squared_diff_mu_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_one_minus_mu_squared_diff_lat_mu__grid_gaussian);

			result.spat_write_file("O_one_minus_mu_squared_diff_mu_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST (1-mu*mu)*d/dmu - max error: " << error_max << std::endl;
		}
		if (true)
		{
			// grad mu
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.grad_lat(h);
			//h = h.spat_truncate();
			h.spat_write_file("O_grad_phi_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_grad_phi__grid_gaussian);

			result.spat_write_file("O_grad_phi_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST GRAD LAT - max error: " << error_max << std::endl;
//exit(1);
		}


		if (true)
		{
			// div lambda
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.div_lon(h);
			//h = h.spat_truncate();
			h.spat_write_file("O_div_lambda_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_div_lambda__grid_gaussian);
			//result = result.spat_truncate();
			result.spat_write_file("O_div_lambda_correct_result.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST DIV LON  - max error: " << error_max << std::endl;
		}


		if (true)
		{
			// div mu
			SPHData h(i_sphConfig);
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			h = op.div_lat(h);
			//h = h.spat_truncate();
			h.spat_write_file("O_div_mu_sph_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(correct_result_div_mu__grid_gaussian);
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
			h.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);

			h = op.div_lon(op.grad_lon(h))
				+ op.div_lat(op.grad_lat(h));

			//h = h.spat_truncate();

			h.spat_write_file("O_laplace_div_grad_result.csv");

			SPHData result(i_sphConfig);
			result.spat_update_lambda_gaussian_grid(test_function__grid_gaussian);
			result = op.laplace(result);
			//result = result.spat_truncate();

			result.spat_write_file("O_laplace_laplace_result.csv");

			(h-result).spat_write_file("O_laplace_laplace_z_diff.csv");

			double error_max = h.spat_reduce_error_max(result);
			std::cout << "TEST LAPLACE  - max error: " << error_max << std::endl;
		}
#endif
	}
};


#endif /* SRC_TESTOPERATORS_HPP_ */
