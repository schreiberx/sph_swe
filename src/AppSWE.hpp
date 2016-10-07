/*
 * AppTestSWE.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTSWE_HPP_
#define SRC_TESTSWE_HPP_

#include <benchmarks/BenchmarkGalewsky.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <rexi/REXI.hpp>
#include <sph/SPHDataComplex.hpp>
#include <sph/SPHSolver.hpp>

#include "SWERexiSPH.hpp"
#include "SWERexiSPHRobert.hpp"


class AppSWE
{
public:
	SimVars &simVars;
	SPHOperators op;
	SPHOperatorsComplex opComplex;

	SPHConfig *sphConfig;

	// Runge-Kutta stuff
	TimesteppingRK timestepping;


	SPHData prog_h;
	SPHData prog_u;
	SPHData prog_v;
	SPHData fdata;


	BenchmarkGalewsky benchmarkGalewsky;

	REXI rexi;



public:
	AppSWE(
			SPHConfig *i_sphConfig,
			SimVars &i_simVars
	)	:
		sphConfig(i_sphConfig),
		prog_h(i_sphConfig),
		prog_u(i_sphConfig),
		prog_v(i_sphConfig),
		fdata(i_sphConfig),
		simVars(i_simVars),
		benchmarkGalewsky(i_simVars)
	{
	}



	void write_output()
	{
		char buffer[1024];

		sprintf(buffer, "prog_h_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.benchmark_scenario_id == 0)
			prog_h.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_h.spat_write_file(buffer);
		std::cout << buffer << " (min: " << prog_h.spat_reduce_min() << ", max: " << prog_h.spat_reduce_max() << ")" << std::endl;

		sprintf(buffer, "prog_u_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.benchmark_scenario_id == 0)
			prog_u.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_u.spat_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_v_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (simVars.benchmark_scenario_id == 0)
			prog_v.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_v.spat_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_eta_t%020.8f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		SPHData vort = op.vort(prog_u, prog_v)/simVars.earth_radius;
		if (simVars.benchmark_scenario_id == 0)
			vort.spat_write_file_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			vort.spat_write_file(buffer);
		std::cout << buffer << std::endl;
	}



	void setup_initial_conditions_gaussian(double i_center_lat = M_PI/3)
	{
		double exp_fac = 10.0;

#if 0
		double center_lon = M_PI/3;
		double center_lat = M_PI/5;
#else
		double center_lon = 0;//M_PI;
		double center_lat = M_PI/3;
		center_lat = i_center_lat;
#endif


		auto initial_condition_h = [&](double lon, double mu, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp(-d*d*exp_fac)*0.1*simVars.h0 + simVars.h0;
		};

		prog_h.spat_update_lambda_gaussian_grid(initial_condition_h);
		prog_u.spat_set_zero();
		prog_v.spat_set_zero();
	}



	SPHData f(SPHData i_sphData)
	{
//		return op.coriolis(i_sphData, simVars.coriolis_omega);
		return fdata*i_sphData*2.0*simVars.coriolis_omega;
	}



	void run()
	{
		// one month runtime
		if (simVars.timecontrol.max_simulation_time == -1)
		{
			simVars.timecontrol.max_simulation_time = 31*60*60*24;

			// 144 h
			//simVars.timecontrol.max_simulation_time = 144*60*60;

			// 200 h
			simVars.timecontrol.max_simulation_time = 200*60*60;
	//		simVars.timecontrol.max_simulation_time = 1;
		}

		simVars.next_output_dt = 0;


		if (simVars.timecontrol.current_timestep_size <= 0)
		{
			// TRY to guess optimal time step size

			// time step size
			if (sphConfig->spat_num_lat < 256)
			{
	//			simVars.viscosity2 = 1e5;
				simVars.timecontrol.current_timestep_size = 0.002*simVars.earth_radius/(double)sphConfig->spat_num_lat;
			}
			else
			{
	//			simVars.viscosity2 = 1e5;
				simVars.timecontrol.current_timestep_size = 0.001*simVars.earth_radius/(double)sphConfig->spat_num_lat;
			}
			//simVars.timecontrol.current_timestep_size = 30;
		}

		if (simVars.output_dt <= 0)
		{
			simVars.output_dt = 60*30;	// output every 1/2 hour
		}


		if (simVars.use_rexi == 1)
		{
			// Override for REXI
			if (simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}

			simVars.use_nonlinear_equations = 0;
		}

		if (simVars.benchmark_scenario_id == 0)
		{
			setup_initial_conditions_gaussian();
		}
		else if (simVars.benchmark_scenario_id == 1)
		{
			benchmarkGalewsky.setup_initial_h(prog_h);
//			prog_h.spat_set_zero();
			benchmarkGalewsky.setup_initial_h_add_bump(prog_h);

			benchmarkGalewsky.setup_initial_u(prog_u);
			benchmarkGalewsky.setup_initial_v(prog_v);
		}
		else if (simVars.benchmark_scenario_id == 2 || simVars.benchmark_scenario_id == 3)
		{
			// Non-dimensional stuff
#if 0
			if (simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}
#endif
			simVars.h0 = 1;
			simVars.gravitation = 1;
			simVars.earth_radius = 1;

			simVars.use_nonlinear_equations = 0;

			if (simVars.benchmark_scenario_id == 2)
			{
				setup_initial_conditions_gaussian(0);
#if 0
	//			prog_h.spat_set_value(simVars.h0);
				prog_v.spat_set_zero();
				prog_u.spat_update_lambda_cogaussian_grid(
						[](double lon, double comu, double &io_data)
						{
							io_data = 1*comu;
						}
					);
				//prog_u.spat_set_zero();
#endif
			}
			else if (simVars.benchmark_scenario_id == 3)
			{
				setup_initial_conditions_gaussian(M_PI/3.0);
//				setup_initial_conditions_gaussian(-M_PI/3.0);
			}
		}


		bool with_coriolis = false;

		std::cout << "Using time step size dt = " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << "Running simulation until t_end = " << simVars.timecontrol.max_simulation_time << std::endl;
		std::cout << "Parameters:" << std::endl;
		std::cout << " + gravity: " << simVars.gravitation << std::endl;
		std::cout << " + earth_radius: " << simVars.earth_radius << std::endl;
		std::cout << " + average height: " << simVars.h0 << std::endl;
		std::cout << " + coriolis_omega: " << simVars.coriolis_omega << std::endl;
		std::cout << " + viscosity D2: " << simVars.viscosity2 << std::endl;
		std::cout << " + use_nonlinear: " << simVars.use_nonlinear_equations << std::endl;
		std::cout << " + Coriolis: " << with_coriolis << std::endl;
		std::cout << std::endl;
		std::cout << " + Benchmark scenario id: " << simVars.benchmark_scenario_id << std::endl;
		std::cout << " + Use robert functions: " << simVars.use_robert_functions << std::endl;
		std::cout << " + Use REXI: " << simVars.use_rexi << std::endl;
		std::cout << " + REXI h: " << simVars.rexi_h << std::endl;
		std::cout << " + REXI M: " << simVars.rexi_M << std::endl;
		std::cout << " + REXI use half poles: " << simVars.rexi_use_half_poles << std::endl;
		std::cout << " + REXI additional modes: " << simVars.rexi_use_extended_modes << std::endl;
		std::cout << std::endl;
		std::cout << " + timestep size: " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << " + output timestep size: " << simVars.output_dt << std::endl;

		std::cout << std::endl;

		fdata.spat_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu;
				}
			);


		// This class is only used in case of added modes
		SPHConfig sphConfigRexiAddedModes;

		// Pointer to SPH configuration for REXI computations
		SPHConfig *sphConfigRexi = nullptr;

		if (simVars.use_rexi == false)
		{
			simVars.timecontrol.current_simulation_time = 0;
			while (simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time)
			{
				if (simVars.timecontrol.current_simulation_time >= simVars.next_output_dt)
				{
					std::cout << std::endl;
					write_output();

					simVars.next_output_dt += simVars.output_dt;
					if (simVars.next_output_dt < simVars.timecontrol.current_simulation_time)
						simVars.next_output_dt = simVars.timecontrol.current_simulation_time;
				}

				double o_dt;
				timestepping.run_rk_timestep(
						this,
						&AppSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						prog_h, prog_u, prog_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.timecontrol.timestepping_runge_kutta_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);

				std::cout << "." << std::flush;

				if (prog_h.isAnyNaNorInf())
				{
					std::cerr << "Instability detected (NaN value in H)" << std::endl;
					assert(false);
					exit(1);
				}

				simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
			}

			std::cout << std::endl;
			write_output();
			std::cout << std::endl;
		}
		else
		{
			if (simVars.rexi_use_extended_modes == 0)
			{
				sphConfigRexi = sphConfig;
			}
			else
			{
				// Add modes only along latitude since these are the "problematic" modes
				sphConfigRexiAddedModes.setupAdditionalModes(
						sphConfig,
						simVars.rexi_use_extended_modes,	// TODO: Extend SPH wrapper to also support m != n to set this guy to 0
						simVars.rexi_use_extended_modes
				);
				sphConfigRexi = &sphConfigRexiAddedModes;
			}

			SPHSolver<double> sphSolver;
			sphSolver.setup(sphConfig, 4);

			rexi.setup(simVars.rexi_h, simVars.rexi_M, 0, simVars.rexi_use_half_poles);

			std::cout << "REXI poles: " << rexi.alpha.size() << std::endl;

			SPHData tmp_prog_phi(sphConfigRexi);
			SPHData tmp_prog_u(sphConfigRexi);
			SPHData tmp_prog_v(sphConfigRexi);

			SPHData accum_prog_phi(sphConfigRexi);
			SPHData accum_prog_u(sphConfigRexi);
			SPHData accum_prog_v(sphConfigRexi);

			std::vector<SWERexiSPHRobert> rexiSPHRobert_vector;
			std::vector<SWERexiSPH> rexiSPH_vector;

			bool use_rexi_preallocaation = true;

			if (use_rexi_preallocaation)
			{
				if (simVars.use_robert_functions)
					rexiSPHRobert_vector.resize(rexi.alpha.size());
				else
					rexiSPH_vector.resize(rexi.alpha.size());


				for (int i = 0; i < rexi.alpha.size(); i++)
				{
					std::complex<double> &alpha = rexi.alpha[i];
					std::complex<double> &beta = rexi.beta_re[i];
					//beta *= 0.9999;

					if (simVars.use_robert_functions)
					{
						rexiSPHRobert_vector[i].setup(
								sphConfigRexi,
								alpha,
								beta,
								simVars.earth_radius,
								simVars.coriolis_omega, //*inv_sqrt_avg_geopo,
								simVars.h0*simVars.gravitation,
								simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
								with_coriolis
						);
					}
					else
					{
						rexiSPH_vector[i].setup(
								sphConfigRexi,
								alpha,
								beta,
								simVars.earth_radius,
								simVars.coriolis_omega, //*inv_sqrt_avg_geopo,
								simVars.h0*simVars.gravitation,
								simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
								with_coriolis
						);
					}
				}
			}

			simVars.timecontrol.current_simulation_time = 0;
			while (simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time)
			{
				if (simVars.timecontrol.current_simulation_time >= simVars.next_output_dt)
				{
					write_output();

					simVars.next_output_dt += simVars.output_dt;
					if (simVars.next_output_dt < simVars.timecontrol.current_simulation_time)
						simVars.next_output_dt = simVars.timecontrol.current_simulation_time;
				}


				{
					// convert to geopotential
					SPHData prog_phi_rexi(sphConfigRexi);
					SPHData prog_u_rexi(sphConfigRexi);
					SPHData prog_v_rexi(sphConfigRexi);

					if (simVars.rexi_use_extended_modes == 0)
					{
						prog_phi_rexi = prog_h*simVars.gravitation;
						prog_u_rexi = prog_u;
						prog_v_rexi = prog_v;
					}
					else
					{
						(prog_h*simVars.gravitation).spec_copyToDifferentModes(prog_phi_rexi);
						prog_u.spec_copyToDifferentModes(prog_u_rexi);
						prog_v.spec_copyToDifferentModes(prog_v_rexi);
					}

					accum_prog_phi.spat_set_zero();
					accum_prog_u.spat_set_zero();
					accum_prog_v.spat_set_zero();

					for (int i = 0; i < rexi.alpha.size(); i++)
					{
						std::complex<double> &alpha = rexi.alpha[i];
						std::complex<double> &beta = rexi.beta_re[i];

						if (simVars.use_robert_functions)
						{
							if (use_rexi_preallocaation)
							{
								rexiSPHRobert_vector[i].solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
							else
							{
								SWERexiSPHRobert rexiSPHRobert;
								rexiSPHRobert.setup(
										sphConfig,
										alpha,
										beta,
										simVars.earth_radius,
										simVars.coriolis_omega, //*inv_sqrt_avg_geopo,
										simVars.h0*simVars.gravitation,
										simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
										with_coriolis
								);

								rexiSPHRobert.solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
						}
						else
						{
							if (use_rexi_preallocaation)
							{
								rexiSPH_vector[i].solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
							else
							{
								SWERexiSPH rexiSPH;
								rexiSPH.setup(
										sphConfig,
										alpha,
										beta,
										simVars.earth_radius,
										simVars.coriolis_omega, //*inv_sqrt_avg_geopo,
										simVars.h0*simVars.gravitation,
										simVars.timecontrol.current_timestep_size, //*sqrt_avg_geopo
										with_coriolis
								);

								rexiSPH.solve(
										prog_phi_rexi, prog_u_rexi, prog_v_rexi,
										tmp_prog_phi, tmp_prog_u, tmp_prog_v,
										opComplex
									);
							}
						}

						accum_prog_phi += tmp_prog_phi;
						accum_prog_u += tmp_prog_u;
						accum_prog_v += tmp_prog_v;
					}

					if (simVars.rexi_use_extended_modes == 0)
					{
						prog_h = accum_prog_phi*(1.0/simVars.gravitation);
						prog_u = accum_prog_u;
						prog_v = accum_prog_v;
					}
					else
					{
						(accum_prog_phi*(1.0/simVars.gravitation)).spec_copyToDifferentModes(prog_h);
						accum_prog_u.spec_copyToDifferentModes(prog_u);
						accum_prog_v.spec_copyToDifferentModes(prog_v);
					}
				}


				std::cout << "." << std::flush;

				/*
				 * Add implicit viscosity
				 */
				if (simVars.viscosity2 != 0)
				{
					double scalar = simVars.viscosity2*simVars.timecontrol.current_timestep_size;
					double r = simVars.earth_radius;

					/*
					 * (1-dt*visc*D2)p(t+dt) = p(t)
					 */
					prog_h = prog_h.spec_solve_helmholtz(1.0, -scalar, r);
					prog_u = prog_u.spec_solve_helmholtz(1.0, -scalar, r);
					prog_v = prog_v.spec_solve_helmholtz(1.0, -scalar, r);
				}

				if (prog_h.isAnyNaNorInf())
				{
					std::cerr << "Instability detected (NaN value in H)" << std::endl;
					assert(false);
					exit(1);
				}

				simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
			}
			write_output();
			std::cout << std::endl;
		}
	}



	// Main routine for method to be used in case of finite differences
	void p_run_euler_timestep_update(
			const SPHData &i_h,	///< prognostic variables
			const SPHData &i_u,	///< prognostic variables
			const SPHData &i_v,	///< prognostic variables

			SPHData &o_h_t,	///< time updates
			SPHData &o_u_t,	///< time updates
			SPHData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (!simVars.use_nonlinear_equations)
		{
			if (simVars.use_robert_functions == 0)
			{
				// linear equations
				o_h_t = -(op.div_lon(i_u)+op.div_lat(i_v))*(simVars.h0/simVars.earth_radius);

				o_u_t = -op.grad_lon(i_h)*(simVars.gravitation/simVars.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.gravitation/simVars.earth_radius);

				if (simVars.coriolis_omega != 0)
				{
					o_u_t += f(i_v);
					o_v_t -= f(i_u);
				}
			}
			else
			{
				// use Robert functions for velocity
				// linear equations
				o_h_t = -(op.robert_div_lon(i_u)+op.robert_div_lat(i_v))*(simVars.h0/simVars.earth_radius);

				o_u_t = -op.robert_grad_lon(i_h)*(simVars.gravitation/simVars.earth_radius);
				o_v_t = -op.robert_grad_lat(i_h)*(simVars.gravitation/simVars.earth_radius);

				if (simVars.coriolis_omega != 0)
				{
					o_u_t += f(i_v);
					o_v_t -= f(i_u);
				}
			}
		}
		else
		{
			assert(simVars.earth_radius > 0);
			assert(simVars.gravitation);

			/*
			 * Height
			 */
			// non-linear equations
			o_h_t = -(op.div_lon(i_h*i_u)+op.div_lat(i_h*i_v))*(1.0/simVars.earth_radius);

			/*
			 * Velocity
			 */
			// linear terms
			o_u_t = -op.grad_lon(i_h)*(simVars.gravitation/simVars.earth_radius);
			o_v_t = -op.grad_lat(i_h)*(simVars.gravitation/simVars.earth_radius);

			if (simVars.coriolis_omega != 0)
			{
				o_u_t += f(i_v);
				o_v_t -= f(i_u);
			}

			// non-linear terms
			o_u_t -= (i_u*op.grad_lon(i_u) + i_v*op.grad_lat(i_u))*(1.0/simVars.earth_radius);
			o_v_t -= (i_u*op.grad_lon(i_v) + i_v*op.grad_lat(i_v))*(1.0/simVars.earth_radius);
		}

		if (simVars.viscosity2 != 0)
		{
			double scalar = simVars.viscosity2/(simVars.earth_radius*simVars.earth_radius);

			o_h_t += op.laplace(i_h)*scalar;
			o_u_t += op.laplace(i_u)*scalar;
			o_v_t += op.laplace(i_v)*scalar;
		}
	}
};


#endif /* SRC_TESTSWE_HPP_ */
