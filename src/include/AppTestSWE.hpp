/*
 * TestSWE.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_TESTSWE_HPP_
#define SRC_TESTSWE_HPP_

#include "../../3rd_party/Adaptive-Integrator/AdaptiveIntegrator.hpp"
#include <BenchmarkGalewsky.hpp>
#include <TimesteppingRK.hpp>
#include <rexi/REXI.hpp>
#include <GaussQuadrature.hpp>
#include <sph/SPHDataComplex.hpp>



class AppTestSWE
{
public:
	SimVars &simVars;
	SPHOperators op;
	SPHConfig *sphConfig;

	// Runge-Kutta stuff
	TimesteppingRK timestepping;


	SPHData prog_h;
	SPHData prog_u;
	SPHData prog_v;
	SPHData fdata;

	// output data every 0.01 simulation seconds
	double output_dt = 0.2;
	double next_output_dt = 0;

	int benchmark_id = 1;

	BenchmarkGalewsky benchmarkGalewsky;

	REXI rexi;



public:
	AppTestSWE(
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
		op.setup(i_sphConfig);
	}



	void write_output()
	{
		char buffer[1024];

		sprintf(buffer, "prog_h_t%010.4f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (benchmark_id == 1)
			prog_h.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_h.spat_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_u_t%010.4f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (benchmark_id == 1)
			prog_u.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_h.spat_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_v_t%010.4f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		if (benchmark_id == 1)
			prog_v.spat_write_file_lon_pi_shifted(buffer);
		else
			prog_h.spat_write_file(buffer);
		std::cout << buffer << std::endl;

		sprintf(buffer, "prog_eta_t%010.4f.csv", simVars.timecontrol.current_simulation_time/(60*60));
		SPHData vort = op.vort(prog_u, prog_v)/simVars.earth_radius;
		if (benchmark_id == 1)
			vort.spat_write_file_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			vort.spat_write_file(buffer);
		std::cout << buffer << std::endl;
	}



	void setup_initial_conditions_gaussian()
	{
		double exp_fac = 10.0;

#if 0
		double center_lon = M_PI/3;
		double center_lat = M_PI/5;
#else
		double center_lon = 0;//M_PI;
		double center_lat = M_PI/3;
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

		// one month runtime
		simVars.timecontrol.max_simulation_time = 31*60*60*24;

		// 144 h
		//simVars.timecontrol.max_simulation_time = 144*60*60;

		// 200 h
		simVars.timecontrol.max_simulation_time = 200*60*60;
//		simVars.timecontrol.max_simulation_time = 1;

		output_dt = 60*30;	// output every 1/2 hour
//		output_dt = 0;

		next_output_dt = 0;

		std::cout << "Using time step size dt = " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << "Running simulation until t_end = " << simVars.timecontrol.max_simulation_time << std::endl;
		std::cout << "Parameters:" << std::endl;
		std::cout << " + gravity: " << simVars.gravitation << std::endl;
		std::cout << " + earth_radius: " << simVars.earth_radius << std::endl;
		std::cout << " + coriolis_omega: " << simVars.coriolis_omega << std::endl;
		std::cout << " + viscosity D2: " << simVars.viscosity2 << std::endl;
		std::cout << std::endl;

		fdata.spat_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu;
				}
			);

		if (benchmark_id == 0)
		{
			setup_initial_conditions_gaussian();
		}
		else if (benchmark_id == 1)
		{
			benchmarkGalewsky.setup_initial_h(prog_h);
//			prog_h.spat_set_zero();
			benchmarkGalewsky.setup_initial_h_add_bump(prog_h);

			benchmarkGalewsky.setup_initial_u(prog_u);
			benchmarkGalewsky.setup_initial_v(prog_v);
		}

		int timestepping_method = 0;

		if (timestepping_method == 0)
		{
			simVars.timecontrol.current_simulation_time = 0;
			while (simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time)
			{
				if (simVars.timecontrol.current_simulation_time >= next_output_dt)
				{
					write_output();

					next_output_dt += output_dt;
					if (next_output_dt < simVars.timecontrol.current_simulation_time)
						next_output_dt = simVars.timecontrol.current_simulation_time;
				}

				double o_dt;
				timestepping.run_rk_timestep(
						this,
						&AppTestSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
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
			write_output();
			std::cout << std::endl;
		}
		else
		{
			rexi.setup(0.2, 256);

			// convert to geopotential
			SPHData phi0 = prog_h*simVars.gravitation;
			SPHData u0 = prog_u;
			SPHData v0 = prog_v;

			SPHData h_accum(prog_h.sphConfig);
			SPHData u_accum(prog_h.sphConfig);
			SPHData v_accum(prog_h.sphConfig);
			h_accum.spat_set_zero();
			u_accum.spat_set_zero();
			v_accum.spat_set_zero();

			SPHDataComplex div0(op.div(u0, v0));
			SPHDataComplex vort0(op.vort(u0, v0));
			SPHDataComplex phi0c(phi0);
			SPHDataComplex fdatac(fdata);

			//SPHDataComplex Fc;

			for (int i = 0; i < rexi.beta_re.size(); i++)
			{
				std::complex<double> alpha = rexi.alpha[i];
				std::complex<double> beta = rexi.beta_re[i];

				SPHDataComplex rhs = div0 - (1.0/alpha)*fdatac*vort0 + (alpha+fdatac*fdatac*(1.0/alpha))*phi0c;
				//rhs -= (1.0/alpha)*Fc;
				//h_accum += ...
			}
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

			double &o_dt,			///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (!simVars.use_nonlinear_equations)
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
