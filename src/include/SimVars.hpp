/*
 * SimVars.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SIMVARS_HPP_
#define SRC_SIMVARS_HPP_


struct SimVars
{
	int spat_res_lon = -1;
	int spat_res_lat = -1;

	int spec_res_m = -1;
	int spec_res_n = -1;

	// constants from Galwesky et al. paper
	double coriolis_omega = 7.292e-5;
	double gravitation = 9.80616;
	double earth_radius = 6.37122e6;

	int program_id = 0;

//	double viscosity2 = 1e-8;
	double viscosity2 = 1e5;
//	double viscosity2 = 1e2;

	// output data every 0.01 simulation seconds
	double output_dt = 0.2;
	double next_output_dt = 0;

	double rexi_h = 0.2;
	int rexi_M = 128;
	bool rexi_use_half_poles = true;

	/*
	 * (phi, u, v) = (geopotential, velocity u, velocity v)
	 *
	 * 0: (phi, u, v)
	 * 1: (phi, U, V) with U=cos(phi)u and V=cos(phi)v
	 */
	int use_robert_functions = 0;


	/*
	 * Use more modes for REXI time step approximation
	 */
	int rexi_use_extended_modes = 0;


	/*
	 * false: Use standard time stepping
	 * true: Use REXI
	 */
	bool use_rexi = false;


	/*
	 * Benchmark scenario / initial conditions, etc.
	 */
	int benchmark_scenario_id = 0;
//	int timestepping_method = 0;

	bool use_nonlinear_equations = true;
	double h0 = 10000.0;

	struct TimeControl
	{
		int timestepping_runge_kutta_order = 4;

		double current_timestep_size = 0;
		double current_simulation_time = 0;
		double max_simulation_time = -1;

	} timecontrol;
};




#endif /* SRC_SIMVARS_HPP_ */
