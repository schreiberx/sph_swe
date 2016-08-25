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
