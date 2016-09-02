/*
 * SimVars.hpp
 *
 *  Created on: ?? Aug 2016
 *      Author: martin
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cassert>

#include <sph/SPHConfig.hpp>
#include <sph/SPHOperators.hpp>
#include <sph/SPHData.hpp>
#include <sph/SPHHelper.hpp>
#include "SimVars.hpp"

//#include "../3rd_party/Adaptive-Integrator/AdaptiveIntegrator.hpp"
#include <sph/SPHConfig.hpp>
#include <sweet/MemBlockAlloc.hpp>

#include "AppSWE.hpp"
#include "AppTestSPHOperators.hpp"
#include "AppTestSPHOperatorsComplex.hpp"
#include "AppTestSPHSolver.hpp"
#include "AppTestSPHSolverComplex.hpp"
#include "AppOutputSphericalHarmonics.hpp"



SimVars simVars;


int setup(
		int i_argc,
		const char *i_argv[],
		SPHConfig *io_sphConfig
)
{
	if (i_argc < 2)
	{
		std::cerr << "Usage: " << i_argv[0] << " T[n]                         # Use T notation" << std::endl;
		std::cerr << "   OR: " << i_argv[0] << " [n] [m]                      # set modes" << std::endl;
		std::cerr << "   OR: " << i_argv[0] << " [n] [m] [res_lon] [res_lat]  # set modes and spatial resolution" << std::endl;
		return -1;
	}

	if (i_argv[1][0] == 'T')
	{
		int spec_order = atoi(&(i_argv[1][1]));
		std::cout << "Using T=" << spec_order << std::endl;

		/*
		 * http://www.ecmwf.int/en/what-horizontal-resolution-data
		 *
		 * Truncation at wave number T
		 *
		 * Hack an Jakob: Eq. (3.2)
		 * Use Triangular truncation M=N=K
		 */

		simVars.spec_res_n = spec_order;
		simVars.spec_res_m = spec_order;

		for (int i = 2; i < i_argc; i++)
		{
			if (i_argv[i][0] == 'F')
			{
				simVars.coriolis_omega = atof(&(i_argv[i][1]));
			}
			else if (i_argv[i][0] == 'P')
			{
				simVars.program_id = atoi(&(i_argv[i][1]));
			}
			else if (i_argv[i][0] == 'N')
			{
				simVars.use_nonlinear_equations = atoi(&(i_argv[i][1]));
			}
			else if (i_argv[i][0] == 'M')
			{
				simVars.timestepping_method = atoi(&(i_argv[i][1]));
			}
			else if (i_argv[i][0] == 'R')
			{
				simVars.timecontrol.current_timestep_size = atof(&(i_argv[i][1]));
			}
		}
	}
	else
	{
		if (i_argc > 1)
			simVars.spec_res_n = atoi(i_argv[1]);

		if (i_argc > 2)
			simVars.spec_res_m = atoi(i_argv[2]);

		if (i_argc > 3)
			simVars.spat_res_lon = atoi(i_argv[3]);

		if (i_argc > 4)
			simVars.spat_res_lat = atoi(i_argv[4]);
	}

	if (simVars.spec_res_m == -1)
		simVars.spec_res_m = simVars.spec_res_n;

#if 0
	// let shtns choose the spatial resolution

	/*
	 * See Hack&Jakob, Eq. (3.10)
	 */
	if (spat_res_lon == -1)
		spat_res_lon = 3*spec_res_m+1;


	/*
	 * See Hack&Jakob, Eq. (3.16)
	 */
	if (spat_res_lat == -1)
		spat_res_lat = (3*spec_res_n+1)/2;

	auto ceilpow2 = [](int x) -> int
	{
		int power = 1;
		while(power < x)
		    power *= 2;

		return power;
	};

	auto ceil16 = [](int x) -> int
	{
		int value = 16;
		while(value < x)
		    value = value+16;

		return value;
	};

	spat_res_lon = ceilpow2(spat_res_lon);
	spat_res_lat = ceil16(spat_res_lat);
#endif

	std::cout << "Running setup..." << std::endl;

	if (simVars.spat_res_lon <= 0)
	{
		io_sphConfig->setupAuto(
				simVars.spec_res_n, simVars.spec_res_m,		// spectral (lon/lat)
				&simVars.spat_res_lon, &simVars.spat_res_lat	// spatial (lon/lat)
		);
	}
	else
	{
		io_sphConfig->setup(
				simVars.spec_res_n, simVars.spec_res_m,		// spectral (lon/lat)
				simVars.spat_res_lon, simVars.spat_res_lat	// spatial (lon/lat)
			);
	}

	std::cout << "Using spectral modes m=" << simVars.spec_res_m << " and n=" << simVars.spec_res_n << std::endl;
	std::cout << "Requested spatial resolution " << simVars.spat_res_lon << " x " << simVars.spat_res_lat << std::endl;
	std::cout << "Using spatial resolution " << io_sphConfig->spat_num_lon << " x " << io_sphConfig->spat_num_lat << std::endl;
	std::cout << "Using scalar for Coriolis " << simVars.coriolis_omega << std::endl;

	return 0;
}


int main(
		int i_argc,
		const char *i_argv[]
)
{

	/*
	 * Initialize NUMA block allocator
	 */
	MemBlockAlloc numaBlockAlloc;

	SPHConfig sphConfig;

	int retval = setup(i_argc, i_argv, &sphConfig);
	if (retval != 0)
		return retval;

	if (simVars.program_id == 0)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppSWE" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppSWE test_swe(&sphConfig, simVars);
		test_swe.run();
	}
	else if (simVars.program_id == 1)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppTestSPHOperators" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppTestSPHOperators test_operators; test_operators.run(&sphConfig);
	}


	else if (simVars.program_id == 2)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppTestSPHOperatorsComplex" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppTestSPHOperatorsComplex test_operators; test_operators.run(&sphConfig);
	}

	else if (simVars.program_id == 3)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppTestSPHSolvers" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppTestSPHSolvers test_operators; test_operators.run(&sphConfig);
	}

	else if (simVars.program_id == 4)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppTestSPHSolversComplex" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppTestSPHSolversComplex test_operators; test_operators.run(&sphConfig);
	}

	else if (simVars.program_id == 5)
	{
		std::cout << "************************************************" << std::endl;
		std::cout << "AppOutputSphericalHarmonics" << std::endl;
		std::cout << "************************************************" << std::endl;

		AppOutputSphericalHarmonics output_spherical_harmonics;
		output_spherical_harmonics.run(&sphConfig);
	}

	sphConfig.shutdown();
}

