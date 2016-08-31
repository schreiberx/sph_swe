/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWE_REXI_SPH_HPP_
#define SRC_SWE_REXI_SPH_HPP_

#include <sph/SPHSolver.hpp>
#include <sph/SPHConfig.hpp>

class SWERexiSPH
{
	SPHConfig *sphConfig;

	SPHSolver<std::complex<double>> sphSolver;

public:
	SWERexiSPH()	:
		sphConfig(nullptr)
	{
	}


	/**
	 * Setup the SWE REXI solver with SPH
	 */
	void setup(
			SPHConfig *i_sphConfig,
			const std::complex<double> &i_alpha
	)
	{
		sphConfig = i_sphConfig;

		sphSolver.setup(sphConfig, 4);


		sphSolver.solver_component_rexi_z1(std::pow(i_alpha, 4.0));
	}
};


#endif /* SRC_SWE_REXI_SPH_HPP_ */
