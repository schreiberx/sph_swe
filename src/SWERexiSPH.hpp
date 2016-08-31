/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEREXISPH_HPP_
#define SRC_SWEREXISPH_HPP_

#include <sph/SPHSolver.hpp>
#include <sph/SPHOperatorsComplex.hpp>
#include <sph/SPHConfig.hpp>

class SWERexiSPH
{
	/// SPH configuration
	SPHConfig *sphConfig;

	/// Solver for given alpha
	SPHSolver<std::complex<double>> sphSolver;

	/// scalar infront of RHS
	std::complex<double> rhs_scalar;


	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

	/// earth radius
	double r;

	/// Coriolis omega
	double coriolis_omega;

	/// 2*\Omega
	double two_omega;

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
			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,
			double i_radius,
			double i_coriolis_omega
	)
	{
		alpha = i_alpha;
		beta = i_beta;

		const std::complex<double> &alpha = i_alpha;
		r = i_radius;
		coriolis_omega = i_coriolis_omega;
		two_omega = 2.0*i_coriolis_omega;

		sphConfig = i_sphConfig;

		sphSolver.setup(sphConfig, 4);

		sphSolver.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
		sphSolver.solver_component_rexi_z2(	alpha*alpha*two_omega*two_omega, r);
		sphSolver.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
		sphSolver.solver_component_rexi_z4(	-alpha*two_omega, r);
		sphSolver.solver_component_rexi_z5(	1.0/alpha*two_omega*two_omega*two_omega, r);
		sphSolver.solver_component_rexi_z6(	2.0*alpha*two_omega*two_omega, r);
		sphSolver.solver_component_rexi_z7(	-alpha*alpha, r);
		sphSolver.solver_component_rexi_z8(	-two_omega*two_omega, r);
	}



	/**
	 * Solve a REXI time step for the given initial conditions
	 */
	void solve(
			const SPHData &i_phi0,
			const SPHData &i_u0,
			const SPHData &i_v0,

			SPHData &o_phi,
			SPHData &o_u,
			SPHData &o_v,

			SPHOperatorsComplex &opComplex,

			double i_timestep_size
	)
	{
#if 0
		SPHDataComplex phi0(i_phi0);
		SPHDataComplex u0(i_u0);
		SPHDataComplex v0(i_v0);

		SPHDataComplex div0(prog_h(opComplex.div(u0, v0));
		SPHDataComplex vort0(opComplex.vort(u0, v0));

//		SPHDataComplex d0(op.div_lat(i_u0) + op.div_lon(i_v0));
//		SPHDataComplex eta0(op.div_lat(i_v0) - op.div_lon(i_u0));

		SPHDataComplex rhs = d0 - (1.0/alpha)*op.mu(eta0) + alpha*i_phi0 + (1.0/alpha)*op.mu2(i_phi0) + (1.0/alpha)*Fc;
		SPHDataComplex b = alpha*alpha*rhs + op.mu2(rhs);
#endif
	}
};


#endif /* SRC_SWEREXISPH_HPP_ */
