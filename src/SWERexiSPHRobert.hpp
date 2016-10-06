/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEREXISPHROBERT_HPP_
#define SRC_SWEREXISPHROBERT_HPP_

#include <sph/SPHSolverComplex.hpp>
#include <sph/SPHOperatorsComplex.hpp>
#include <sph/SPHConfig.hpp>



/**
 * REXI solver for SWE based on Robert function formulation
 */
class SWERexiSPHRobert
{
	/// SPH configuration
	SPHConfig *sphConfig;

	/// Solver for given alpha
	SPHSolverComplex<std::complex<double>> sphSolverPhi;
	SPHSolverComplex<std::complex<double>> sphSolverVel;

	/// scalar infront of RHS
	std::complex<double> rhs_scalar;

	/// REXI alpha
	std::complex<double> alpha;

	/// REXI beta
	std::complex<double> beta;

	bool include_coriolis_effect;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double ir;

	/// Coriolis omega
	double coriolis_omega;

	/// 2*\Omega
	double two_omega;

	/// Average geopotential
	double avg_geopotential;

public:
	SWERexiSPHRobert()	:
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
			double i_coriolis_omega,
			double i_avg_geopotential,
			double i_timestep_size,
			bool i_include_coriolis_effect = true
	)
	{
		include_coriolis_effect = i_include_coriolis_effect;
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta/timestep_size;

		r = i_radius;
		ir = 1.0/r;

		coriolis_omega = i_coriolis_omega;
		two_omega = 2.0*coriolis_omega;
		avg_geopotential = i_avg_geopotential;

		sphConfig = i_sphConfig;

		sphSolverPhi.setup(sphConfig, 1);
		// alpha^2
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha), r);

		// -avg_geopotential*laplace
		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential, r);

		// (-avg_geopotential)*(-2)*Phi
		sphSolverPhi.solver_component_rexi_z1(	-2.0*avg_geopotential/(r*r)*(-2.0)	, r);

		if (include_coriolis_effect)
		{
		}

		// not necessary yet
		sphSolverVel.setup(sphConfig, 2);
		if (include_coriolis_effect)
		{
		}
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

			SPHOperatorsComplex &op
	)
	{
#if 0
		SPHDataComplex mu(i_phi0.sphConfig);
		mu.spat_update_lambda_gaussian_grid(
				[&](double lon, double mu, std::complex<double> &o_data)
				{
					o_data = mu;
				}
			);
#endif

		SPHDataComplex phi0(i_phi0);
		SPHDataComplex u0(i_u0);
		SPHDataComplex v0(i_v0);

#if 0
		phi0 *= (1.0/timestep_size);
		u0 *= (1.0/timestep_size);
		v0 *= (1.0/timestep_size);
#endif

		SPHDataComplex div0(ir*op.robert_div(u0, v0));
		SPHDataComplex eta0(ir*op.robert_vort(u0, v0));

		SPHDataComplex phi(sphConfig);
		SPHDataComplex u(sphConfig);
		SPHDataComplex v(sphConfig);

		if (include_coriolis_effect)
		{
			assert(false);
			std::cerr << "NOT YET IMPLEMENTED!" << std::endl;
		}
		else
		{
#if 1

			SPHDataComplex rhs = avg_geopotential*div0 + alpha*i_phi0;
			phi = rhs.spec_solve_helmholtz(alpha*alpha, -avg_geopotential, r);
//			phi = sphSolverPhi.solve(rhs);

			u = (1.0/alpha) * (u0 + ir*op.robert_grad_lon(phi));
			v = (1.0/alpha) * (v0 + ir*op.robert_grad_lat(phi));

#else

			SPHDataComplex rhs = avg_geopotential*div0 + alpha*i_phi0;
			phi = rhs.spec_solve_helmholtz(alpha*alpha, -avg_geopotential, r);

			// same solver, but without solving inverse problem
			u = (1.0/alpha) * (u0 + ir*op.grad_lon(phi));
			v = (1.0/alpha) * (v0 + ir*op.grad_lat(phi));

#endif
		}

		phi *= beta;
		u *= beta;
		v *= beta;

		phi.spat_RealToSPHData(o_phi);
		u.spat_RealToSPHData(o_u);
		v.spat_RealToSPHData(o_v);
	}
};


#endif /* SRC_SWEREXISPHROBERT_HPP_ */
