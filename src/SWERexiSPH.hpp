/*
 * SWE_REXI_SPH.hpp
 *
 *  Created on: 30 Aug 2016
 *      Author: martin
 */

#ifndef SRC_SWEREXISPH_HPP_
#define SRC_SWEREXISPH_HPP_

#include <sph/SPHSolverComplex.hpp>
#include <sph/SPHOperatorsComplex.hpp>
#include <sph/SPHConfig.hpp>



class SWERexiSPH
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

		sphSolverPhi.setup(sphConfig, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);

		if (include_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
			sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
			sphSolverPhi.solver_component_rexi_z4(	-avg_geopotential*alpha*two_omega, r);
			sphSolverPhi.solver_component_rexi_z5(	avg_geopotential/alpha*two_omega*two_omega*two_omega, r);
			sphSolverPhi.solver_component_rexi_z6(	avg_geopotential*2.0*two_omega*two_omega, r);
		}
		sphSolverPhi.solver_component_rexi_z7(	-avg_geopotential*alpha*alpha, r);
		if (include_coriolis_effect)
		{
			sphSolverPhi.solver_component_rexi_z8(	-avg_geopotential*two_omega*two_omega, r);
		}

		sphSolverVel.setup(sphConfig, 2);
		sphSolverVel.solver_component_rexi_z1(	alpha*alpha, r);
		if (include_coriolis_effect)
		{
			sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
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
		SPHDataComplex mu(i_phi0.sphConfig);
		mu.spat_update_lambda_gaussian_grid(
				[&](double lon, double mu, std::complex<double> &o_data)
				{
					o_data = mu;
				}
			);

		SPHDataComplex phi0(i_phi0);
		SPHDataComplex u0(i_u0);
		SPHDataComplex v0(i_v0);


		SPHDataComplex div0(ir*op.div(u0, v0));
		SPHDataComplex eta0(ir*op.vort(u0, v0));

		SPHDataComplex phi(sphConfig);
		SPHDataComplex u(sphConfig);
		SPHDataComplex v(sphConfig);

		if (include_coriolis_effect)
		{

			SPHDataComplex Fck = two_omega*ir*op.grad_lon(mu)*(-(alpha*alpha*u0 - two_omega*two_omega*op.mu2(u0)) + 2.0*alpha*two_omega*op.mu(v0));
			SPHDataComplex foo = div0 - two_omega*(1.0/alpha)*op.mu(eta0) + alpha*i_phi0 + two_omega*two_omega*(1.0/alpha)*op.mu2(i_phi0);
			SPHDataComplex rhs = alpha*alpha*foo + two_omega*two_omega*op.mu2(foo) + (1.0/alpha)*Fck;

			phi = sphSolverPhi.solve(rhs);

			SPHDataComplex a = u0 + ir*op.grad_lon(phi);
			SPHDataComplex b = v0 + ir*op.grad_lat(phi);

			SPHDataComplex rhsa = alpha*a - op.mu(b);
			SPHDataComplex rhsb = op.mu(a) + alpha*b;

			u = sphSolverVel.solve(rhsa);
			v = sphSolverVel.solve(rhsb);
		}
		else
		{
			SPHDataComplex foo = div0 + alpha*i_phi0;
			SPHDataComplex rhs = alpha*alpha*foo;

			phi = sphSolverPhi.solve(rhs);

			SPHDataComplex a = u0 + ir*op.grad_lon(phi);
			SPHDataComplex b = v0 + ir*op.grad_lat(phi);

			SPHDataComplex rhsa = alpha*a;
			SPHDataComplex rhsb = alpha*b;

			u = sphSolverVel.solve(rhsa);
			v = sphSolverVel.solve(rhsb);
		}


		phi *= beta;
		u *= beta;
		v *= beta;

		phi.spat_RealToSPHData(o_phi);
		u.spat_RealToSPHData(o_u);
		v.spat_RealToSPHData(o_v);
	}
};


#endif /* SRC_SWEREXISPH_HPP_ */
