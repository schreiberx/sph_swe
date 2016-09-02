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

	/// timestep size
	double timestep_size;

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
			double i_coriolis_omega,
			double i_timestep_size
	)
	{
		timestep_size = i_timestep_size;

		alpha = i_alpha/timestep_size;
		beta = i_beta;

		const std::complex<double> &alpha = i_alpha;
		r = i_radius;
		coriolis_omega = i_coriolis_omega;
		two_omega = 2.0*i_coriolis_omega;

		sphConfig = i_sphConfig;

		sphSolverPhi.setup(sphConfig, 4);
		sphSolverPhi.solver_component_rexi_z1(	(alpha*alpha)*(alpha*alpha), r);
		sphSolverPhi.solver_component_rexi_z2(	2.0*two_omega*two_omega*alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z3(	(two_omega*two_omega)*(two_omega*two_omega), r);
		sphSolverPhi.solver_component_rexi_z4(	-alpha*two_omega, r);
		sphSolverPhi.solver_component_rexi_z5(	1.0/alpha*two_omega*two_omega*two_omega, r);
		sphSolverPhi.solver_component_rexi_z6(	2.0*two_omega*two_omega, r);
		sphSolverPhi.solver_component_rexi_z7(	-alpha*alpha, r);
		sphSolverPhi.solver_component_rexi_z8(	-two_omega*two_omega, r);

		sphSolverVel.setup(sphConfig, 2);
		sphSolverVel.solver_component_rexi_z2(	two_omega*two_omega, r);
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

		SPHDataComplex div0(op.div(u0, v0));
		SPHDataComplex eta0(op.vort(u0, v0));


		SPHDataComplex Fck = two_omega*op.grad_lon(mu)*(-(alpha*alpha*u0 - two_omega*two_omega*op.mu2(u0)) + 2.0*alpha*two_omega*op.mu(v0));
		SPHDataComplex foo = div0 - two_omega*(1.0/alpha)*op.mu(eta0) + alpha*i_phi0 + two_omega*two_omega*(1.0/alpha)*op.mu2(i_phi0);
		SPHDataComplex rhs = alpha*alpha*foo + two_omega*two_omega*op.mu2(foo) + (1.0/alpha)*Fck;

		sphSolverPhi.solve(rhs);
		SPHDataComplex &phi = rhs;

		SPHDataComplex a = u0 + op.grad_lon(phi);
		SPHDataComplex b = v0 + op.grad_lat(phi);

		SPHDataComplex rhsa = alpha*a - op.mu(b);
		SPHDataComplex rhsb = op.mu(a) + alpha*b;

		sphSolverVel.solve(rhsa);
		SPHDataComplex &u = rhsa;

		sphSolverVel.solve(rhsb);
		SPHDataComplex &v = rhsb;

		phi.spat_RealToSPHData(o_phi);
		u.spat_RealToSPHData(o_u);
		v.spat_RealToSPHData(o_v);

		o_phi = o_phi*timestep_size;
		o_u = o_u*timestep_size;
		o_v = o_v*timestep_size;
	}
};


#endif /* SRC_SWEREXISPH_HPP_ */
