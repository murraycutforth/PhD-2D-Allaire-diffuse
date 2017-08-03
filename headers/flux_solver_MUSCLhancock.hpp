/*
 *	DESCRIPTION:	An implementation of the second order accurate
 *			MUSCL-Hancock flux solver applied to the 
 *			five equation system of Allaire. The family
 * 			of methods described as MUSCL-Hancock is quite
 * 			wide, and a couple of different variations are
 * 			implemented here for comparison purposes.
 * 
 * 			1)	- Conserved variable reconstruction
 * 				- Monotone - centred slope limiting
 * 				- Boundary flux time evolution with a
 * 				  separate half step evolution formula
 * 				  for volume fractions
 * 
 * 			2)	- Primitive variable reconstruction
 * 				- Minmod slope limiting
 * 				- Jacobian-based half step evolution 
 *  
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		25/07/2017
 */


#ifndef MUSCLHANCOCK_ALLAIRE_H
#define MUSCLHANCOCK_ALLAIRE_H

#include "flux_solver_base.hpp"
#include "riemann_solver_PVRS.hpp"
#include "misc.hpp"
#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>

vectype unlimited_central_diff (const vectype& diff_L, const vectype& diff_R)
{
	/*
	 * Central difference with no limiting. Reproduces Fromm's method.
	 */
	 
	vectype slope = 0.5 * (diff_L + diff_R);
	return slope;
}

vectype MC_limited_slope (const vectype& diff_L, const vectype& diff_R)
{
	/*
	 * Monotone central limited slope. Less diffusive than minmod.
	 */
	 
	vectype slope (6);
	double sigmai;
	
	for (int k=0; k<6; k++)
	{
		if (diff_L(k) > 0.0 && diff_R(k) > 0.0) sigmai = 1.0;
		else if (diff_L(k) < 0.0 && diff_R(k) < 0.0) sigmai = - 1.0;
		else sigmai = 0.0;
		
		slope(k) = sigmai * std::min(0.5 * fabs(diff_L(k) + diff_R(k)), 
						std::min(2.0 * fabs(diff_L(k)), 2.0 * fabs(diff_R(k))));
	}
	
	return slope;
}

vectype minmod_limited_slope (const vectype& diff_L, const vectype& diff_R)
{
	/*
	 * Minmod slope limiter.
	 */
	 
	vectype slope (6);
	
	for (int k=0; k<6; k++)
	{
		if (diff_R(k) > 0.0)
		{
			slope(k) = std::max(0.0, std::min(diff_L(k), diff_R(k)));
		}
		else
		{
			slope(k) = std::min(0.0, std::max(diff_L(k), diff_R(k)));
		}
	}
	
	return slope;
}

class flux_solver_MUSCLHANCOCK_1 : public flux_solver_base {
	
	public:
	
	vectype diff_L;
	vectype diff_C;
	vectype diff_R;
	vectype UL_L;
	vectype UL_R;
	vectype UR_L;
	vectype UR_R;
	vectype UL;
	vectype UR;
	
	flux_solver_MUSCLHANCOCK_1 (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		flux_solver_base(RS_ptr, params, gamma1, gamma2, pinf1, pinf2),
		diff_L	(6),
		diff_C	(6),
		diff_R	(6),
		UL_L	(6),
		UL_R	(6),
		UR_L	(6),
		UR_R	(6),
		UL	(6),
		UR	(6)
	{}
		
	
	void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double& z_star)
	{
		assert(stencil.size() == 4);
				
		diff_L = stencil[1] - stencil[0];
		diff_C = stencil[2] - stencil[1];
		diff_R = stencil[3] - stencil[2];
		
		UL_L = stencil[1] - 0.5 * minmod_limited_slope(diff_L, diff_C);
		UL_R = stencil[1] + 0.5 * minmod_limited_slope(diff_L, diff_C);
		UR_L = stencil[2] - 0.5 * minmod_limited_slope(diff_C, diff_R);
		UR_R = stencil[2] + 0.5 * minmod_limited_slope(diff_C, diff_R);
		
		
		// Use zero slope if reconstructed pressure or density goes negative
		
		if ( is_physical_state(gamma1, gamma2, pinf1, pinf2, UL_L)
			&& is_physical_state(gamma1, gamma2, pinf1, pinf2, UL_R)
			&& is_physical_state(gamma1, gamma2, pinf1, pinf2, UR_L)
			&& is_physical_state(gamma1, gamma2, pinf1, pinf2, UR_R)
		   )
		{}
		else
		{
			UL_L = stencil[1];
			UL_R = stencil[1];
			UR_L = stencil[2];
			UR_R = stencil[2];
			
			std::cout << "[MUSCLHANCOCK1] Unphysical state created by slope extrapolation." << std::endl;
		}
		
		
		// Advance conserved variables by half time step
		
		UL = UL_R + 0.5 * (dt / dx) * (flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UL_L) - flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UL_R));
		UR = UR_L + 0.5 * (dt / dx) * (flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UR_L) - flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UR_R));
		
		
		// Advance volume fractions by half time step
		
		double uL_L = UL_L(2) / (UL_L(0) + UL_L(1));
		double uL_R = UL_R(2) / (UL_R(0) + UL_R(1));
		double uR_L = UR_L(2) / (UR_L(0) + UR_L(1));
		double uR_R = UR_R(2) / (UR_R(0) + UR_R(1));
		UL(5) = UL_R(5) - 0.5 * (dt / dx) * (uL_R * UL_R(5) - uL_L * UL_L(5) - UL_R(5) * (uL_R - uL_L));
		UR(5) = UR_L(5) - 0.5 * (dt / dx) * (uR_R * UR_R(5) - uR_L * UR_L(5) - UR_R(5) * (uR_R - uR_L));
		
		
		// Use zero slope if reconstructed pressure or density goes negative
		
		if ( is_physical_state(gamma1, gamma2, pinf1, pinf2, UL)
			&& is_physical_state(gamma1, gamma2, pinf1, pinf2, UR)
		   )
		{}
		else
		{
			
			
			double rhoL = UL(0) + UL(1);
			double uL = UL(2) / rhoL;
			double vL = UL(3) / rhoL;
			double eL = UL(4) / rhoL - 0.5 * (uL * uL + vL * vL);
			double zL = UL(5);
			double rhoR = UR(0) + UR(1);
			double uR = UR(2) / rhoR;
			double vR = UR(3) / rhoR;
			double eR = UR(4) / rhoR - 0.5 * (uR * uR + vR * vR);
			double zR = UR(5);
			
			if (rhoL < 0 || rhoR < 0 || UL(0) < 0 || UL(1) < 0 || UR(0) < 0 || UR(1) < 0) std::cout << "bad density" << std::endl;
			if (zL < 0 || zL > 1 || zR < 0 || zL > 1) std::cout << "bad volume fraction" << std::endl;
			if (eL < 0 || eR < 0) std::cout << "bad internal energy" << std::endl;
			
			std::cout << "[MUSCLHANCOCK1] Unphysical state created by time evolution." << std::endl;
			UL = stencil[1];
			UR = stencil[2];
		}
		
		flux = RS_ptr->solve_RP(UL, UR, &u_star);
		
		if (u_star > 0.0)
		{
			z_star = UL(5);
		}
		else
		{
			z_star = UR(5);
		}
	}
	
	std::shared_ptr<flux_solver_base> clone ()
	{
		return std::make_shared<flux_solver_MUSCLHANCOCK_1>(RS_ptr->clone(), params, gamma1, gamma2, pinf1, pinf2);
	}
};

class flux_solver_MUSCLHANCOCK_2 : public flux_solver_base {
	
	public:
	
	vectype diff_L;
	vectype diff_C;
	vectype diff_R;
	vectype del_L;
	vectype del_R;
	vectype prim_L_n;
	vectype prim_R_n;
	Matrix6d A;
	vectype UL;
	vectype UR;
	PVRS_riemann_solver pvrs;
	std::function<vectype(const vectype&, const vectype&)> slopelimiter;
	
	flux_solver_MUSCLHANCOCK_2 (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		flux_solver_base(RS_ptr, params, gamma1, gamma2, pinf1, pinf2),
		diff_L	(6),
		diff_C	(6),
		diff_R	(6),
		del_L	(6),
		del_R	(6),
		A	(),
		UL	(6),
		UR	(6),
		pvrs	(params, gamma1, gamma2, pinf1, pinf2),
		slopelimiter (minmod_limited_slope)
	{}
		
	
	void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double& z_star)
	{
		assert(stencil.size() == 4);
		prim_L_n = conserved_to_primitives(gamma1, gamma2, pinf1, pinf2, stencil[1]);
		prim_R_n = conserved_to_primitives(gamma1, gamma2, pinf1, pinf2, stencil[2]);
		
		// Compute limited slopes of primitive variables using slope limiter
		
		/*
				
		diff_L = prim_L_n - conserved_to_primitives(gamma1, gamma2, pinf1, pinf2, stencil[0]);
		diff_C = prim_R_n - prim_L_n;
		diff_R = conserved_to_primitives(gamma1, gamma2, pinf1, pinf2, stencil[3]) - prim_R_n;
		
		del_L = slopelimiter(diff_L, diff_C);
		del_R = slopelimiter(diff_C, diff_R);
		
		*/
		
		
		// OR compute slopes using characteristic limiting
		
		
		
		std::vector<vectype> jumpsL (3, vectype(6));
		std::vector<vectype> jumpsR (3, vectype(6));
		
		pvrs.compute_primitive_jumps(stencil[0], stencil[1], jumpsL);
		pvrs.compute_primitive_jumps(stencil[1], stencil[2], jumpsR);
		del_L = Eigen::VectorXd::Zero(6);
		
		for (int k=0; k<3; k++)
		{
			del_L += slopelimiter(jumpsL[k], jumpsR[k]);
		}
		
		pvrs.compute_primitive_jumps(stencil[1], stencil[2], jumpsL);
		pvrs.compute_primitive_jumps(stencil[2], stencil[3], jumpsR);
		del_R = Eigen::VectorXd::Zero(6);
		
		for (int k=0; k<3; k++)
		{
			del_R += slopelimiter(jumpsL[k], jumpsR[k]);
		}
		
		
		
		
		
		// Compute left/right state to be supplied to Riemann solver
		
		A_primitive_vars(gamma1, gamma2, pinf1, pinf2, prim_L_n, A);
		
		UL = prim_L_n + 0.5 * del_L - 0.5 * (dt / dx) * A * del_L;
		
		A_primitive_vars(gamma1, gamma2, pinf1, pinf2, prim_R_n, A);
		
		UR = prim_R_n - 0.5 * del_R - 0.5 * (dt / dx) * A * del_R;
		
		
		// Convert back to conserved variables
		
		UL = primitives_to_conserved(gamma1, gamma2, pinf1, pinf2, UL);
		UR = primitives_to_conserved(gamma1, gamma2, pinf1, pinf2, UR);

		
		// Use zero slope if reconstructed pressure or density goes negative
		
		if ( is_physical_state(gamma1, gamma2, pinf1, pinf2, UL)
			&& is_physical_state(gamma1, gamma2, pinf1, pinf2, UR)
		   )
		{}
		else
		{
			
			
			double rhoL = UL(0) + UL(1);
			double uL = UL(2) / rhoL;
			double vL = UL(3) / rhoL;
			double eL = UL(4) / rhoL - 0.5 * (uL * uL + vL * vL);
			double zL = UL(5);
			double rhoR = UR(0) + UR(1);
			double uR = UR(2) / rhoR;
			double vR = UR(3) / rhoR;
			double eR = UR(4) / rhoR - 0.5 * (uR * uR + vR * vR);
			double zR = UR(5);
			
			if (rhoL < 0 || rhoR < 0 || UL(0) < 0 || UL(1) < 0 || UR(0) < 0 || UR(1) < 0) std::cout << "bad density" << std::endl;
			if (zL < 0 || zL > 1 || zR < 0 || zL > 1) std::cout << "bad volume fraction" << std::endl;
			if (eL < 0 || eR < 0) std::cout << "bad internal energy" << std::endl;
			
			std::cout << "[MUSCLHANCOCK2] Unphysical state created by time evolution." << std::endl;
			UL = primitives_to_conserved(gamma1, gamma2, pinf1, pinf2, stencil[1]);
			UR = primitives_to_conserved(gamma1, gamma2, pinf1, pinf2, stencil[2]);
		}
		
		flux = RS_ptr->solve_RP(UL, UR, &u_star);
		z_star = 0.5 * (1.0 + std::copysign(1.0, u_star)) * UL(5) + 0.5 * (1.0 - std::copysign(1.0, u_star)) * UR(5);
	}
	
	std::shared_ptr<flux_solver_base> clone ()
	{
		return std::make_shared<flux_solver_MUSCLHANCOCK_2>(RS_ptr->clone(), params, gamma1, gamma2, pinf1, pinf2);
	}
};

#endif
