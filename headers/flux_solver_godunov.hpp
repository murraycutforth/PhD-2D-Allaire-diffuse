/*
 *	DESCRIPTION:	An implementation of Godunov's first order
 * 			flux solver applied to the five equation system
 * 			of Allaire.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		25/07/2017
 */


#ifndef GODUNOV_ALLAIRE_H
#define GODUNOV_ALLAIRE_H

#include "flux_solver_base.hpp"
#include <cassert>
#include <cmath>

class flux_solver_godunov : public flux_solver_base {
	
	public:
	
	flux_solver_godunov (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		flux_solver_base(RS_ptr, params, gamma1, gamma2, pinf1, pinf2)
	{}
	
	void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double& u_star, double& z_star)
	{
		assert(stencil.size() == 2);
				
		flux = RS_ptr->solve_RP(stencil[0], stencil[1], &u_star);
				
		z_star = ((1.0 + std::copysign(1.0, u_star)) / 2.0 ) * stencil[0](5)
			 + ((1.0 - std::copysign(1.0, u_star)) / 2.0 ) * stencil[1](5);
	}
};

#endif
