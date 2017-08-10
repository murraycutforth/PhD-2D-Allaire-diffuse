/*
 *	DESCRIPTION:	An abstract base class for solvers used to find
 * 			the flux used in the conservative update
 *			in the hyperbolic five equation system of Allaire.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef FS_ALLAIRE_H
#define FS_ALLAIRE_H

#include "typedefs.hpp"
#include "riemann_solver_base.hpp"
#include "stiffened_gas_eos.hpp"
#include <vector>
#include <memory>

class flux_solver_base {
	
	public:
	
	// Riemann solver
	
	std::shared_ptr<riemann_solver_base> RS_ptr;
	
	
	// Constants used in computation
	
	const sim_info params;
	binarySGparams eosparams;
	
	
	flux_solver_base (std::shared_ptr<riemann_solver_base> RS_ptr, sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		RS_ptr (RS_ptr),
		params (params),
		eosparams (gamma1, gamma2, pinf1, pinf2)
	{}
	
	virtual void flux_computation (const std::vector<vectype>& stencil, vectype& flux, double dt, double dx, double& u_star, double& z_star) =0;
	
	virtual std::shared_ptr<flux_solver_base> clone () =0;
};

#endif
