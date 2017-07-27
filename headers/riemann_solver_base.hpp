/*
 *	DESCRIPTION:	An abstract base class for Riemann solvers used
 *			in the hyperbolic five equation system of Allaire.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef RS_ALLAIRE_2D
#define RS_ALLAIRE_2D

#include "sim_info.hpp"
#include "typedefs.hpp"
#include <memory>

class riemann_solver_base {
	
	public:
	
	// Constants used in computation
	
	const sim_info params;
	const double gamma1;
	const double gamma2;
	const double pinf1;
	const double pinf2;
	
	
	riemann_solver_base (sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		params (params),
		gamma1 (gamma1),
		gamma2 (gamma2),
		pinf1 (pinf1),
		pinf2 (pinf2)
	{}
	
	virtual vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr) =0;
	
	virtual std::shared_ptr<riemann_solver_base> clone () =0;
	
};	

#endif
