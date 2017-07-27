/*
 *	DESCRIPTION:	Miscellaneous useful functions
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef MISC_ALLAIRE_H
#define MISC_ALLAIRE_H

#include "typedefs.hpp"
#include "mixturemodel.hpp"

inline vectype flux_conserved_var (double gamma1, double gamma2, double pinf1, double pinf2, const vectype& U)
{
	/*
	 * The flux of conserved variables in the x-direction
	 */
	 
	vectype flux (6);
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, U(5));
	
	flux(0) = U(0) * u;
	flux(1) = U(1) * u;
	flux(2) = U(2) * u + p;
	flux(3) = U(3) * u;
	flux(4) = (U(4) + p) * u;
	flux(5) = 0.0;
	
	return flux;
}

inline bool is_physical_state (double gamma1, double gamma2, double pinf1, double pinf2, vectype& U)
{
	if (fabs(U(0)) < 1e-12) U(0) = 0.0;
	if (fabs(U(1)) < 1e-12) U(1) = 0.0;
	if (fabs(U(5)) < 1e-12) U(5) = 0.0;
	if (fabs(U(5) - 1.0) < 1e-12) U(5) = 1.0;
	
	
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double z = U(5);
	
	assert(U(0) >= 0.0 && U(1) >= 0.0 && e >= 0.0 && z >= 0.0 && z <= 1.0);
	return U(0) >= 0.0 && U(1) >= 0.0 && e >= 0.0 && z >= 0.0 && z <= 1.0;
}


#endif
