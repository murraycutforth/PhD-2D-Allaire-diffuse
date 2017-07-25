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
	double e = U(4) / rho - (u * u + v * v);
	double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, U(5));
	
	flux(0) = U(0) * u;
	flux(1) = U(1) * u;
	flux(2) = U(2) * u + p;
	flux(3) = U(3) * u;
	flux(4) = (U(4) + p) * u;
	flux(5) = 0.0;
	
	return flux;
}


#endif
