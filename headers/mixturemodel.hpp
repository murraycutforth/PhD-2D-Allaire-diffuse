/*
 *	DESCRIPTION:	Functions which define the mixture pressure
 * 			and sound speed in the five equation Allaire
 * 			model.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */

#ifndef MIXTUREMODEL_ALLAIRE_H
#define MIXTUREMODEL_ALLAIRE_H

#include <cmath>

namespace allairemodel {

	inline double mixture_pressure (double gamma1, double gamma2, double pinf1, double pinf2, double rho, double e, double z)
	{
		double eps = z / (gamma1 - 1.0) + (1.0 - z) / (gamma2 - 1.0);
		double k = (z * gamma1 * pinf1 / (gamma1 - 1.0) + (1.0 - z) * gamma2 * pinf2 / (gamma2 - 1.0)) / eps;
		return rho * e / eps - k;
	}
	
	inline double mixture_soundspeed (double gamma1, double gamma2, double pinf1, double pinf2, double rho, double p, double z)
	{
		double eps = z / (gamma1 - 1.0) + (1.0 - z) / (gamma2 - 1.0);
		double k = (z * gamma1 * pinf1 / (gamma1 - 1.0) + (1.0 - z) * gamma2 * pinf2 / (gamma2 - 1.0)) / eps;
		return sqrt( (k + p * ((1.0 / eps) + 1.0) ) / rho);
	}

}

#endif