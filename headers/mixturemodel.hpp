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

#include "stiffened_gas_eos.hpp"
#include <cmath>

namespace allairemodel {

	inline double mixture_pressure (const binarySGparams& eosparams, double rho, double e, double z)
	{
		double eps = z / (eosparams.gamma1 - 1.0) + (1.0 - z) / (eosparams.gamma2 - 1.0);
		double k = (z * eosparams.gamma1 * eosparams.pinf1 / (eosparams.gamma1 - 1.0) 
				+ (1.0 - z) * eosparams.gamma2 * eosparams.pinf2 / (eosparams.gamma2 - 1.0)) / eps;
		return rho * e / eps - k;
	}
	
	inline double mixture_specific_ie (const binarySGparams& eosparams, double rho, double p, double z)
	{
		double eps = z / (eosparams.gamma1 - 1.0) + (1.0 - z) / (eosparams.gamma2 - 1.0);
		double k = (z * eosparams.gamma1 * eosparams.pinf1 / (eosparams.gamma1 - 1.0) 
				+ (1.0 - z) * eosparams.gamma2 * eosparams.pinf2 / (eosparams.gamma2 - 1.0)) / eps;
		return (p + k) * eps / rho;
	}
	
	inline double mixture_soundspeed (const binarySGparams& eosparams, double rho, double p, double z)
	{
		double eps = z / (eosparams.gamma1 - 1.0) + (1.0 - z) / (eosparams.gamma2 - 1.0);
		double k = (z * eosparams.gamma1 * eosparams.pinf1 / (eosparams.gamma1 - 1.0) 
				+ (1.0 - z) * eosparams.gamma2 * eosparams.pinf2 / (eosparams.gamma2 - 1.0)) / eps;
		return sqrt( (k + p * ((1.0 / eps) + 1.0) ) / rho);
	}
	
	inline double shock_speed_jump (const binarySGparams& eosparams, double p, double rho, double p_star, double z)
	{
		/*
		 * Such that shock speed S = u +/- shock_speed_jump()
		 */
		double eps = z / (eosparams.gamma1 - 1.0) + (1.0 - z) / (eosparams.gamma2 - 1.0);
		double k = (z * eosparams.gamma1 * eosparams.pinf1 / (eosparams.gamma1 - 1.0) 
				+ (1.0 - z) * eosparams.gamma2 * eosparams.pinf2 / (eosparams.gamma2 - 1.0)) / eps;
		double p_avg = 0.5 * (p_star + p);
		return sqrt((p_star + k + p_avg / eps) / rho);
	}
}

#endif
