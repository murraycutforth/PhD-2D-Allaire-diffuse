/*
 *	DESCRIPTION:	Miscellaneous useful functions for the Allaire system
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef MISC_ALLAIRE_H
#define MISC_ALLAIRE_H

#include "typedefs.hpp"
#include "mixturemodel.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cassert>

inline double specific_ie (const vectype& U)
{
	/*
	 * The specific internal energy from vector of conserved vars
	 */
	
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	return (U(4) / rho) - 0.5 * (u*u + v*v);
}

inline vectype flux_conserved_var (const binarySGparams& eosparams, const vectype& U)
{
	/*
	 * The flux of conserved variables in the x-direction
	 */
	 
	vectype flux (6);
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double p = allairemodel::mixture_pressure(eosparams, rho, e, U(5));
	
	flux(0) = U(0) * u;
	flux(1) = U(1) * u;
	flux(2) = U(2) * u + p;
	flux(3) = U(3) * u;
	flux(4) = (U(4) + p) * u;
	flux(5) = 0.0;
	
	return flux;
}

inline vectype flux_primitive_var (const binarySGparams& eosparams, const vectype& W)
{
	vectype flux (6);
	double rho = W(0) + W(1);
	double e = allairemodel::mixture_specific_ie(eosparams, rho, W(4), W(5));
	double E = rho * e + 0.5 * rho * (W(2) * W(2) + W(3) * W(3));
	
	flux(0) = W(0) * W(2);
	flux(1) = W(1) * W(2);
	flux(2) = rho * W(2) * W(2) + W(4);
	flux(3) = rho * W(3) * W(2);
	flux(4) = (E + W(4)) * W(2);
	flux(5) = 0.0;
	
	return flux;
}

inline vectype conserved_to_primitives (const binarySGparams& eosparams, const vectype& U)
{
	vectype prims (6);
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double p = allairemodel::mixture_pressure(eosparams, rho, e, U(5));
	
	prims(0) = U(0);
	prims(1) = U(1);
	prims(2) = u;
	prims(3) = v;
	prims(4) = p;
	prims(5) = U(5);
	
	return prims;
}

inline vectype primitives_to_conserved (const binarySGparams& eosparams, const vectype& W)
{
	vectype conserved (6);
	double rho = W(0) + W(1);
	double e = allairemodel::mixture_specific_ie(eosparams, rho, W(4), W(5));
	
	conserved(0) = W(0);
	conserved(1) = W(1);
	conserved(2) = rho * W(2);
	conserved(3) = rho * W(3);
	conserved(4) = rho * e + 0.5 * rho * (W(2) * W(2) + W(3) * W(3));
	conserved(5) = W(5);
	
	return conserved;
}

inline vectype primitives_to_conserved (const binarySGparams& eosparams, double rhoz1, double rhoz2, double u, double v, double p, double z)
{
	vectype W (6);
	W << rhoz1, rhoz2, u, v, p, z;
	return primitives_to_conserved(eosparams, W);
}

inline void A_primitive_vars (const binarySGparams& eosparams, const vectype& W, Matrix6d& A)
{
	/*
	 * Use vector of primitive variables (W) to set the value of the
	 * matrix A which is the Jacobian of the system in primitive
	 * variable form.
	 */
	 
	double rho = W(0) + W(1);
	double c = allairemodel::mixture_soundspeed(eosparams, rho, W(4), W(5));
	
	A = W(2) * Eigen::Matrix<double, 6, 6>::Identity();
	A(0,2) = W(0);
	A(1,2) = W(1);
	A(2,4) = 1.0 / rho;
	A(4,2) = rho * c * c;
}

inline bool is_physical_state (const binarySGparams& eosparams, const vectype& U)
{
	double rho = U(0) + U(1);
	double u = U(2) / rho;
	double v = U(3) / rho;
	double e = U(4) / rho - 0.5 * (u * u + v * v);
	double z = U(5);
		
	return U(0) >= 0.0 && U(1) >= 0.0 && e >= 0.0 && z >= 0.0 && z <= 1.0;
}


inline bool is_physical_state (const binarySGparams& eosparams, vectype& U)
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
		
	return U(0) >= 0.0 && U(1) >= 0.0 && e >= 0.0 && z >= 0.0 && z <= 1.0;
}


#endif
