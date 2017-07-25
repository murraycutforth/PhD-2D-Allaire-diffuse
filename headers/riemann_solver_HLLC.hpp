/*
 *	DESCRIPTION:	The HLLC approximate Riemann solver for the 
 * 			five-equation system of Allaire.
 * 
 * 	REFERENCE:	Garrick, Owkes, Regele - JCP - 2017
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		24/07/2017
 */


#ifndef HLLC_ALLAIRE_2D
#define HLLC_ALLAIRE_2D

#include "riemann_solver_base.hpp"
#include "mixturemodel.hpp"
#include "misc.hpp"
#include <cassert>
#include <cmath>
#include <algorithm>

class HLLC_riemann_solver : public riemann_solver_base {
	
	public:
	
	HLLC_riemann_solver (sim_info params, double gamma1, double gamma2, double pinf1, double pinf2)
	:
		riemann_solver_base (params, gamma1, gamma2, pinf1, pinf2)
	{}
	
	vectype solve_RP (const vectype& UL, const vectype& UR, double* u_star_ptr = nullptr, double* p_star_ptr = nullptr)
	{
		assert(UL.rows() == UR.rows());
		assert(UL.rows() == 6);
		
		double zL = UL(5);
		double rhoL = UL(0) + UL(1);
		double uL = UL(2) / rhoL;
		double vL = UL(3) / rhoL;
		double eL = UL(4) / rhoL - 0.5 * (uL * uL + vL * vL);
		double pL = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rhoL, eL, zL);
		double cL = allairemodel::mixture_soundspeed(gamma1, gamma2, pinf1, pinf2, rhoL, pL, zL);
		
		double zR = UR(5);
		double rhoR = UR(0) + UR(1);
		double uR = UR(2) / rhoR;
		double vR = UR(3) / rhoR;
		double eR = UR(4) / rhoR - 0.5 * (uR * uR + vR * vR);
		double pR = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rhoR, eR, zR);
		double cR = allairemodel::mixture_soundspeed(gamma1, gamma2, pinf1, pinf2, rhoR, pR, zR);
		
		double u_avg = 0.5 * (uL + uR);
		double c_avg = 0.5 * (cL + cR);
		
		double SL = std::min(u_avg - c_avg, uL - cL);
		double SR = std::max(u_avg + c_avg, uR + cR);
		
		double u_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR * (SR - uR));
		double p_star = pL + rhoL * (SL - uL) * (u_star - uL);
		
		double factorL = (SL - uL) / (SL - u_star);
		double factorR = (SR - uR) / (SR - u_star);
		
		
		// Use C-style pointer to return u_star and p_star as well if needed
		
		if (u_star_ptr) *u_star_ptr = u_star;
		if (p_star_ptr) *p_star_ptr = p_star;
		
		vectype U_starL (6);
		U_starL(0) = factorL * UL(0);
		U_starL(1) = factorL * UL(1);
		U_starL(2) = factorL * rhoL * u_star;
		U_starL(3) = factorL * UL(3);
		U_starL(4) = factorL * (UL(4) + (u_star - uL) * (rhoL * u_star + (pL / (SL - uL))));
		U_starL(5) = 0.0;
		
		vectype U_starR (6);
		U_starR(0) = factorR * UR(0);
		U_starR(1) = factorR * UR(1);
		U_starR(2) = factorR * rhoR * u_star;
		U_starR(3) = factorR * UR(3);
		U_starR(4) = factorR * (UR(4) + (u_star - uR) * (rhoR * u_star + (pR / (SR - uR))));
		U_starR(5) = 0.0;
		
		vectype fluxL = flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UL);
		vectype fluxR = flux_conserved_var(gamma1, gamma2, pinf1, pinf2, UR);
		
		return ((1.0 + std::copysign(1.0, u_star)) / 2.0)
			* (fluxL + std::min(0.0, SL) * (U_starL - UL))
			+ ((1.0 - std::copysign(1.0, u_star)) / 2.0)
			* (fluxR + std::max(0.0, SR) * (U_starR - UR));
	}
	
};

#endif
