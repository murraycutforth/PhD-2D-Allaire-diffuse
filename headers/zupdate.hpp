/*
 *	DESCRIPTION:	An abstract base class for the formulae used to
 *			update the volume fractions in the diffuse
 *			interface model. The volume fractions obey:
 * 
 * 				Dz/Dt = 0
 * 
 * 			where D/Dt is the material derivative with velocity u:
 * 
 * 				D/Dt = partial / partial t + u . grad
 * 
 * 			The two discretisations are the first order
 * 			accurate upwind method, and a second order
 * 			method which uses star states obtained from the
 * 			MUSCL slope reconstruction of z and u.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		25/07/2017
 */


#ifndef ZUPDATE_H
#define ZUPDATE_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <iostream>

class zupdate_base {
	
public:
	
	virtual double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R) =0;
	
	virtual std::shared_ptr<zupdate_base> clone () =0;
};


class zupdate_upwind : public zupdate_base {
	
public:
	
	double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R)
	{
		double newz = z - (dt / dx) * ((zR - z) * std::min(0.0, u_star_R) + (z - zL) * std::max(0.0, u_star_L));
				
		return std::max(std::min(1.0, newz), 0.0);
	}
	
	std::shared_ptr<zupdate_base> clone ()
	{
		return std::make_shared<zupdate_upwind>();
	}
};


class zupdate_secondorder : public zupdate_base {
	
public:
	
	double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R)
	{
		double newz = z - (dt / dx) * (u_star_R * z_star_R - u_star_L * z_star_L - z * (u_star_R - u_star_L));
				
		return std::max(std::min(1.0, newz), 0.0);
	}
	
	std::shared_ptr<zupdate_base> clone ()
	{
		return std::make_shared<zupdate_secondorder>();
	}
};

#endif
