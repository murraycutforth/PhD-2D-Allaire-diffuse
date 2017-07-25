/*
 *	DESCRIPTION:	An abstract base class for the formulae used to
 *			update the volume fractions in the diffuse
 *			interface model.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		25/07/2017
 */


#ifndef ZUPDATE_H
#define ZUPDATE_H

#include <algorithm>

class zupdate_base {
	
	public:
	
	virtual double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R) =0;
};


class zupdate_upwind : public zupdate_base {
	
	public:
	
	double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R)
	{
		return z - (dt / dx) * ((zR - z) * std::min(0.0, u_star_R) + (z - zL) * std::max(0.0, u_star_L));
	}
};


class zupdate_secondorder : public zupdate_base {
	
	public:
	
	double zupdate (double dx, double dt, double zL, double z, double zR, double u_star_L, double u_star_R, double z_star_L, double z_star_R)
	{
		return z - (dt / dx) * (u_star_L * z_star_L - u_star_R * z_star_R - z * (u_star_R - u_star_L));
	}
};

#endif
