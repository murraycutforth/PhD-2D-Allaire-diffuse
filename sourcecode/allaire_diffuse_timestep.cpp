#include "allaire_diffuse.hpp"
#include "mixturemodel.hpp"
#include "stiffened_gas_eos.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>

double allaire_diffuse :: compute_dt (const gridtype& grid, const sim_info& params, int n, double t)
{
	double maxu = 0.0, maxv = 0.0, mass1sum = 0.0, mass2sum = 0.0;
	
	
	#pragma omp parallel for schedule(dynamic) reduction(+:mass1sum, mass2sum)
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0) + grid[i][j](1);
			double u = grid[i][j](2) / rho;
			double v = grid[i][j](3) / rho;
			double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
			double z = grid[i][j](5);
			double p = allairemodel::mixture_pressure(eosparams, rho, e, z);
			double c = allairemodel::mixture_soundspeed(eosparams, rho, p, z);
			maxu = std::max(maxu, fabs(u) + c);
			maxv = std::max(maxv, fabs(v) + c);
			mass1sum += grid[i][j](0);
			mass2sum += grid[i][j](1);
		}
	}
	
	
	// Record fluid masses at this time step
	
	time.push_back(t);
	mass1.push_back(mass1sum * params.dx * params.dy);
	mass2.push_back(mass2sum * params.dx * params.dy);
	
	double CFL = params.CFL;
	
	if (n < 5) CFL = std::min(CFL, 0.2);
	
	double dt = CFL * std::min(params.dx / maxu, params.dy / maxv);
	
	if (t + dt > params.T) dt = params.T - t;
	
	return dt;	
}
