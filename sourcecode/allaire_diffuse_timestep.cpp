#include "allaire_diffuse.hpp"
#include "mixturemodel.hpp"
#include "stiffened_gas_eos.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>

double allaire_diffuse :: compute_dt (const gridtype& grid, const sim_info& params, int n, double t)
{
	double maxu = 0.0, maxv = 0.0;
	
	
	#pragma omp parallel for schedule(dynamic)
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0) + grid[i][j](1);
			double u = grid[i][j](2) / rho;
			double v = grid[i][j](3) / rho;
			double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
			double z = grid[i][j](5);
			double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, z);
			double c = allairemodel::mixture_soundspeed(gamma1, gamma2, pinf1, pinf2, rho, p, z);
			maxu = std::max(maxu, fabs(u) + c);
			maxv = std::max(maxv, fabs(v) + c);
		}
	}
	
	double CFL = params.CFL;
	
	if (n < 5) CFL = std::min(CFL, 0.2);
	
	double dt = CFL * std::min(params.dx / maxu, params.dy / maxv);
	
	if (t + dt > params.T) dt = params.T - t;
	
	return dt;	
}
