#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include "misc.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>


void allaire_diffuse :: update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t)
{
	// Storage for flux and velocity across each internal cell edge
	
	rowtype fluxes (params.Nx + 2 * params.numGC, vectype(6));
	std::vector<double> u_stars (params.Nx + 2 * params.numGC, 0.0);
	std::vector<double> z_stars (params.Nx + 2 * params.numGC, 0.0);
		
	
	// Compute cell edge fluxes across each edge of this row
	
	std::vector<vectype> stencil (2*params.stclsize);
		
	for (int j=params.numGC; j < params.Nx + params.numGC + 1; j++)
	{
		for (int l = j - params.stclsize; l <= j + params.stclsize - 1; l++)
		{
			stencil[l - j + params.stclsize] = grid[i][l];
		}

		FS_ptr->flux_computation(stencil, fluxes[j], dt, params.dx, u_stars[j], z_stars[j]);
	}
		
	
	// Conservative update formula
	
	for (int j=params.numGC; j < params.Nx + params.numGC; j++)
	{
		future_grid[i][j] = grid[i][j] + (dt / params.dx) * (fluxes[j] - fluxes[j+1]);
	}
	
	
	// Volume fraction update formula
	
	for (int j=params.numGC; j < params.Nx + params.numGC; j++)
	{
		future_grid[i][j](5) = zupdate_ptr->zupdate(params.dx, dt, grid[i][j-1](5), grid[i][j](5), grid[i][j+1](5), 
							    u_stars[j], u_stars[j+1], z_stars[j], z_stars[j+1]);
							    
		assert(is_physical_state(eosparams, future_grid[i][j]));
	}
	
}
