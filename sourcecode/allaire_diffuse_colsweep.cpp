#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include "misc.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>

void allaire_diffuse :: update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t)
{
	// Storage for flux and velocity across each internal cell edge
	
	rowtype fluxes (params.Ny + 2 * params.numGC, vectype(6));
	rowtype transformed_row (params.Ny + 2 * params.numGC, vectype(6));
	std::vector<double> u_stars (params.Ny + 2 * params.numGC, 0.0);
	std::vector<double> z_stars (params.Ny + 2 * params.numGC, 0.0);
	
	
	// Rotate problem by -pi/2 so that x-direction methods may be applied
	
	double temp;
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		transformed_row[i] = grid[i][j];
		temp = transformed_row[i](2);
		transformed_row[i](2) = transformed_row[i](3);
		transformed_row[i](3) = - temp;
	}
	
	
	// Compute fluxes across each edge of transformed row
	
	std::vector<vectype> stencil (2*params.stclsize);
	
	for (int i=params.numGC; i<params.Ny + params.numGC + 1; i++)
	{
		for (int l = i - params.stclsize; l <= i + params.stclsize - 1; l++)
		{
			stencil[l - i + params.stclsize] = transformed_row[l];
		}
		
		FS_ptr->flux_computation(stencil, fluxes[i], dt, params.dx, u_stars[i], z_stars[i]);
	}
	
	
	// Conservative update formula
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		future_grid[i][j] = transformed_row[i] + (dt / params.dy) * (fluxes[i] - fluxes[i+1]);
	}
	
	
	// Volume fraction update formula
		
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		future_grid[i][j](5) = zupdate_ptr->zupdate(params.dy, dt, transformed_row[i-1](5), transformed_row[i](5), transformed_row[i+1](5), 
							    u_stars[i], u_stars[i+1], z_stars[i], z_stars[i+1]);
		
		assert(is_physical_state(eosparams, future_grid[i][j]));
	}
	
	
	// Rotate this column of future_grid by pi/2 back into original frame
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		temp = future_grid[i][j](2);
		future_grid[i][j](2) = - future_grid[i][j](3);
		future_grid[i][j](3) = temp;
	}
}
