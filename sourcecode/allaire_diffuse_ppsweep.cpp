#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>

void allaire_diffuse :: set_boundary_conditions (gridtype& grid, const sim_info& params)
{
	// TODO: Implement setting boundary conditions here
}

void allaire_diffuse :: pre_sweep (gridtype& grid, const sim_info& params)
{
	set_boundary_conditions(grid, params);
}
	
void allaire_diffuse :: post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params)
{
	for (int i=0; i<params.Ny; i++)
	{
		for (int j=0; j<params.Nx; j++)
		{
			grid[i][j] = future_grid[i][j];
		}
	}
}
