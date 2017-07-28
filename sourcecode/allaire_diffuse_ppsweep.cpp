#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include <cmath>
#include <memory>
#include <string>

void allaire_diffuse :: set_boundary_conditions (gridtype& grid, const sim_info& params)
{
	
	if (params.BC_L == "transmissive")
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=0; j<params.numGC; j++)
			{
				grid[i][j] = grid[i][2 * params.numGC - j - 1];
			}
		}
	}
	else if (params.BC_L == "reflective")
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=0; j<params.numGC; j++)
			{
				grid[i][j] = grid[i][2 * params.numGC - j - 1];
				grid[i][j](2) = - grid[i][j](2);
			}
		}
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid BC_L");
	}
	
	if (params.BC_R == "transmissive")
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=0; j<params.numGC; j++)
			{
				grid[i][j + params.Nx + params.numGC] = grid[i][params.Nx + params.numGC - 1 - j];
			}
		}
	}
	else if (params.BC_R == "reflective")
	{
		for (int i=params.numGC; i<params.Ny + params.numGC; i++)
		{
			for (int j=0; j<params.numGC; j++)
			{
				grid[i][j + params.Nx + params.numGC] = grid[i][params.Nx + params.numGC - 1 - j];
				grid[i][j + params.Nx + params.numGC](2) = - grid[i][j + params.Nx + params.numGC](2);
			}
		}
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid BC_R");
	}
	
	if (params.BC_B == "transmissive")
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			for (int i=0; i<params.numGC; i++)
			{
				grid[i][j] = grid[2 * params.numGC - i - 1][j];
			}
		}
	}
	else if (params.BC_B == "reflective")
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			for (int i=0; i<params.numGC; i++)
			{
				grid[i][j] = grid[2 * params.numGC - i - 1][j];
				grid[i][j](3) = - grid[i][j](3);
			}
		}
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid BC_B");
	}
	
	if (params.BC_T == "transmissive")
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			for (int i=0; i<params.numGC; i++)
			{
				grid[i + params.Ny + params.numGC][j] = grid[params.Ny + params.numGC - 1 - i][j];
			}
		}
	}
	else if (params.BC_T == "reflective")
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			for (int i=0; i<params.numGC; i++)
			{
				grid[i + params.Ny + params.numGC][j] = grid[params.Ny + params.numGC - 1 - i][j];
				grid[i + params.Ny + params.numGC][j](3) = - grid[i + params.Ny + params.numGC][j](3);
			}
		}
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid BC_T");
	}	
}

void allaire_diffuse :: pre_sweep (gridtype& grid, const sim_info& params)
{
	set_boundary_conditions(grid, params);
}
	
void allaire_diffuse :: post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params)
{	
	grid.swap(future_grid);
}
