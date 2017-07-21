#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include <iostream>
#include <cmath>
#include <memory>

std::shared_ptr<gridtype> allaire_diffuse :: set_ICs (settings_file SF, sim_info& params)
{
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.T = 0.9;
	
	// Initialise a 2D grid of vectors in R^5
	gridtype testgrid (SF.Ny, rowtype(SF.Nx, vectype(3)));
	
	// Loop over grid and set element i,j to i*Nx + j except for last component
	for (int i=0; i<SF.Ny; i++)
	{
		for (int j=0; j<SF.Nx; j++)
		{
			testgrid[i][j](0) = i;
			testgrid[i][j](1) = j;
			testgrid[i][j](2) = 0.0;
		}
	}
		
	return std::make_shared<gridtype>(testgrid);
}
	
void allaire_diffuse :: output (const gridtype& grid, const sim_info& params, int n, double t)
{
	// Print info to screen for verification
	
	std::cout << std::endl << "[allaire_diffuse] Outputting..." << std::endl;
	std::cout << "[allaire_diffuse] Current state of grid is: " << std::endl;
	
	for (int i=0; i<params.Ny; i++)
	{
		for (int j=0; j<params.Nx; j++)
		{
			std::cout << "(" << grid[i][j](0) << "," << grid[i][j](1) << "," << grid[i][j](2) << ")" << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
	
double allaire_diffuse :: compute_dt (const gridtype& grid, const sim_info& params, double t)
{
	return 0.5;
}
	
void allaire_diffuse :: pre_sweep (gridtype& grid, const sim_info& params)
{
	// This is where we set boundary conditions
	
	for (int i=0; i<params.Ny; i++)
	{
		grid[i][params.Nx-1](2) = 7.0;
	}
}
	
void allaire_diffuse :: update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t)
{
	// Iterate along row and do very expensive operation on each element
	
	int N = 1000000.0;
		
	for (int j=1; j<params.Nx - 1; j++)
	{
		double sum = 0.0;
		
		for (int k=0; k<N; k++)
		{
			sum += sin(k) + pow(cos(k)*cos(k), 1.0 / sqrt(k+1));
		}
		
		future_grid[i][j](2) = sum / N + grid[i][j - 1](2) + grid[i][j + 1](2);
	}
}
	
void allaire_diffuse :: update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t)
{
	// Iterate down column and do very expensive operation on each element
		
	int N = 1000000.0;
		
	for (int i=1; i<params.Ny - 1; i++)
	{
		double sum = 0.0;
		
		for (int k=0; k<N; k++)
		{
			sum += sin(k) + pow(cos(k)*cos(k), 1.0 / sqrt(k+1));
		}
		
		future_grid[i][j](2) = sum / N + grid[i - 1][j](2) + grid[i + 1][j](2);
	}
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
