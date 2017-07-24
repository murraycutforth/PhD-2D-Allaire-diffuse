#include "allaire_diffuse.hpp"
#include "stiffened_gas_eos.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>

std::shared_ptr<gridtype> allaire_diffuse :: set_ICs (settings_file SF, sim_info& params)
{
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.CFL = SF.CFL;
	params.outputname = SF.basename;
	
	if (SF.riemann_solver == "HLLC")
	{
		// TODO: pointer to abstract Riemann solver class
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid riemann_solver in settings file.");
	}
	
	if (SF.flux_solver == "Godunov")
	{
		params.numGC = 1;
		// TODO: pointer to abstract flux solver class?
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid flux_solver in settings file.");
	}
	
	if (SF.volfracupdater == "Upwind")
	{
		// TODO: pointer to abstract volume fraction update class?
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid volfracupdater in settings file.");
	}
	
		
	
	gridtype ICgrid (SF.Ny + 2 * params.numGC, rowtype(SF.Nx + 2 * params.numGC, vectype(6)));
	
	if (SF.test_case == "ST1")
	{

		// EOS parameters:
		
		gamma1 = 1.4;
		gamma2 = 2.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		
		// Computational domain and boundary conditions:
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.14;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
		
		
		// Initial fluid states:
		
		double rho1 = 1.0;
		double p1 = 1.0;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 0.125;
		double p2 = 0.1;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u = 0.0;
		double v = 0.0;
		double z;
		
		for (int i=0; i<SF.Ny; i++)
		{
			for (int j=0; j<SF.Nx; j++)
			{
				if (j < SF.Nx / 2)
				{
					// Fluid 1 here
					z = 1.0;
				}
				else
				{
					// Fluid 2 here
					z = 0.0;
				}
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid test_case in settings file.");
	}
		
	return std::make_shared<gridtype>(ICgrid);
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
			std::cout << "(" << grid[i][j](0) << "," << grid[i][j](1) << "," << grid[i][j](2) << "," << grid[i][j](3) << "," << grid[i][j](4) << "," << grid[i][j](5) << ")" << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
