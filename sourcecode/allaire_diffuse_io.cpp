#include "allaire_diffuse.hpp"
#include "stiffened_gas_eos.hpp"
#include "flux_solver_godunov.hpp"
#include "mixturemodel.hpp"
#include "riemann_solver_HLLC.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <fstream>
#include <string>

std::shared_ptr<gridtype> allaire_diffuse :: set_ICs (settings_file SF, sim_info& params)
{
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.CFL = SF.CFL;
	params.outputname = SF.basename;
	
	std::shared_ptr<riemann_solver_base> RS_ptr = nullptr;
	
	
	// EOS parameters, computational domain and boundary conditions
	
	if (SF.test_case == "ST1")
	{
		gamma1 = 1.4;
		gamma2 = 2.4;
		pinf1 = 0.0;
		pinf2 = 0.0;

		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.14;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid test_case in settings file.");
	}
	
	
	// Number of ghost cells needed
	
	if (SF.flux_solver == "Godunov")
	{
		params.numGC = 1;
		params.stclsize = 1;
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid flux_solver in settings file.");
	}
	
	
	// Riemann solver object
	
	if (SF.riemann_solver == "HLLC")
	{
		RS_ptr = std::make_shared<HLLC_riemann_solver>(params, gamma1, gamma2, pinf1, pinf2);
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid riemann_solver in settings file.");
	}
	
	
	// Flux solver object
	
	if (SF.flux_solver == "Godunov")
	{
		FS_ptr = std::make_shared<flux_solver_godunov>(RS_ptr, params, gamma1, gamma2, pinf1, pinf2);
		zupdate_ptr = std::make_shared<zupdate_upwind>();
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid flux_solver in settings file.");
	}
			
	
	gridtype ICgrid (SF.Ny + 2 * params.numGC, rowtype(SF.Nx + 2 * params.numGC, vectype(6)));
	
	
	// Initial states
	
	if (SF.test_case == "ST1")
	{
		double rho1 = 1.0;
		double p1 = 1.0;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 0.125;
		double p2 = 0.1;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u = 0.0;
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				if (j < params.Nx / 2)
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

void allaire_diffuse :: vtk_output (const gridtype& grid, const sim_info& params, int n, double t)
{
	std::string filename5 = params.outputname + "-state-" + std::to_string(n) + ".vtk";
	std::ofstream outfile5;
	outfile5.open(filename5);
	outfile5 << "# vtk DataFile Version 3.0\n";
	outfile5 << "Fluid state\n";
	outfile5 << "ASCII\n";
	outfile5 << "DATASET RECTILINEAR_GRID\n";

	outfile5 << "FIELD FieldData 2\n";
	outfile5 << "TIME 1 1 double\n";
	outfile5 << t << std::endl;
	outfile5 << "CYCLE 1 1 int\n";
	outfile5 << n << std::endl;

	outfile5 << "DIMENSIONS " + std::to_string(params.Nx+1) + " " + std::to_string(params.Ny+1) + " 1\n";
	outfile5 << "X_COORDINATES " + std::to_string(params.Nx+1) + " double\n";
	for (int i=0; i<params.Nx+1; i++)
	{
		outfile5 << params.x0 + double(i)*params.dx << " ";
	}
	outfile5 << std::endl;
	outfile5 << "Y_COORDINATES " + std::to_string(params.Ny+1) + " double\n";
	for (int j=0; j<params.Ny+1; j++)
	{
		outfile5 << params.y0 + double(j)*params.dy << " ";
	}
	outfile5 << std::endl;
	outfile5 << "Z_COORDINATES 1 double\n";
	outfile5 << "0\n";
	outfile5 << "CELL_DATA " + std::to_string(params.Nx*params.Ny) + "\n";

	outfile5 << "VECTORS vfield double\n";
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0) + grid[i][j](1);
			double u = grid[i][j](2) / rho;
			double v = grid[i][j](3) / rho;
			outfile5 << u << " " << v << " " << 0.0  << "\n";
		}
	}

	outfile5 << "FIELD FieldData 1\n";
	outfile5 << "Density 1 " + std::to_string(params.Nx*params.Ny) + " double\n";
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0) + grid[i][j](1);
			outfile5 << rho << "\n";
		}
	}
	
	outfile5 << "FIELD FieldData 1\n";
	outfile5 << "Pressure 1 " + std::to_string(params.Nx*params.Ny) + " double\n";
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
			outfile5 << p << "\n";
		}
	}
	
	outfile5 << "FIELD FieldData 1\n";
	outfile5 << "Specificinternalenergy 1 " + std::to_string(params.Nx*params.Ny) + " double\n";
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho = grid[i][j](0) + grid[i][j](1);
			double u = grid[i][j](2) / rho;
			double v = grid[i][j](3) / rho;
			double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
			outfile5 << e << "\n";
		}
	}
	
	outfile5 << "FIELD FieldData 1\n";
	outfile5 << "Volumefraction 1 " + std::to_string(params.Nx*params.Ny) + " double\n";
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			outfile5 << grid[i][j](5) << "\n";
		}
	}
}

void allaire_diffuse :: gnuplot_lineout (const gridtype& grid, const sim_info& params, int n, double t)
{
	std::string filename = params.outputname + "-lineoutx-" + std::to_string(n) + ".dat";
	std::ofstream outfile;
	outfile.open(filename);
	
	int i = params.Ny / 2;
	
	for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	{
		double x = params.x0 + double(j - params.numGC)*params.dx;
		double rho = grid[i][j](0) + grid[i][j](1);
		double u = grid[i][j](2) / rho;
		double v = grid[i][j](3) / rho;
		double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
		double z = grid[i][j](5);
		double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, z);
			
		outfile << x << " " << rho << " " << u << " " << v << " " << e << " " << p << " " << z << std::endl;
	}
	
	outfile.close();
}
	
void allaire_diffuse :: output (const gridtype& grid, const sim_info& params, int n, double t)
{
	if (n==0 || t==params.T)
	{
		vtk_output(grid, params, n, t);
		gnuplot_lineout(grid, params, n, t);
	}
}
