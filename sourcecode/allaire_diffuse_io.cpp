#include "allaire_diffuse.hpp"
#include "stiffened_gas_eos.hpp"
#include "flux_solver_godunov.hpp"
#include "flux_solver_MUSCLhancock.hpp"
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
	params.output_freq = SF.output_freq;
	
	std::shared_ptr<riemann_solver_base> RS_ptr = nullptr;
	
	
	// EOS parameters, computational domain and boundary conditions
	
	if (SF.test_case == "ST1_x" || SF.test_case == "ST1_y")
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
	else if (SF.test_case == "ST2_x" || SF.test_case == "ST2_y")
	{
		gamma1 = 4.4;
		gamma2 = 1.4;
		pinf1 = 600000000.0;
		pinf2 = 0.0;

		params.x0 = -2.0;
		params.y0 = -2.0;
		params.dx = 4.0/params.Nx;
		params.dy = 4.0/params.Ny;
		params.T = 0.0009;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "TTC4_x_pure" || SF.test_case == "TTC4_y_pure")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.035;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "TTC1_x_pure" || SF.test_case == "TTC1_y_pure")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.25;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "TTC2_x_pure" || SF.test_case == "TTC2_y_pure")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.15;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "TTC3_x_pure" || SF.test_case == "TTC3_y_pure")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.012;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "TTC5_x_pure" || SF.test_case == "TTC5_y_pure")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 1.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 0.035;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "circular_explosion")
	{
		gamma1 = 1.4;
		gamma2 = 1.4;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 2.0/params.Nx;
		params.dy = 2.0/params.Ny;
		params.T = 0.25;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		// See "Practical techniques in ghost fluid method..", Xu, Communications in computational physics, 2016
		
		gamma1 = 7.15;
		gamma2 = 1.4;
		pinf1 = 3309.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 12.0/params.Nx;
		params.dy = 12.0/params.Ny;
		params.T = 0.05;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "underwater_explosion")
	{
		gamma1 = 1.4;
		gamma2 = 7.15;
		pinf1 = 0.0;
		pinf2 = 3.309e8;
		
		params.x0 = -5.0;
		params.y0 = -5.0;
		params.dx = 10.0/params.Nx;
		params.dy = 10.0/params.Ny;
		params.T = 0.003;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		gamma1 = 1.4;
		gamma2 = 1.667;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 325.0/params.Nx;
		params.dy = 89.0/params.Ny;
		params.T = 250.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (SF.test_case == "shocked_SF6")
	{
		gamma1 = 1.4;
		gamma2 = 1.076;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 0.45/params.Nx;
		params.dy = 0.2/params.Ny;
		params.T = 0.25;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "reflective";
		params.BC_B = "reflective";
	}
	else if (SF.test_case == "RMI_SF6")
	{
		// From "A volume of fluid method based ghost fluid method for compressible multi-fluid flows" - Computers & Fluids - 2014
		
		gamma1 = 1.4;
		gamma2 = 1.093;
		pinf1 = 0.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 4.0/params.Nx;
		params.dy = 1.0/params.Ny;
		params.T = 10.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (SF.test_case == "tin_air_implosion")
	{
		gamma1 = 3.27;
		gamma2 = 1.4;
		pinf1 = 149500.0;
		pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = -25.0;
		params.dx = 25.0/params.Nx;
		params.dy = 50.0/params.Ny;
		params.T = 0.08;
		params.BC_L = "reflective";
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
	else if (SF.flux_solver == "MUSCL1" || SF.flux_solver == "MUSCL2")
	{
		params.numGC = 2;
		params.stclsize = 2;
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
	else if (SF.flux_solver == "MUSCL1")
	{
		FS_ptr = std::make_shared<flux_solver_MUSCLHANCOCK_1>(RS_ptr, params, gamma1, gamma2, pinf1, pinf2);
		zupdate_ptr = std::make_shared<zupdate_secondorder>();
	}
	else if (SF.flux_solver == "MUSCL2")
	{
		FS_ptr = std::make_shared<flux_solver_MUSCLHANCOCK_2>(RS_ptr, params, gamma1, gamma2, pinf1, pinf2);
		zupdate_ptr = std::make_shared<zupdate_secondorder>();
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid flux_solver in settings file.");
	}
			
	
	gridtype ICgrid (SF.Ny + 2 * params.numGC, rowtype(SF.Nx + 2 * params.numGC, vectype(6)));
	
	
	// Initial states
	
	if (SF.test_case == "ST1_x" || SF.test_case == "ST1_y")
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
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "ST1_x")
				{
					if (cc(0) < 0.5)
					{
						z = 1.0;
					}
					else
					{
						z = 0.0;
					}
				}
				else
				{
					if (cc(1) < 0.5)
					{
						z = 1.0;
					}
					else
					{
						z = 0.0;
					}
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
	else if (SF.test_case == "ST2_x" || SF.test_case == "ST2_y")
	{
		double rho1 = 1000.0;
		double p1 = 1000000000.0;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 50.0;
		double p2 = 100000.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u = 0.0;
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				
				if (SF.test_case == "ST2_x")
				{
					if (cc(0) < 0.7)
					{
						// Fluid 1 here
						z = 1.0;
					}
					else
					{
						// Fluid 2 here
						z = 0.0;
					}
				}
				else
				{
					if (cc(1) < 0.7)
					{
						// Fluid 1 here
						z = 1.0;
					}
					else
					{
						// Fluid 2 here
						z = 0.0;
					}
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
	else if (SF.test_case == "TTC5_x_pure" || SF.test_case == "TTC5_y_pure")
	{
		double rho1 = 5.99924;
		double p1 = 460.894;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 5.99242;
		double p2 = 46.0950;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u1 = 19.5975;
		double u2 = -6.19633;
		double v = 0.0;
		double z = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "TTC5_x_pure")
				{
					if (cc(0) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
				else
				{					
					if (cc(1) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho1;
						ICgrid[i][j](3) = u1 * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho2;
						ICgrid[i][j](3) = u2 * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
			}
		}
	}
	else if (SF.test_case == "TTC1_x_pure" || SF.test_case == "TTC1_y_pure")
	{
		double rho1 = 1.0;
		double p1 = 1.0;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 0.125;
		double p2 = 0.1;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u1 = 0.0;
		double u2 = 0.0;
		double v = 0.0;
		double z = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "TTC1_x_pure")
				{
					if (cc(0) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
				else
				{					
					if (cc(1) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
			}
		}
	}
	else if (SF.test_case == "TTC2_x_pure" || SF.test_case == "TTC2_y_pure")
	{
		double rho1 = 1.0;
		double p1 = 0.4;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 1.0;
		double p2 = 0.4;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u1 = -2.0;
		double u2 = 2.0;
		double v = 0.0;
		double z = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "TTC2_x_pure")
				{
					if (cc(0) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
				else
				{					
					if (cc(1) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho1;
						ICgrid[i][j](3) = u1 * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho2;
						ICgrid[i][j](3) = u2 * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
			}
		}
	}
	else if (SF.test_case == "TTC3_x_pure" || SF.test_case == "TTC3_y_pure")
	{
		double rho1 = 1.0;
		double p1 = 1000.0;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 1.0;
		double p2 = 0.01;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u1 = 0.0;
		double u2 = 0.0;
		double v = 0.0;
		double z = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "TTC3_x_pure")
				{
					if (cc(0) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
				else
				{					
					if (cc(1) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho1;
						ICgrid[i][j](3) = u1 * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho2;
						ICgrid[i][j](3) = u2 * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
			}
		}
	}
	else if (SF.test_case == "TTC4_x_pure" || SF.test_case == "TTC4_y_pure")
	{
		double rho1 = 1.0;
		double p1 = 0.1;
		double e1 = eos::specific_ie(gamma1, pinf1, p1, rho1);
		double rho2 = 1.0;
		double p2 = 100.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double u1 = 0.0;
		double u2 = 0.0;
		double v = 0.0;
		double z = 1.0;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				if (SF.test_case == "TTC4_x_pure")
				{
					if (cc(0) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u1 * rho1;
						ICgrid[i][j](3) = v * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = u2 * rho2;
						ICgrid[i][j](3) = v * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
				else
				{					
					if (cc(1) < 0.5)
					{
						ICgrid[i][j](0) = z * rho1;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho1;
						ICgrid[i][j](3) = u1 * rho1;
						ICgrid[i][j](4) = rho1 * e1 + 0.5 * rho1 * (u1*u1 + v*v);
						ICgrid[i][j](5) = z;
					}
					else
					{
						ICgrid[i][j](0) = z * rho2;
						ICgrid[i][j](1) = 0.0;
						ICgrid[i][j](2) = v * rho2;
						ICgrid[i][j](3) = u2 * rho2;
						ICgrid[i][j](4) = rho2 * e2 + 0.5 * rho2 * (u2*u2 + v*v);
						ICgrid[i][j](5) = z;
					}
				}
			}
		}
	}
	else if (SF.test_case == "circular_explosion")
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
				// Set z as fraction of area inside circle of radius 0.4 at (1.0, 1.0)
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos(0) -= 1.0;
						samplepos(1) -= 1.0;
						
						if (samplepos.norm() <= 0.4) numinside++;
					}
				}

				z = double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		double rho_preshock = 1.0;
		double p_preshock = 1.0;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.3764;
		double p_postshock = 1.5698;
		double u_postshock = -0.394;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 0.138;
		double p2 = 1.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) < 225.0)
				{
					rho1 = rho_preshock;
					u = u_preshock;
					e1 = e_preshock;
				}
				else
				{
					rho1 = rho_postshock;
					u = u_postshock;
					e1 = e_postshock;
				}
				
				// Set z as fraction of area inside circle of radius 25 at (175, 44.5)
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos(0) -= 175.0;
						samplepos(1) -= 44.5;
						
						if (samplepos.norm() <= 25.0) numinside++;
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		double rho_preshock = 1.0;
		double p_preshock = 1.0;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.31;
		double p_postshock = 19000.0;
		double u_postshock = 67.32;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 0.0012;
		double p2 = 1.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) < 2.4)
				{
					rho1 = rho_preshock;
					u = u_preshock;
					e1 = e_preshock;
				}
				else
				{
					rho1 = rho_postshock;
					u = u_postshock;
					e1 = e_postshock;
				}
				
				// Set z as fraction of area inside circle of radius 3 at (6, 6)
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos(0) -= 6.0;
						samplepos(1) -= 6.0;
						
						if (samplepos.norm() <= 3.0) numinside++;
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "underwater_explosion")
	{
		double rho_preshock = 1.0;
		double p_preshock = 1.0e5;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1270.0;
		double p_postshock = 8.29e8;
		double u_postshock = 0.0;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 1000.0;
		double p2 = 1.0e5;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) > -0.1)
				{
					rho1 = rho_preshock;
					u = u_preshock;
					e1 = e_preshock;
				}
				else
				{
					rho1 = rho_postshock;
					u = u_postshock;
					e1 = e_postshock;
				}
				
				// Set z as fraction of area inside circle of radius 1 at (0, 0) or above line y=2.5
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						samplepos(0) -= 0.0;
						samplepos(1) -= 0.0;
						
						if (samplepos.norm() <= 1.0 || samplepos(1) > 2.5)
						{ 
							numinside++;
						}
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "shocked_SF6")
	{
		double rho_preshock = 1.153;
		double p_preshock = 9.6856;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.6672;
		double p_postshock = 16.3256;
		double u_postshock = 1.33273;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 5.805;
		double p2 = 9.6856;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) > 0.05)
				{
					rho1 = rho_preshock;
					u = u_preshock;
					e1 = e_preshock;
				}
				else
				{
					rho1 = rho_postshock;
					u = u_postshock;
					e1 = e_postshock;
				}
				
				
				// Set z as fraction of area inside rectangular region [0.1, 0.25] x [0.0, 0.1]
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						if (samplepos(0) >= 0.1 && samplepos(0) <= 0.25
							&& samplepos(1) >= 0.0 && samplepos(1) <= 0.1) 
						{
							numinside++;
						}
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "RMI_SF6")
	{
		double rho_preshock = 1.0;
		double p_preshock = 1.0;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.411;
		double p_postshock = 1.628;
		double u_postshock = -0.39;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 5.04;
		double p2 = 1.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) < 3.2)
				{
					rho1 = rho_preshock;
					u = u_preshock;
					e1 = e_preshock;
				}
				else
				{
					rho1 = rho_postshock;
					u = u_postshock;
					e1 = e_postshock;
				}
				
				
				// Set z as fraction of area inside sinusoidal interface
				
				double x1 = 2.9;
				double eps = 0.1;
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double pi = atan(1.0) * 4.0;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						double interfacex = x1 - eps * sin(2.0 * pi * (samplepos(1) + 0.25));
						
						
						if (interfacex > samplepos(0)) 
						{
							numinside++;
						}
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
				ICgrid[i][j](0) = z * rho1;
				ICgrid[i][j](1) = (1.0 - z) * rho2;
				ICgrid[i][j](2) = u * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](3) = v * (z * rho1 + (1.0 - z) * rho2);
				ICgrid[i][j](4) = z * rho1 * e1 + (1.0 - z) * rho2 * e2 + 0.5 * (z * rho1 + (1.0 - z) * rho2) * (u*u + v*v);
				ICgrid[i][j](5) = z;
			}
		}		
	}
	else if (SF.test_case == "tin_air_implosion")
	{
		double rho_preshock = 7.28;
		double p_preshock = 1.0;
		double e_preshock = eos::specific_ie(gamma1, pinf1, p_preshock, rho_preshock);
		double rho_postshock = 11.84;
		double p_postshock = 1000000.0;
		double e_postshock = eos::specific_ie(gamma1, pinf1, p_postshock, rho_postshock);
		double rho2 = 0.001;
		double p2 = 1.0;
		double e2 = eos::specific_ie(gamma2, pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				// Set tin state as fraction of area inside circle of radius 24
				
				int numsamples = 10;
				int totalnumsamples = numsamples*numsamples;
				int numinside = 0;
				double delx = params.dx/numsamples;
				double dely = params.dy/numsamples;
				
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);
				Eigen::Vector2d BL;
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						if (samplepos.norm() <= 24.0) numinside++;
					}
				}
				
				double insideratio = double(numinside) / totalnumsamples;
				
				double rho1, u = 0.0, e1;
				rho1 = insideratio * rho_preshock + (1.0 - insideratio) * rho_postshock;
				e1 = insideratio * e_preshock + (1.0 - insideratio) * e_postshock;
				
				
				// Set z as fraction of area inside circle of radius 20 at (0, 0)
				
				numsamples = 10;
				totalnumsamples = numsamples*numsamples;
				numinside = 0;
				delx = params.dx/numsamples;
				dely = params.dy/numsamples;
				
				cc = params.cellcentre_coord(i, j);
				BL(0) = cc(0) - 0.5 * params.dx;
				BL(1) = cc(1) - 0.5 * params.dy;
				
				for (int a=0; a<numsamples; a++)
				{
					for (int b=0; b<numsamples; b++)
					{
						Eigen::Vector2d samplepos;
						samplepos(0) = BL(0) + (a + 0.5) * delx;
						samplepos(1) = BL(1) + (b + 0.5) * dely;
						
						double theta = atan2(samplepos(1), samplepos(0));
						theta -= atan(1) * 2;
						double r_interface = 20.0 + 0.4 * cos(31 * theta) + 0.4 * cos(25 * theta) + 0.4 * cos(39 * theta);
						
						if (samplepos.norm() <= r_interface) numinside++;
					}
				}

				z = 1.0 - double(numinside)/totalnumsamples;
				
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
	
	std::cout << "[allaire_diffuse] Output to vtk complete" << std::endl;
}

void allaire_diffuse :: gnuplot_lineout (const gridtype& grid, const sim_info& params, int n, double t)
{
	std::string filename = params.outputname + "-lineoutx-" + std::to_string(n) + ".dat";
	std::ofstream outfile;
	outfile.open(filename);
	
	int i = params.Ny / 2;
	
	for (int j=0; j<params.Nx + 2 * params.numGC; j++)
	{
		double x = params.x0 + double(j - params.numGC)*params.dx + 0.5*params.dx;
		double rho = grid[i][j](0) + grid[i][j](1);
		double u = grid[i][j](2) / rho;
		double v = grid[i][j](3) / rho;
		double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
		double z = grid[i][j](5);
		double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, z);
			
		outfile << x << " " << rho << " " << u << " " << v << " " << e << " " << p << " " << z << std::endl;
	}
	
	outfile.close();
	
	filename = params.outputname + "-lineouty-" + std::to_string(n) + ".dat";
	std::ofstream outfile2;
	outfile2.open(filename);
	
	int j = params.Nx / 2;
	
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		double y = params.y0 + double(i - params.numGC)*params.dy + 0.5*params.dy;
		double rho = grid[i][j](0) + grid[i][j](1);
		double u = grid[i][j](2) / rho;
		double v = grid[i][j](3) / rho;
		double e = grid[i][j](4) / rho - 0.5 * (u * u + v * v);
		double z = grid[i][j](5);
		double p = allairemodel::mixture_pressure(gamma1, gamma2, pinf1, pinf2, rho, e, z);
			
		outfile2 << y << " " << rho << " " << u << " " << v << " " << e << " " << p << " " << z << std::endl;
	}
	
	outfile2.close();
	
	std::cout << "[allaire_diffuse] Output to gnuplot complete" << std::endl;
}
	
void allaire_diffuse :: output (const gridtype& grid, const sim_info& params, int n, double t)
{
	if (n==0 || t==params.T)
	{
		vtk_output(grid, params, n, t);
		gnuplot_lineout(grid, params, n, t);
	}
	else if (params.output_freq != 0.0)
	{
		if (n % params.output_freq == 0)
		{
			vtk_output(grid, params, n, t);
		}
	}
}
