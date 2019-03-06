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

void allaire_diffuse :: set_parameters (std::string test_case, sim_info& params, binarySGparams& eosparams)
{
	if (test_case == "ST1_x" || test_case == "ST1_y")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 2.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;

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
	else if (test_case == "ST2_x" || test_case == "ST2_y")
	{
		eosparams.gamma1 = 4.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 600000000.0;
		eosparams.pinf2 = 0.0;

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
	else if (test_case == "TTC4_x_pure" || test_case == "TTC4_y_pure")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "TTC1_x_pure" || test_case == "TTC1_y_pure")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "TTC2_x_pure" || test_case == "TTC2_y_pure")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "TTC3_x_pure" || test_case == "TTC3_y_pure")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "TTC5_x_pure" || test_case == "TTC5_y_pure")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "circular_explosion")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "underwater_shocked_bubble")
	{
		// See "Practical techniques in ghost fluid method..", Xu, Communications in computational physics, 2016
		
		eosparams.gamma1 = 7.15;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 3309.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "underwater_explosion")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 7.15;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 3.309e8;
		
		params.x0 = -5.0;
		params.y0 = -5.0;
		params.dx = 10.0/params.Nx;
		params.dy = 10.0/params.Ny;
		params.T = 0.005;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else if (test_case == "underwater_explosion_modified")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 7.15;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 3.309e8;
		
		params.x0 = -5.0;
		params.y0 = -5.0;
		params.dx = 10.0/params.Nx;
		params.dy = 10.0/params.Ny;
		params.T = 0.010;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (test_case == "shocked_helium_bubble")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.667;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 325.0/params.Nx;
		params.dy = 89.0/params.Ny;
		params.T = 280.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (test_case == "shocked_R22_bubble")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.249;

		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 0.445/params.Nx;
		params.dy = 0.089/params.Ny;
		params.T = 0.001080;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (test_case == "shocked_SF6")
	{
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.076;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "RMI_SF6")
	{
		// From "A volume of fluid method based ghost fluid method for compressible multi-fluid flows" - Computers & Fluids - 2014
		
		eosparams.gamma1 = 1.4;
		eosparams.gamma2 = 1.093;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 4.0/params.Nx;
		params.dy = 0.5/params.Ny;
		params.T = 10.0;
		params.BC_L = "transmissive";
		params.BC_T = "reflective";
		params.BC_R = "transmissive";
		params.BC_B = "reflective";
	}
	else if (test_case == "tin_air_implosion")
	{
		eosparams.gamma1 = 3.27;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 149500.0;
		eosparams.pinf2 = 0.0;
		
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
	else if (test_case == "TSTM")
	{
		eosparams.gamma1 = 1.5;
		eosparams.gamma2 = 1.4;
		eosparams.pinf1 = 0.0;
		eosparams.pinf2 = 0.0;
		
		params.x0 = 0.0;
		params.y0 = 0.0;
		params.dx = 7.0/params.Nx;
		params.dy = 3.0/params.Ny;
		params.T = 8.0;
		params.BC_L = "transmissive";
		params.BC_T = "transmissive";
		params.BC_R = "transmissive";
		params.BC_B = "transmissive";
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid test_case in settings file.");
	}
}

void allaire_diffuse :: set_halfspace_IC (const vectype& U_under, const double a, const double b, const double c, gridtype& grid, const sim_info& params)
{
	/*
	 * Set every cell where a*x + b*y <= c to the state U_under.
	 */
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			
			if (cc(0)*a + cc(1)*b <= c)
			{
				grid[i][j] = U_under;
			}
		}
	}
}

void allaire_diffuse :: set_planar_IC (const vectype& U_under, const vectype& U_over, const double a, const double b, const double c, gridtype& grid, const sim_info& params)
{
	/*
	 * Set every cell where a*x + b*y <= c to the state U_under,
	 * otherwise U_over. Useful for planar initial conditions.
	 */
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			
			if (cc(0)*a + cc(1)*b <= c)
			{
				grid[i][j] = U_under;
			}
			else
			{
				grid[i][j] = U_over;
			}
		}
	}
}


void allaire_diffuse :: set_circular_IC (const vectype& W_in, const vectype& W_out, const Eigen::Vector2d& centre, const double R, gridtype& grid, const sim_info& params, const int N)
{
	/*
	 * Set state according to weighting of volume fraction inside/outside
	 * the circle, according to Monte Carlo estimate using N^2 samples
	 */
	 
	for (int i=0; i<params.Ny + 2 * params.numGC; i++)
	{
		for (int j=0; j<params.Nx + 2 * params.numGC; j++)
		{
			int totalnumsamples = N*N;
			int numinside = 0;
			double delx = params.dx/N;
			double dely = params.dy/N;
			
			Eigen::Vector2d cc = params.cellcentre_coord(i, j);
			Eigen::Vector2d BL;
			BL(0) = cc(0) - 0.5 * params.dx;
			BL(1) = cc(1) - 0.5 * params.dy;
			Eigen::Vector2d samplepos;
			
			for (int a=0; a<N; a++)
			{
				for (int b=0; b<N; b++)
				{
					samplepos(0) = BL(0) + (a + 0.5) * delx;
					samplepos(1) = BL(1) + (b + 0.5) * dely;
					
					samplepos -= centre;
					
					if (samplepos.norm() <= R) numinside++;
				}
			}

			double frac = double(numinside)/totalnumsamples;
			
			vectype W = frac * W_in + (1.0 - frac) * W_out;
			grid[i][j] = primitives_to_conserved(eosparams, W);
		}
	}
}


std::shared_ptr<gridtype> allaire_diffuse :: set_ICs (settings_file SF, sim_info& params)
{
	params.Nx = SF.Nx;
	params.Ny = SF.Ny;
	params.CFL = SF.CFL;
	params.outputname = SF.basename;
	params.output_freq = SF.output_freq;
	set_parameters(SF.test_case, params, eosparams);
	
	
	// Number of ghost cells needed
	
	if (SF.flux_solver == "Godunov")
	{
		params.numGC = 1;
		params.stclsize = 1;
	}
	else if (SF.flux_solver == "MUSCL")
	{
		params.numGC = 2;
		params.stclsize = 2;
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid flux_solver in settings file.");
	}
	
	
	// Riemann solver object
	
	std::shared_ptr<riemann_solver_base> RS_ptr = nullptr;	
	
	if (SF.riemann_solver == "HLLC")
	{
		RS_ptr = std::make_shared<HLLC_riemann_solver>(params, eosparams.gamma1, eosparams.gamma2, eosparams.pinf1, eosparams.pinf2);
	}
	else
	{
		assert(!"[allaire_diffuse] Invalid riemann_solver in settings file.");
	}
	
	
	// Flux solver object
	
	if (SF.flux_solver == "Godunov")
	{
		FS_ptr = std::make_shared<flux_solver_godunov>(RS_ptr, params, eosparams.gamma1, eosparams.gamma2, eosparams.pinf1, eosparams.pinf2);
		zupdate_ptr = std::make_shared<zupdate_upwind>();
	}
	else if (SF.flux_solver == "MUSCL")
	{
		FS_ptr = std::make_shared<flux_solver_MUSCLHANCOCK>(RS_ptr, params, eosparams.gamma1, eosparams.gamma2, eosparams.pinf1, eosparams.pinf2);
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
		vectype Uleft = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 0.0, 0.125, 0.0, 0.0, 0.1, 0.0);
		
		if (SF.test_case == "ST1_x")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "ST1_y")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "ST2_x" || SF.test_case == "ST2_y")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 1000.0, 0.0, 0.0, 0.0, 1000000000.0, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 0.0, 50.0, 0.0, 0.0, 100000.0, 0.0);
		
		if (SF.test_case == "ST2_x")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.7, ICgrid, params);
		}
		else if (SF.test_case == "ST2_y")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.7, ICgrid, params);
		}
	}
	else if (SF.test_case == "TTC5_x_pure" || SF.test_case == "TTC5_y_pure")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 5.99924, 0.0, 19.5975, 0.0, 460.894, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 5.99242, 0.0, -6.19633, 0.0, 46.0950, 1.0);
		
		if (SF.test_case == "TTC5_x_pure")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "TTC5_y_pure")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "TTC1_x_pure" || SF.test_case == "TTC1_y_pure")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 0.125, 0.0, 0.0, 0.0, 0.1, 1.0);
		
		if (SF.test_case == "TTC1_x_pure")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "TTC1_y_pure")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "TTC2_x_pure" || SF.test_case == "TTC2_y_pure")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 1.0, 0.0, -2.0, 0.0, 0.4, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 1.0, 0.0, 2.0, 0.0, 0.4, 1.0);
		
		if (SF.test_case == "TTC2_x_pure")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "TTC2_y_pure")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "TTC3_x_pure" || SF.test_case == "TTC3_y_pure")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 1000.0, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 0.01, 1.0);
		
		if (SF.test_case == "TTC3_x_pure")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "TTC3_y_pure")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "TTC4_x_pure" || SF.test_case == "TTC4_y_pure")
	{
		vectype Uleft = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 0.1, 1.0);
		vectype Uright = primitives_to_conserved(eosparams, 1.0, 0.0, 0.0, 0.0, 100.0, 1.0);
		
		if (SF.test_case == "TTC4_x_pure")
		{
			set_planar_IC(Uleft, Uright, 1.0, 0.0, 0.5, ICgrid, params);
		}
		else if (SF.test_case == "TTC4_y_pure")
		{
			set_planar_IC(Uleft, Uright, 0.0, 1.0, 0.5, ICgrid, params);
		}
	}
	else if (SF.test_case == "circular_explosion")
	{
		vectype Win (6);
		Win << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
		vectype Wout (6);
		Wout << 0.0, 0.125, 0.0, 0.0, 0.1, 0.0;
		
		Eigen::Vector2d centre;
		centre << 1.0, 1.0;
		
		double R = 0.4;
		
		set_circular_IC(Win, Wout, centre, R, ICgrid, params, 10);		
	}
	else if (SF.test_case == "shocked_helium_bubble")
	{
		vectype W_helium (6);
		W_helium << 0.0, 0.138, 0.0, 0.0, 1.0, 0.0;
		
		vectype W_preshock (6);
		W_preshock << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
		
		vectype W_postshock (6);
		W_postshock << 1.3764, 0.0, -0.394, 0.0, 1.5698, 1.0;
		
		Eigen::Vector2d centre;
		centre << 175.0, 44.5;
		
		double R = 25.0;
		
		set_circular_IC(W_helium, W_preshock, centre, R, ICgrid, params, 10);
		
		set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), -1.0, 0.0, -225.0, ICgrid, params);	
	}
	else if (SF.test_case == "shocked_R22_bubble")
	{
		vectype W_R22 (6);
		W_R22 << 0.0, 3.863, 0.0, 0.0, 1.01325e5, 0.0;
		
		vectype W_preshock (6);
		W_preshock << 1.225, 0.0, 0.0, 0.0, 1.01325e5, 1.0;
		
		vectype W_postshock (6);
		W_postshock << 1.686, 0.0, -113.5, 0.0, 1.59e5, 1.0;
		
		Eigen::Vector2d centre;
		centre << 0.225, 0.0445;
		
		double R = 0.025;
		
		set_circular_IC(W_R22, W_preshock, centre, R, ICgrid, params, 25);
		
		set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), -1.0, 0.0, -0.275, ICgrid, params);	
	}
	else if (SF.test_case == "underwater_shocked_bubble")
	{
		vectype W_air (6);
		W_air << 0.0, 0.0012, 0.0, 0.0, 1.0, 0.0;
		
		vectype W_preshock (6);
		W_preshock << 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
		
		vectype W_postshock (6);
		W_postshock << 1.31, 0.0, 67.32, 0.0, 19000.0, 1.0;
		
		Eigen::Vector2d centre;
		centre << 6.0, 6.0;
		
		double R = 3.0;
		
		set_circular_IC(W_air, W_preshock, centre, R, ICgrid, params, 10);
		
		set_halfspace_IC(primitives_to_conserved(eosparams, W_postshock), 1.0, 0.0, 2.4, ICgrid, params);
	}
	else if (SF.test_case == "underwater_explosion" || SF.test_case == "underwater_explosion_modified")
	{
		vectype W_water (6);
		W_water << 0.0, 1000.0, 0.0, 0.0, 1.0e5, 0.0;
		
		vectype W_air (6);
		W_air << 1.0, 0.0, 0.0, 0.0, 1.0e5, 1.0;
		
		vectype W_airbubble (6);
		W_airbubble << 1270.0, 0.0, 0.0, 0.0, 8.29e8, 1.0;
		
		Eigen::Vector2d centre;
		centre << 0.0, 0.0;
		
		double R = 1.0;
		
		set_circular_IC(W_airbubble, W_water, centre, R, ICgrid, params, 10);
		
		set_halfspace_IC(primitives_to_conserved(eosparams, W_air), 0.0, -1.0, -2.5, ICgrid, params);
	}
	else if (SF.test_case == "shocked_SF6")
	{
		double rho_preshock = 1.153;
		double p_preshock = 9.6856;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.6672;
		double p_postshock = 16.3256;
		double u_postshock = 1.33273;
		double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
		double rho2 = 5.805;
		double p2 = 9.6856;
		double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
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
	else if (SF.test_case == "TSTM")
	{
		double rho_preshock = 0.125;
		double p_preshock = 0.1;
		double u_preshock = 0.0;
		double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.0;
		double p_postshock = 1.0;
		double u_postshock = 0.0;
		double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
		double rho2 = 1.0;
		double p2 = 0.1;
		double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
		double v = 0.0;
		double z;
		
		for (int i=0; i<params.Ny + 2 * params.numGC; i++)
		{
			for (int j=0; j<params.Nx + 2 * params.numGC; j++)
			{
				Eigen::Vector2d cc = params.cellcentre_coord(i, j);	
				
				double rho1, u, e1;
				
				if (cc(0) > 1.0)
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
				
				
				// Set z as fraction of area inside rectangular region [1.0, 8.0] x [0.0, 1.5]
				
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
						
						if (samplepos(0) >= 1.0 && samplepos(0) <= 8.0
							&& samplepos(1) >= -1.0 && samplepos(1) <= 1.5) 
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
		double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
		double rho_postshock = 1.411;
		double p_postshock = 1.628;
		double u_postshock = -0.39;
		double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
		double rho2 = 5.04;
		double p2 = 1.0;
		double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
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
				double eps = 0.2;
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
		double e_preshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_preshock, rho_preshock);
		double rho_postshock = 11.84;
		double p_postshock = 1000000.0;
		double e_postshock = eos::specific_ie(eosparams.gamma1, eosparams.pinf1, p_postshock, rho_postshock);
		double rho2 = 0.001;
		double p2 = 1.0;
		double e2 = eos::specific_ie(eosparams.gamma2, eosparams.pinf2, p2, rho2);
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
						double r_interface = 20.0 + 0.4 * cos(22 * theta) + 0.4 * cos(17 * theta) + 0.3 * cos(29 * theta);
						
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
			double p = allairemodel::mixture_pressure(eosparams, rho, e, z);
			outfile5 << p << "\n";
		}
	}
	
	outfile5 << "FIELD FieldData 1\n";
	outfile5 << "Soundspeed 1 " + std::to_string(params.Nx*params.Ny) + " double\n";
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
			outfile5 << c << "\n";
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


void allaire_diffuse :: gnuplot_schlieren (const gridtype& grid, const sim_info& params, int n, double t)
{

	std::string filename2 = params.outputname + "-schlieren-" + std::to_string(t) + ".dat";
	std::ofstream outfile2;
	outfile2.open(filename2);

	std::vector<double> allgrads;
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			double rho_L = grid[i][j-1](0) + grid[i][j-1](1);
			double rho_R = grid[i][j+1](0) + grid[i][j+1](1);
			double rho_T = grid[i+1][j](0) + grid[i+1][j](1);
			double rho_B = grid[i-1][j](0) + grid[i-1][j](1);

			double gradx = (rho_R - rho_L) / (2.0 * params.dx);
			double grady = (rho_T - rho_B) / (2.0 * params.dy);

			allgrads.push_back(sqrt(gradx*gradx + grady*grady));
		}
	}

	double maxgrad = *std::max_element(allgrads.begin(), allgrads.end());

	int counter = 0;

	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);

			outfile2 << CC(0) << " " << CC(1) << " " << allgrads[counter] / maxgrad << std::endl;
			counter++;
		}
		outfile2 << std::endl;
	}

	outfile2.close();
	std::cout << "[allaire_diffuse] Schlieren output to gnuplot complete" << std::endl;
}


void allaire_diffuse :: gnuplot_output (const gridtype& grid, const sim_info& params, int n, double t)
{
	std::string filename = params.outputname + "-state-" + std::to_string(t) + ".dat";
	std::ofstream outfile;
	outfile.open(filename);
	
	for (int i=params.numGC; i<params.Ny + params.numGC; i++)
	{
		for (int j=params.numGC; j<params.Nx + params.numGC; j++)
		{
			Eigen::Vector2d CC = params.cellcentre_coord(i, j);
			
			double rho = grid[i][j](0) + grid[i][j](1);
			double u = grid[i][j](2) / rho;
			double v = grid[i][j](3) / rho;
			double E = grid[i][j](4);
			double e = E / rho - 0.5 * (u * u + v * v);
			double z = grid[i][j](5);
			double p = allairemodel::mixture_pressure(eosparams, rho, e, z);
			//double c = allairemodel::mixture_soundspeed(eosparams, rho, p, z);
			
			outfile << CC(0) << " " << CC(1) << " " << rho << " " << u << " " << v << " " << e << " " << z << " " << p << std::endl;
		}
		outfile << std::endl;
	}
	
	outfile.close();
	std::cout << "[allaire_diffuse] State output to gnuplot complete" << std::endl;
	
	
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
		double p = allairemodel::mixture_pressure(eosparams, rho, e, z);
			
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
		double p = allairemodel::mixture_pressure(eosparams, rho, e, z);
			
		outfile2 << y << " " << rho << " " << u << " " << v << " " << e << " " << p << " " << z << std::endl;
	}
	
	outfile2.close();
	
	std::cout << "[allaire_diffuse] Lineout to gnuplot complete" << std::endl;
}

void allaire_diffuse :: gnuplot_masschange (const sim_info& params)
{
	std::string filename = params.outputname + "-masschange.dat";
	std::ofstream outfile;
	outfile.open(filename);
		
	for (unsigned int k=0; k<time.size(); k++)
	{
		outfile << time[k] << " " << (mass1[k] - mass1[0])/mass1[0] << " " << (mass2[k] - mass2[0])/mass2[0] << std::endl;
	}
	
	outfile.close();
	
	std::cout << "[allaire_diffuse] Mass change output complete" << std::endl;
}
	
void allaire_diffuse :: output (const gridtype& grid, const sim_info& params, int n, double t)
{
	const double outputinterval = params.T / params.output_freq;
	static double lastoutputtime = 0.0;


	// Hard-code some specific output times :(
	
	static std::vector<double> output_times {0.000115, 0.000175, 0.000247, 0.000307, 0.000378, 0.000402, 0.000477, 0.001080};

	if (!output_times.empty())
	{
		if (t > output_times[0])
		{
			gnuplot_schlieren(grid, params, n, t);
			output_times.erase(output_times.begin());
		}
	}


	if (n==0 || t==params.T)
	{
		//vtk_output(grid, params, n, t);
		gnuplot_output(grid, params, n, t);
		gnuplot_schlieren(grid, params, n, t);
		gnuplot_lineout(grid, params, n, t);
		
		if (t == params.T)
		{
			gnuplot_masschange(params);
		}
			
	}
	else if (params.output_freq != 0.0)
	{
		if (t - lastoutputtime > outputinterval)
		{
			//vtk_output(grid, params, n, t);
			gnuplot_output(grid, params, n, t);
			gnuplot_schlieren(grid, params, n, t);
			lastoutputtime = t;
		}
	}
}
