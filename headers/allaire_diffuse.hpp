/*
 *	DESCRIPTION:	An implementation of the 2D five-equation two
 *			fluid method of Allaire. 
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		21/07/2017
 */


#ifndef ALLAIREDIFFUSEPROBLEM_H
#define ALLAIREDIFFUSEPROBLEM_H

#include "problem_base.hpp"
#include "sim_info.hpp"
#include "typedefs.hpp"
#include "settings_file.hpp"
#include "flux_solver_base.hpp"
#include "stiffened_gas_eos.hpp"
#include "zupdate.hpp"
#include <vector>

class allaire_diffuse : public problem_base {
	
	public:
	
	// EOS parameters of two stiffened gas fluids
	
	binarySGparams eosparams;
	
	
	// Algorithm for updating conservative variables
	
	std::shared_ptr<flux_solver_base> FS_ptr;
	
	
	// Algorithm for updating diffuse volume fractions
	
	std::shared_ptr<zupdate_base> zupdate_ptr;
	
	
	// Storage for mass of each fluid at each time step
	
	std::vector<double> time;
	std::vector<double> mass1;
	std::vector<double> mass2;
	
	
	// Constructors
	
	allaire_diffuse ()
	{}
	
	allaire_diffuse (double gamma1, double gamma2, double pinf1, double pinf2, std::shared_ptr<flux_solver_base> FS_ptr, std::shared_ptr<zupdate_base> zupdate_ptr)
	:
		eosparams (gamma1, gamma2, pinf1, pinf2),
		FS_ptr (FS_ptr),
		zupdate_ptr (zupdate_ptr)
	{}
		
	
	// Functions specific to this problem
	
	void set_boundary_conditions (gridtype& grid, const sim_info& params);
	
	void vtk_output (const gridtype& grid, const sim_info& params, int n, double t);
	
	void gnuplot_lineout (const gridtype& grid, const sim_info& params, int n, double t);
	
	void gnuplot_masschange (const sim_info& params);
	
	double get_rho (const vectype& U)
	{
		return U(0) + U(1);
	}
	
	double get_u (const vectype& U)
	{
		return U(2) / get_rho(U);
	}
	
	double get_v (const vectype& U)
	{
		return U(3) / get_rho(U);
	}
	
	double get_e (const vectype& U)
	{
		return U(4) / get_rho(U) - 0.5 * (get_u(U) * get_u(U) + get_v(U) * get_v(U));
	}
	
	
	// Over-ride all pure virtual member functions of problem_base
	
	std::shared_ptr<gridtype> set_ICs (settings_file SF, sim_info& params);
	
	void output (const gridtype& grid, const sim_info& params, int n, double t);
	
	double compute_dt (const gridtype& grid, const sim_info& params, int n, double t);
	
	void pre_sweep (gridtype& grid, const sim_info& params);
	
	void update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t);
	
	void update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t);
	
	void post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params);
	
	
	// Clone function which returns a shared_ptr<problem_base> object for deep copy of problem between threads
	
	std::shared_ptr<problem_base> clone ()
	{
		return std::make_shared<allaire_diffuse>(eosparams.gamma1, eosparams.gamma2, eosparams.pinf1, eosparams.pinf2, FS_ptr->clone(), zupdate_ptr->clone());
	}
};

#endif
