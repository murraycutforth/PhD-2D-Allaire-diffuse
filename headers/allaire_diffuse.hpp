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

class allaire_diffuse : public problem_base {
	
	public:
	
	// EOS parameters of two fluids
	
	double gamma1;
	double gamma2;
	double pinf1;
	double pinf2;
	
	
	// Functions specific to this problem
	
	void set_boundary_conditions (gridtype& grid, const sim_info& params);
	
	
	// Over-ride all pure virtual member functions of problem_base
	
	std::shared_ptr<gridtype> set_ICs (settings_file SF, sim_info& params);
	
	void output (const gridtype& grid, const sim_info& params, int n, double t);
	
	double compute_dt (const gridtype& grid, const sim_info& params, double t);
	
	void pre_sweep (gridtype& grid, const sim_info& params);
	
	void update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t);
	
	void update_col (const gridtype& grid, gridtype& future_grid, const sim_info& params, int j, double dt, double t);
	
	void post_sweep (gridtype& grid, gridtype& future_grid, const sim_info& params);
	
};

#endif
