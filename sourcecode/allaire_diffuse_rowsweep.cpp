#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>

void allaire_diffuse :: update_row (const gridtype& grid, gridtype& future_grid, const sim_info& params, int i, double dt, double t)
{
	// Iterate along every edge bordering a real cell
		// Depending on stencil size (make this a member of params) fill up a vector with states and pass to abstract solver class (which has members of flux, u_star, p_star, riemann solver, EOS params)
	
	// Iterate along every real cell
		// Use stored fluxes (in allaire_diffuse class) to update conserved variables
		// Call volume fraction update method somehow..
	
	
}
