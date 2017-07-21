#include "allaire_diffuse.hpp"
#include "stiffened_gas_eos.hpp"
#include <cmath>
#include <algorithm>

double allaire_diffuse :: compute_dt (const gridtype& grid, const sim_info& params, double t)
{
	return 0.5;
}
