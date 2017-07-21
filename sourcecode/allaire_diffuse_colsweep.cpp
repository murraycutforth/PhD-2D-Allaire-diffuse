#include "allaire_diffuse.hpp"
#include "sim_info.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <string>

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
