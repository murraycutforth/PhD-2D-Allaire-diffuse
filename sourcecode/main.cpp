/*
 *	DESCRIPTION:	Unit test of the HCL framework which runs a
 *			simulation.
 * 
 * 	AUTHOR:		Murray Cutforth
 * 	DATE:		21/07/2017
 */

#include "settings_file.hpp"
#include "problem_base.hpp"
#include "simulation.hpp"
#include "allaire_diffuse.hpp"
#include <iostream>
#include <memory>
#include <omp.h>

int main()
{
	std::cout << "[main] Beginning simulation.." << std::endl;
	
	settings_file SF;
	SF.read_settings_file("./settingsfile.txt");
	std::shared_ptr<problem_base> problem = std::make_shared<allaire_diffuse>();
	
	simulation sim (SF, problem);
	
	double starttime = omp_get_wtime();
	
	int status = sim.run_simulation();
	
	double endtime = omp_get_wtime();
	
	std::cout << "[main] Simulation exited with status = " << status << std::endl;
	std::cout << "[main] The wall clock runtime was " << endtime - starttime << "s." << std::endl;
	
	return 0;
}
