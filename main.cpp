#include <iostream>
#include <cmath>
#include <cfloat>
#include <omp.h> // add OpenMP support in Configuration Properties
#define NUM_THREADS omp_get_max_threads() // number of threads computing in parallel

#include "common.h"
#include "synthData.h"
#include "relpose4p3v.h"



int main()
{
	Hist hRt;
	hRt.iniHist(-15, 2);
	double speed=0;
	Timer timer;
	
	// --------------------------------------------------------------------------------------------
	
	int ntrials=0;
	
	for (int i=1; i<=NTRIALS; ++i)
	{
		double data[_V][3][_P];
		Camera cam_gt, cam_est;
		
		synthData(data, cam_gt); // generate synthetic data
		
		timer.start(); // start timer
		
		if (!relative4p3v(data, cam_est)) continue; // run solver

		speed+=timer.stop(); // stop timer

		hRt.updateHist(numError(cam_est.Rt, cam_gt.Rt)); // update error histogram

		++ntrials;
	}

	// results
	std::cout.precision(4);
	std::cout << "\nNumber of successful trials: " << ntrials << "\n\n";
	std::cout << "Average runtime (ms): " << speed*((1e+3)/(double)ntrials) << "\n\n";
	std::cout << "Median numerical error: " << hRt.getMedian(0.5) << "\n\n";
	std::cout << "Mean numerical error: " << hRt.getMean() << "\n\n";
	std::cout << "Fails (%): " << hRt.getFails(-1.) << "\n\n";
	std::cout << "Numerical error distribution:\n"; hRt.printHist();

	return 0;
}