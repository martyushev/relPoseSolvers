#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <omp.h> // add OpenMP support in Configuration Properties
#define NUM_THREADS omp_get_max_threads() // number of threads

#include "common.h"
#include "math.h"
#include "synthData.h"
#include "relpose4p3v.h"
#include "stats.h"



int main()
{
	Stats stats;

	char dataOut[100];
	sprintf_s(dataOut, "%sdata_%i_%i.txt", FOLDER_OUT, SCENE, MOTION);
	std::ofstream exportData(dataOut, std::ios::out);
	exportData.precision(DBL_DIG);

	for (int i=1; i<=NTRIALS; ++i)
	{
		double data[NVIEWS][3][NPOINTS];
		Camera cam_gt, cam_est;
		
		synthData(data,cam_gt); // generate synthetic data and ground truth cameras
		
		stats.timer.start(); // start timer
		
		if (!relative4p3v(data,cam_est)) continue; // run solver

		stats.totalTime+=stats.timer.stop(); // stop timer

		stats.updateStats(cam_est,cam_gt); // update statistics

		exportData << stats.errNum << "\n"; // export errors to file
	}
	exportData.close();

	stats.printStats(); // print results

	//int n;
	//std::cin >> n;

	return 0;
}