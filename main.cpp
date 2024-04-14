#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <omp.h> // add OpenMP support in Configuration Properties
#define NUM_THREADS omp_get_max_threads() // number of threads

#include "common.h"
#include "math.h"
#include "synthData.h"
#include "p4p.h"
#include "r4p3v.h"
#include "r4p3v_ns.h"
#include "stats.h"



int main()
{
	Stats stats;

	char dataOut[100];
	sprintf_s(dataOut, "%sdata_%i_%i.txt", FOLDER_OUT, SCENE, MOTION);
	std::ofstream exportData(dataOut,std::ios::out);
	exportData.precision(DBL_DIG);

	for (int i=1; i<=NTRIALS; ++i)
	{
		double q[NVIEWS][NPOINTS][3];
		Camera cam_gt, cam_est;
		
		synthData(q,cam_gt); // generate synthetic image points and ground truth cameras
		
		stats.timer.start(); // start timer
		
		if (!r4p3v_ns(q,cam_est)) continue; // run solver

		stats.totalTime+=stats.timer.stop(); // stop timer

		stats.updateStats(cam_est,cam_gt); // update statistics

		//exportData << stats.errNum << "\n"; // export errors to file
	}
	exportData.close();

	stats.printStats(); // print results

	return 0;
}