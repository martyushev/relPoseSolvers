#include <iostream>
#include <cmath>
#include <cfloat>

#include "common.h"
#include "synthData.h"
#include "relpose4p3v.h"



int main()
{
	Hist hRt;
	double speed=0, err=0;
	Timer timer;
	hRt.iniHist(-15, 2);

	// ---------------------------------------------------------------------------------------------------------------------

	std::cout.precision(DBL_DIG);
	int ntrials=0;

	for (int i=1; i<=NTRIALS; ++i)
	{
		double data[_V][3][_P];
		Camera cam_gt, cam_est;

		synthData(data, cam_gt); // synthetic data

		timer.start(); // start timer

		if (!relative4p3v(data, cam_est)) continue; // run solver

		double T=timer.stop(); // stop timer
		speed+=T;

		double E=0; // numerical error
		for (int k=0; k<12; ++k)
		{
			const double t1=cam_gt.Rt[0][k]-cam_est.Rt[0][k], t2=cam_gt.Rt[1][k]-cam_est.Rt[1][k];
			E+=t1*t1+t2*t2;
		}
		E=0.5*log10(E);
		hRt.updateHist(E);

		++ntrials;
	}

	// results
	std::cout << "\nNumber of successful trials: " << ntrials << "\n\n";
	std::cout << "Average runtime: " << speed*((1e+6)/(double)ntrials) << " microseconds per call" << "\n\n";
	std::cout << "Median numerical error: " << pow(10, hRt.getMedian(0.5)) << "\n\n";

	return 0;
}