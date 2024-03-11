#include <windows.h>



struct Timer
{
	__int64 pf;
	__int64 timestarted;
	__int64 timefinished;
	double pfreq;

	Timer()
	{
		QueryPerformanceFrequency((LARGE_INTEGER *)&pf);
		pfreq=1./(double)pf;
	}

	void start()
	{
		QueryPerformanceCounter((LARGE_INTEGER *)&timestarted);
	}

	double stop()
	{
		QueryPerformanceCounter((LARGE_INTEGER *)&timefinished);
		return (timefinished-timestarted)*pfreq;
	}
};



struct Hist
{
	double u[NCOLS+1];
	int n[NCOLS];

	void iniHist(const double, const double);
	void updateHist(const double &);
	void printHist();
	double getMedian(const double);
	double getMean();
	double getFails(const double);
};



// initialize histogram, uN > u0
void Hist::iniHist(const double u0, const double uN)
{
	const double d=(uN-u0)/(double)NCOLS;
	u[0]=u0;
	u[NCOLS]=uN;
	n[0]=0;
	for (int i=1; i<NCOLS; ++i)
	{
		u[i]=u[i-1]+d;
		n[i]=0;
	}
}



// update histogram with new value x
void Hist::updateHist(const double &x)
{
	int m=NCOLS;
	if (x<u[0]) m=0;
	for (int i=0; i<NCOLS; ++i)
	{
		if (x>=u[i] && x<u[i+1])
		{
			m=i;
			break;
		}
	}
	if (x>=u[NCOLS]) m=NCOLS-1;
	++n[m];
}



// print histogram
void Hist::printHist()
{
	std::cout.precision(3);
	for (int i=0; i<NCOLS; ++i)
		std::cout << u[i] << "\t" << n[i] << "\n";
	std::cout << "\n";
}



// compute median (q=0.5) or lower quartile (q=0.25) on groupped data
double Hist::getMedian(const double q)
{
	int n_tot=n[0]; // total frequency
	for (int i=1; i<NCOLS; ++i) n_tot+=n[i];
	const double Q=q*(double)n_tot; // fraction of the total frequency
	
	int m; // find class median, which is the first class with the value of cumulative frequency equal at least Q
	double F=0.; // cumulative frequency of the class median
	for (int i=0; i<NCOLS; ++i)
	{
		F+=(double)n[i];
		if (F>=Q)
		{
			m=i;
			break;
		}
	}
	return u[m]+(Q-F+(double)n[m])/(double)n[m]*(u[m+1]-u[m]);
}



// compute mean on groupped data
double Hist::getMean()
{
	int n_tot=0; // total frequency
	double sum_fx=0.;
	for (int i=0; i<NCOLS; ++i)
	{
		n_tot+=n[i];
		sum_fx+=0.5*(u[i]+u[i+1])*(double)n[i];
	}
	return sum_fx/(double)n_tot;
}



// compute fails on groupped data
double Hist::getFails(const double q)
{
	int n_tot=0; // total frequency
	int n_fail=0;
	for (int i=0; i<NCOLS; ++i)
	{
		n_tot+=n[i];
		if (u[i]>q) n_fail+=n[i];
	}
	return 100*(double)n_fail/(double)n_tot;
}



// compute numerical error
double numError(const double P[2][12], const double Pgt[2][12])
{
	double err=0.;
	for (int k=0; k<12; ++k)
	{
		const double t1=P[0][k]-Pgt[0][k], t2=P[1][k]-Pgt[1][k];
		err+=t1*t1+t2*t2;
	}
	return 0.5*log10(err);
}



// compute rotational error
void rotError(const double P[2][12], const double Pgt[2][12], double err[2])
{
	for (int i=0; i<2; ++i)
	{
		double t=0;
		for (int j=0; j<9; ++j) t+=P[i][j]*Pgt[i][j];
		t=0.5*(t-1);
		err[i] = log10(acos(t));
	}
}



// compute translational error
void translError(const double P[2][12], const double Pgt[2][12], double err[2])
{
	for (int i=0; i<2; ++i)
	{
		const double fac1=Pgt[i][9]*Pgt[i][9]+Pgt[i][10]*Pgt[i][10]+Pgt[i][11]*Pgt[i][11];
		const double fac2=P[i][9]*P[i][9]+P[i][10]*P[i][10]+P[i][11]*P[i][11];
		const double t=(Pgt[i][9]*P[i][9]+Pgt[i][10]*P[i][10]+Pgt[i][11]*P[i][11])/sqrt(fac1*fac2);
		err[i]=log10(acos(t));
	}
}



struct Stats
{
	Timer timer;
	Hist hNum, hRot, hTra;
	int ntrials;
	double totalTime, errNum;

	Stats()
	{
		hNum.iniHist(-15,2);
		ntrials=0;
		totalTime=0.;
	}

	void updateStats(const Camera &, const Camera &);
	void printStats();
};



void Stats::updateStats(const Camera &cam_est, const Camera &cam_gt)
{
	errNum=numError(cam_est.Rt,cam_gt.Rt);
	hNum.updateHist(errNum); // update error histogram
	++ntrials;
}



void Stats::printStats()
{
	std::cout.precision(4);
	std::cout << "\nNumber of threads: " << NUM_THREADS << "\n\n";
	std::cout << "Number of successful trials: " << ntrials << "\n\n";
	std::cout << "Average runtime (ms): " << totalTime*((1e+3)/(double)ntrials) << "\n\n";
	std::cout << "Median numerical error: " << hNum.getMedian(0.5) << "\n\n";
	std::cout << "Mean numerical error: " << hNum.getMean() << "\n\n";
	std::cout << "Fails (%): " << hNum.getFails(-2.) << "\n\n";
	std::cout << "Numerical error distribution:\n"; hNum.printHist();
}