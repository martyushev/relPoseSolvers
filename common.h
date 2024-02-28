#include <windows.h>

#define PI 3.1415926535897932385
#define RATIO 0.6180339887498948482
#define _V 3 // number of views
#define _V1 (_V-1)
#define _P 4 // number of points
#define MAXIT 30 // number of iterations for golden section search
#define MAXLM 50 // max number of local minima (> 1), number of cost function evaluations is 2*MAXLM-1
#define NCOLS 100 // number of columns in histogram
#define NTRIALS 100000 // number of trials


struct Camera
{
	double Err;
	double Rt[_V1][12];
};


struct Timer
{
	__int64 pf;
	__int64 timestarted;
	__int64 timefinished;
	double pfreq;

	Timer()
	{
		QueryPerformanceFrequency((LARGE_INTEGER *)&pf);
		pfreq=1.0/(double)pf;
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
	for (int i=1; i<NCOLS; ++i)
		n_tot+=n[i];
	const double Q=q*(double)n_tot; // fraction of the total frequency

    int m; // find class median, which is the first class with the value of cumulative frequency equal at least Q
    double F=0.0; // cumulative frequency of the class median
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
	double sum_fx=0.0;
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
	double err=0;
	for (int k=0; k<12; ++k)
	{
		const double t1=P[0][k]-Pgt[0][k], t2=P[1][k]-Pgt[1][k];
		err+=t1*t1+t2*t2;
	}
	return 0.5*log10(err);
}