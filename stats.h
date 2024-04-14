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
	double u[NBINS+1];
	int n[NBINS];

	void iniHist(const double &, const double &);
	void updateHist(const double &);
	void printHist();
	double getMedian(const double &);
	double getMean();
	double getFails(const double &);
};



// initialize histogram, uN > u0
void Hist::iniHist(const double &u0, const double &uN)
{
	const double d=(uN-u0)/(double)NBINS;
	u[0]=u0;
	u[NBINS]=uN;
	n[0]=0;
	for (int i=1; i<NBINS; ++i)
	{
		u[i]=u[i-1]+d;
		n[i]=0;
	}
}



// update histogram with new value x
void Hist::updateHist(const double &x)
{
	int m=NBINS;
	if (x<u[0]) m=0;
	for (int i=0; i<NBINS; ++i)
	{
		if (x>=u[i] && x<u[i+1])
		{
			m=i;
			break;
		}
	}
	if (x>=u[NBINS]) m=NBINS-1;
	++n[m];
}



// print histogram
void Hist::printHist()
{
	std::cout.precision(3);
	for (int i=0; i<NBINS; ++i)
		std::cout << u[i] << "\t" << n[i] << "\n";
	std::cout << "\n";
}



// compute median (q=0.5) or lower quartile (q=0.25) on groupped data
double Hist::getMedian(const double &q)
{
	int n_tot=n[0]; // total frequency
	for (int i=1; i<NBINS; ++i) n_tot+=n[i];
	const double Q=q*(double)n_tot; // fraction of the total frequency
	
	int m; // find class median, which is the first class with the value of cumulative frequency equal at least Q
	double F=0; // cumulative frequency of the class median
	for (int i=0; i<NBINS; ++i)
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
	double sum_fx=0;
	for (int i=0; i<NBINS; ++i)
	{
		n_tot+=n[i];
		sum_fx+=0.5*(u[i]+u[i+1])*(double)n[i];
	}
	return sum_fx/(double)n_tot;
}



// compute fails on groupped data
double Hist::getFails(const double &q)
{
	int n_tot=0; // total frequency
	int n_fail=0;
	for (int i=0; i<NBINS; ++i)
	{
		n_tot+=n[i];
		if (u[i]>q) n_fail+=n[i];
	}
	n_fail+=NTRIALS-n_tot;
	return 100.*(double)n_fail/(double)NTRIALS;
}



// compute numerical error
double numError(const double P[NVIEWS][12], const double Pgt[NVIEWS][12])
{
	double err=0;
	for (int k=0; k<12; ++k)
	{
		const double t1=P[1][k]-Pgt[1][k], t2=P[2][k]-Pgt[2][k];
		err+=t1*t1+t2*t2;
	}
	return 0.5*log10(err);
}



// compute rotational errors
void rotError(const double P[NVIEWS][12], const double Pgt[NVIEWS][12], double err[NVIEWS])
{
	for (int i=1; i<NVIEWS; ++i)
	{
		double Rt[9];
		transpose(P[i],Rt);
		double R[9];
		mult(Rt,Pgt[i],R);
		const double cosR=0.5*(R[0]+R[4]+R[8]-1.);
		const double v[3]={R[1]*R[5]-R[2]*(R[4]-1.), R[2]*R[3]+(1.-R[0])*R[5], (R[0]-1.)*(R[4]-1.)-R[3]*R[1]}; // rotation axis
		const double fac=1./sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		const double sinR=0.5*fac*(v[0]*(R[7]-R[5])+v[1]*(R[2]-R[6])+v[2]*(R[3]-R[1]));
		err[i]=log10(atan2(fabs(sinR),cosR));
	}
}



// compute translational errors
void traError(const double P[NVIEWS][12], const double Pgt[NVIEWS][12], double err[NVIEWS])
{
	for (int i=1; i<NVIEWS; ++i)
	{
		const double fac1=Pgt[i][9]*Pgt[i][9]+Pgt[i][10]*Pgt[i][10]+Pgt[i][11]*Pgt[i][11];
		const double fac2=P[i][9]*P[i][9]+P[i][10]*P[i][10]+P[i][11]*P[i][11];
		const double cosT=(Pgt[i][9]*P[i][9]+Pgt[i][10]*P[i][10]+Pgt[i][11]*P[i][11])/sqrt(fac1*fac2);
		const double x[3]={-P[i][10]*Pgt[i][11]+P[i][11]*Pgt[i][10], P[i][9]*Pgt[i][11]-P[i][11]*Pgt[i][9], -P[i][9]*Pgt[i][10]+P[i][10]*Pgt[i][9]};
		const double sinT=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
		const double t1=atan2(sinT,cosT), t2=atan2(sinT,-cosT);
		err[i]=(t1<t2)? t1:t2;
		err[i]=log10(err[i]);
	}
}



struct Stats
{
	Timer timer;
	Hist hNum, hRot[NVIEWS], hTra[NVIEWS];
	int ntrials;
	double totalTime, errNum, errRot[NVIEWS], errTra[NVIEWS];

	Stats()
	{
		hNum.iniHist(-15,2);
		for (int i=1; i<NVIEWS; ++i)
		{
			hRot[i].iniHist(-15,2);
			hTra[i].iniHist(-15,2);
		}
		ntrials=0;
		totalTime=0;
	}

	void updateStats(const Camera &, const Camera &);
	void printStats();
};



void Stats::updateStats(const Camera &cam_est, const Camera &cam_gt)
{
	errNum=numError(cam_est.Rt,cam_gt.Rt);
	hNum.updateHist(errNum);
	rotError(cam_est.Rt,cam_gt.Rt,errRot);
	traError(cam_est.Rt,cam_gt.Rt,errTra);
	for (int i=1; i<NVIEWS; ++i)
	{
		hRot[i].updateHist(errRot[i]);
		hTra[i].updateHist(errTra[i]);
	}
	++ntrials;
}



void Stats::printStats()
{
	std::cout.precision(4);
	std::cout << "\nNumber of threads: " << NUM_THREADS << "\n\n";
	std::cout << "Number of successful trials: " << ntrials << "\n\n";
	std::cout << "Fails (%): " << hNum.getFails(-2) << "\n\n";
	std::cout << "Average runtime (ms): " << totalTime*((1e+3)/(double)ntrials) << "\n\n";
	std::cout << "Median numerical error: " << hNum.getMedian(0.5) << "\n\n";
	std::cout << "Average numerical error: " << hNum.getMean() << "\n\n";
	std::cout << "Median rotational error for 2nd camera (degrees): " << (180./PI)*pow(10.,hRot[1].getMedian(0.5)) << "\n\n";
	std::cout << "Median translational error for 2nd camera (degrees): " << (180./PI)*pow(10.,hTra[1].getMedian(0.5)) << "\n\n";
	std::cout << "Median rotational error for 3rd camera (degrees): " << (180./PI)*pow(10.,hRot[2].getMedian(0.5)) << "\n\n";
	std::cout << "Median translational error for 3rd camera (degrees): " << (180./PI)*pow(10.,hTra[2].getMedian(0.5)) << "\n\n";
	std::cout << "Numerical error distribution:\n"; hNum.printHist();
}