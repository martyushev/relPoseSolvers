#define PI 3.1415926535897932385
#define RATIO 0.6180339887498948482

#include "solvePoly.h"



// uniformly distributed random numbers from [a, b]
double rnd(const double &a, const double &b)
{
	return (b-a)*rand()/(double)RAND_MAX+a;
}



// normally distributed random numbers with mathematical expectation=a, stdandard deviation=sigma, Box-Muller method
void rndn(const double a, const double sigma, double x[2])
{
	double s0, s1;
	s0=rnd(1.e-15,1.);
	s1=rnd(1.e-15,2.*PI);
	s0=sigma*sqrt(-2.*log(s0));
	x[0]=s0*cos(s1)+a;
	x[1]=s0*sin(s1)+a;
}



// multiply two 3-by-3 matrices C=A*B
void mult(const double A[9], const double B[9], double C[9])
{
	C[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	C[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	C[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	C[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	C[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	C[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	C[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	C[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	C[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}



// find the basis of the null space of M-by-N (M < N) matrix A using specifically tailored QR factorization
// result is N-M N-vectors Q
// note that matrix A changes during the computation
template <const int M, const int N>
void nullQR(double A[M][N], double Q[N-M][N])
{
	const int NM=N-M, M1=M-1;

	// construct Hauseholder vectors
	for (int j=0; j<M; ++j)
	{
		double t=0;
		for (int i=j+1; i<N; ++i) t+=A[j][i]*A[j][i];
		double mu=sqrt(A[j][j]*A[j][j]+t);
		mu=(A[j][j]<0)? 1./(A[j][j]-mu):1./(A[j][j]+mu);

		for (int i=j+1; i<N; ++i) A[j][i]*=mu;

		double beta=-2./(1.+t*(mu*mu));
		for (int k=j+1; k<M; ++k)
		{
			double w=A[k][j];
			for (int i=j+1; i<N; ++i) w+=A[k][i]*A[j][i];
			w*=beta;
			for (int i=j+1; i<N; ++i) A[k][i]+=A[j][i]*w;
		}
	}

	// multiply M Householder matrices, inverse order is more efficient
	// we only need last N-M columns of the resulting matrix
	for (int i=0; i<NM; ++i) memset(Q[i],0,sizeof(double)*M1);
	
	// start from the Mth matrix
	double beta=1;
	for (int i=M; i<N; ++i) beta+=A[M1][i]*A[M1][i];

	beta=-2./beta;
	for (int k=0; k<NM; ++k)
	{
		Q[k][M1]=A[M1][k+M]*beta;
		for (int i=M; i<N; ++i)
			Q[k][i]=(k==i-M)? 1.+A[M1][i]*Q[k][M1]:A[M1][i]*Q[k][M1];
	}
	
	// multiply by the remaining M-1 matrices
	for (int j=M1-1; j>=0; --j)
	{
		double beta=1;
		for (int i=j+1; i<N; ++i) beta+=A[j][i]*A[j][i];
		beta=-2./beta;
		for (int k=0; k<NM; ++k)
		{
			double w=Q[k][j];
			for (int i=j+1; i<N; ++i) w+=Q[k][i]*A[j][i];
			w*=beta;
			Q[k][j]+=w;
			for (int i=j+1; i<N; ++i) Q[k][i]+=A[j][i]*w;
		}
	}
}



// compute adjoint to 3-by-3 matrix A, i.e. det(A) A^{-1}
void adjoint(const double A[9], double B[9], const bool sym)
{
	B[0]=A[4]*A[8]-A[5]*A[7];
	B[1]=-A[1]*A[8]+A[2]*A[7];
	B[2]=A[1]*A[5]-A[2]*A[4];
	B[4]=A[0]*A[8]-A[2]*A[6];
	B[5]=-A[0]*A[5]+A[2]*A[3];
	B[8]=A[0]*A[4]-A[1]*A[3];
	if (sym)
	{
		B[3]=B[1];
		B[6]=B[2];
		B[7]=B[5];
		
	}
	else
	{
		B[3]=-A[3]*A[8]+A[5]*A[6];
		B[6]=A[3]*A[7]-A[4]*A[6];
		B[7]=-A[0]*A[7]+A[1]*A[6];
	}
}



// transpose of 3-by-3 matrix A
void transpose(const double A[9], double B[9])
{
	B[0]=A[0];
	B[1]=A[3];
	B[2]=A[6];
	B[3]=A[1];
	B[4]=A[4];
	B[5]=A[7];
	B[6]=A[2];
	B[7]=A[5];
	B[8]=A[8];
}



// triangulate scene point Q_i using the pair of cameras [I 0] and Rt = [R t]
void triang2v(const double q[NVIEWS][NPOINTS][3], const double Rt[12], const int &i, double Q[NPOINTS][3])
{
	const double t1=(Rt[3]*q[0][i][0]+Rt[4]*q[0][i][1]+Rt[5]*q[0][i][2])*q[1][i][2]-(Rt[6]*q[0][i][0]+Rt[7]*q[0][i][1]+Rt[8]*q[0][i][2])*q[1][i][1];
	const double t2=(Rt[0]*q[0][i][0]+Rt[1]*q[0][i][1]+Rt[2]*q[0][i][2])*q[1][i][2]-(Rt[6]*q[0][i][0]+Rt[7]*q[0][i][1]+Rt[8]*q[0][i][2])*q[1][i][0];
	double Q3;
	if (fabs(t1)>fabs(t2))
	{
		const double t=Rt[11]*q[1][i][1]-Rt[10]*q[1][i][2];
		Q3=t/t1;
	}
	else
	{
		const double t=Rt[11]*q[1][i][0]-Rt[9]*q[1][i][2];
		Q3=t/t2;
	}

	Q[i][0]=q[0][i][0]*Q3;
	Q[i][1]=q[0][i][1]*Q3;
	Q[i][2]=q[0][i][2]*Q3;
}



// check the cheirality constraint
bool cheirality(const double q[NVIEWS][NPOINTS][3], double Rt[12], double Q[NPOINTS][3])
{
	triang2v(q,Rt,0,Q); // triangulate 1st scene point
	
	const double c1=Q[0][2], c2=Rt[6]*Q[0][0]+Rt[7]*Q[0][1]+Rt[8]*Q[0][2]+Rt[11];
	if (c1>0 && c2>0); // Rt2 is correct
	else if (c1<0 && c2<0) // R2 is correct, t2=-t2
	{
		Rt[9]=-Rt[9];
		Rt[10]=-Rt[10];
		Rt[11]=-Rt[11];
	}
	else
	{ // in this case R2 = Ht*R2, where Ht = -I + 2*(t2*t2')/(t2'*t2)
		const double b=2./(Rt[9]*Rt[9]+Rt[10]*Rt[10]+Rt[11]*Rt[11]);
		for (int j=0; j<3; ++j)
		{
			const int j3=j+3, j6=j+6;
			const double w=b*(Rt[9]*Rt[j]+Rt[10]*Rt[j3]+Rt[11]*Rt[j6]);
			Rt[j]=Rt[9]*w-Rt[j];
			Rt[j3]=Rt[10]*w-Rt[j3];
			Rt[j6]=Rt[11]*w-Rt[j6];
		}
		triang2v(q,Rt,0,Q);
		const double c1=Q[0][2], c2=Rt[6]*Q[0][0]+Rt[7]*Q[0][1]+Rt[8]*Q[0][2]+Rt[11];
		if (c1>0 && c2>0); // Rt2 is correct
		else // R2 is correct, t2=-t2
		{
			Rt[9]=-Rt[9];
			Rt[10]=-Rt[10];
			Rt[11]=-Rt[11];
		}
	}
	
	for (int i=0; i<NPOINTS; ++i)
		triang2v(q,Rt,i,Q); // triangulate all scene points
	
	// check that the rest of scene points are in front of the first two cameras
	// if any of them is not, then the solution is dropped
	for (int i=1; i<NPOINTS; ++i)
		if (Q[i][2]<0 || Rt[6]*Q[i][0]+Rt[7]*Q[i][1]+Rt[8]*Q[i][2]+Rt[11]<0) return 0;

	return 1;
}



double getVol(const double q[NPOINTS][3], const int ind[NPOINTS])
{
	const double m1=q[ind[0]][1]*q[ind[1]][2]-q[ind[0]][2]*q[ind[1]][1];
	const double m2=q[ind[1]][0]*q[ind[0]][2]-q[ind[0]][0]*q[ind[1]][2];
	const double m3=q[ind[0]][0]*q[ind[1]][1]-q[ind[1]][0]*q[ind[0]][1];
	return fabs(q[ind[2]][0]*m1+q[ind[2]][1]*m2+q[ind[2]][2]*m3);
}



// permute image points to improve numerics
void iniPerm(const double q0[NVIEWS][NPOINTS][3], double q[NVIEWS][NPOINTS][3])
{
	int m=0, ind[4][NPOINTS]={{2,3,1,0}, {2,0,3,1}, {0,1,3,2}, {2,1,0,3}};
	double maxVol=getVol(q0[0],ind[0]);
	for (int i=1; i<4; ++i)
	{
		const double vol=getVol(q0[0],ind[i]);
		if (vol<maxVol)
		{
			m=i;
			maxVol=vol;
		}
	}

	for (int j=0; j<NVIEWS; ++j)
		for (int i=0; i<NPOINTS; ++i) memcpy(q[j][i],q0[j][ind[m][i]],sizeof(double)*3);
}



void getK(double q[NVIEWS][NPOINTS][3], double K[18])
{
	for (int j=0; j<2; ++j)
		for (int i=0; i<NPOINTS; ++i)
		{
			const double fac=1./q[j][i][2];
			for (int k=0; k<3; ++k) q[j][i][k]*=fac;
		}
	double fac=1./q[2][3][2];
	for (int k=0; k<3; ++k) q[2][3][k]*=fac;
	
	for (int i=0; i<3; ++i)
	{
		const double fac=1./sqrt(q[2][i][0]*q[2][i][0]+q[2][i][1]*q[2][i][1]+q[2][i][2]*q[2][i][2]);
		for (int k=0; k<3; ++k) q[2][i][k]*=fac;
	}
	
	const double m1=q[2][0][1]*q[2][1][2]-q[2][0][2]*q[2][1][1];
	const double m2=q[2][1][0]*q[2][0][2]-q[2][0][0]*q[2][1][2];
	const double m3=q[2][0][0]*q[2][1][1]-q[2][1][0]*q[2][0][1];

	K[17]=1./(q[2][2][0]*m1+q[2][2][1]*m2+q[2][2][2]*m3);
	// inverse of matrix
	K[0]=(q[2][1][1]*q[2][2][2]-q[2][1][2]*q[2][2][1])*K[17];
	K[1]=(q[2][0][2]*q[2][2][1]-q[2][0][1]*q[2][2][2])*K[17];
	K[2]=m1*K[17];
	K[3]=(q[2][2][0]*q[2][1][2]-q[2][1][0]*q[2][2][2])*K[17];
	K[4]=(q[2][0][0]*q[2][2][2]-q[2][2][0]*q[2][0][2])*K[17];
	K[5]=m2*K[17];
	K[6]=(q[2][1][0]*q[2][2][1]-q[2][2][0]*q[2][1][1])*K[17];
	K[7]=(q[2][2][0]*q[2][0][1]-q[2][0][0]*q[2][2][1])*K[17];
	K[8]=m3*K[17];

	K[9]=2.*(q[2][1][0]*q[2][2][0]+q[2][1][1]*q[2][2][1]+q[2][1][2]*q[2][2][2]);
	K[10]=2.*(q[2][0][0]*q[2][2][0]+q[2][0][1]*q[2][2][1]+q[2][0][2]*q[2][2][2]);
	K[11]=2.*(q[2][0][0]*q[2][1][0]+q[2][0][1]*q[2][1][1]+q[2][0][2]*q[2][1][2]);
	K[12]=K[9]*K[9];
	K[13]=K[10]*K[10];
	K[14]=K[11]*K[11];
	K[15]=K[9]*K[11];
	K[16]=K[12]+K[14];
}