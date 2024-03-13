#include "auxArrays.h"
#include "p4p.h"



void svdE(const double E[9], double U[9], double Vt[9])
{
	const double E0[3]={E[1]*E[5]-E[2]*E[4], E[2]*E[3]-E[0]*E[5], E[0]*E[4]-E[1]*E[3]},
		E1[3]={E[1]*E[8]-E[2]*E[7], E[2]*E[6]-E[0]*E[8], E[0]*E[7]-E[1]*E[6]},
		E2[3]={E[4]*E[8]-E[5]*E[7], E[5]*E[6]-E[3]*E[8], E[3]*E[7]-E[4]*E[6]};

	double fac0=1./sqrt(E0[0]*E0[0]+E0[1]*E0[1]+E0[2]*E0[2]),
		fac1=1./sqrt(E1[0]*E1[0]+E1[1]*E1[1]+E1[2]*E1[2]),
		fac2=1./sqrt(E2[0]*E2[0]+E2[1]*E2[1]+E2[2]*E2[2]);

	if (fac0<=fac1 && fac0<=fac2)
		for(int i=6; i<9; ++i) Vt[i]=E0[i-6]*fac0;
	else if (fac1<=fac0 && fac1<=fac2)
		for(int i=6; i<9; ++i) Vt[i]=E1[i-6]*fac1;
	else
		for(int i=6; i<9; ++i) Vt[i]=E2[i-6]*fac2;

	fac0=1./sqrt(E[0]*E[0]+E[1]*E[1]+E[2]*E[2]);
	Vt[0]=E[0]*fac0;
	Vt[1]=E[1]*fac0;
	Vt[2]=E[2]*fac0;
	Vt[3]=Vt[7]*Vt[2]-Vt[8]*Vt[1];
	Vt[4]=Vt[8]*Vt[0]-Vt[6]*Vt[2];
	Vt[5]=Vt[6]*Vt[1]-Vt[7]*Vt[0];

	U[0]=E[0]*Vt[0]+E[1]*Vt[1]+E[2]*Vt[2];
	U[3]=E[3]*Vt[0]+E[4]*Vt[1]+E[5]*Vt[2];
	U[6]=E[6]*Vt[0]+E[7]*Vt[1]+E[8]*Vt[2];
	U[1]=E[0]*Vt[3]+E[1]*Vt[4]+E[2]*Vt[5];
	U[4]=E[3]*Vt[3]+E[4]*Vt[4]+E[5]*Vt[5];
	U[7]=E[6]*Vt[3]+E[7]*Vt[4]+E[8]*Vt[5];
		
	fac0=1./sqrt(U[0]*U[0]+U[3]*U[3]+U[6]*U[6]);
	fac1=1./sqrt(U[1]*U[1]+U[4]*U[4]+U[7]*U[7]);
	U[0]*=fac0;
	U[3]*=fac0;
	U[6]*=fac0;
	U[1]*=fac1;
	U[4]*=fac1;
	U[7]*=fac1;
	U[2]=U[3]*U[7]-U[6]*U[4];
	U[5]=U[6]*U[1]-U[0]*U[7];
	U[8]=U[0]*U[4]-U[3]*U[1];
}



// find all candidates for the epipoles
int getEpipoles(const auxArrays &S, const double u[2], double B1[6], double e[4][2][3])
{
	double V1[9], G1[6];
	if (fabs(u[1])>fabs(u[0]))
	{
		const double r=u[0]/u[1];
		for (int i=0; i<6; ++i)
		{
			B1[i]=S.B[0][i]*r+S.B[1][i];
			V1[i]=((S.V[0][i]*r+S.V[1][i])*r+S.V[2][i])*r+S.V[3][i];
			G1[i]=(((S.G[0][i]*r+S.G[1][i])*r+S.G[2][i])*r+S.G[3][i])*r+S.G[4][i];
		}
		for (int i=6; i<9; ++i) V1[i]=((S.V[0][i]*r+S.V[1][i])*r+S.V[2][i])*r+S.V[3][i];
	}
	else
	{
		const double r=u[1]/u[0];
		for (int i=0; i<6; ++i)
		{
			B1[i]=S.B[1][i]*r+S.B[0][i];
			V1[i]=((S.V[3][i]*r+S.V[2][i])*r+S.V[1][i])*r+S.V[0][i];
			G1[i]=(((S.G[4][i]*r+S.G[3][i])*r+S.G[2][i])*r+S.G[1][i])*r+S.G[0][i];
		}
		for (int i=6; i<9; ++i) V1[i]=((S.V[3][i]*r+S.V[2][i])*r+S.V[1][i])*r+S.V[0][i];
	}

	const double f=1./(B1[0]*G1[1]-B1[1]*G1[0]), fd=2.*f, fh=0.5*f;
    const double a=(G1[1]*B1[2]-B1[1]*G1[2])*fd;
    const double c=(B1[0]*G1[2]-G1[0]*B1[2])*f, c2=c*c, ac=a*c, cd=2.*c;
	const double g[3]={(G1[1]*B1[3]-B1[1]*G1[3])*f, (G1[1]*B1[4]-B1[1]*G1[4])*fd, (G1[1]*B1[5]-B1[1]*G1[5])*f};
	const double d[3]={(B1[0]*G1[3]-G1[0]*B1[3])*fh, (B1[0]*G1[4]-G1[0]*B1[4])*f, (B1[0]*G1[5]-G1[0]*B1[5])*fh};

	double p[5], v[4];
	p[4]=d[0]*d[0]+g[0];
    p[3]=d[0]*(-a+2.*d[1])+cd*g[0]+g[1];
    p[2]=d[0]*(-ac+2.*d[2])+cd*g[1]+d[1]*(d[1]-a)+c2*g[0]+g[2];
    p[1]=d[1]*(-ac+2.*d[2])+cd*g[2]-d[2]*a+c2*g[1];
    p[0]=d[2]*(-ac+d[2])+c2*g[2];

	const int nv=solveQuartic(p,v); // find all real roots of the quartic polynomial

	for (int k=0; k<nv; ++k)
	{
		const double v1=v[k], u1=-(d[0]*v1*v1+d[1]*v1+d[2])/(v1+c);
        e[k][0][0]=u1;
		e[k][0][1]=v1;
		e[k][0][2]=1.;
		e[k][1][0]=V1[0]*e[k][0][0]+V1[1]*e[k][0][1]+V1[2];
		e[k][1][1]=V1[3]*e[k][0][0]+V1[4]*e[k][0][1]+V1[5];
		e[k][1][2]=V1[6]*e[k][0][0]+V1[7]*e[k][0][1]+V1[8];
		for (int j=0; j<2; ++j)
		{
			const double fac=1./sqrt(e[k][j][0]*e[k][j][0]+e[k][j][1]*e[k][j][1]+e[k][j][2]*e[k][j][2]);
			e[k][j][0]*=fac;
			e[k][j][1]*=fac;
			e[k][j][2]*=fac;
		}
	}

	return nv;
}



// triangulate scene point Q_i using the pair of cameras [I 0] and Rt = [R t]
void triang2v(const double A[NVIEWS+1][3][NPOINTS], const double Rt[12], const int &i, double Q[NPOINTS][3])
{
	const double t1=Rt[11]*A[1][1][i]-Rt[10]*A[1][2][i];
	const double t2=Rt[11]*A[1][0][i]-Rt[9]*A[1][2][i];
	const double t3=((Rt[3]*A[0][0][i]+Rt[4]*A[0][1][i]+Rt[5]*A[0][2][i])*A[1][2][i]-(Rt[6]*A[0][0][i]+Rt[7]*A[0][1][i]+Rt[8]*A[0][2][i])*A[1][1][i]);
	const double t4=((Rt[0]*A[0][0][i]+Rt[1]*A[0][1][i]+Rt[2]*A[0][2][i])*A[1][2][i]-(Rt[6]*A[0][0][i]+Rt[7]*A[0][1][i]+Rt[8]*A[0][2][i])*A[1][0][i]);
	const double Q3=(fabs(t3)>fabs(t4))? t1/t3:t2/t4;

	Q[i][0]=A[0][0][i]*Q3;
	Q[i][1]=A[0][1][i]*Q3;
	Q[i][2]=A[0][2][i]*Q3;
}



// check the cheirality constraint
bool cheirality(const double A[NVIEWS+1][3][NPOINTS], double Rt[12], double Q[NPOINTS][3])
{
	triang2v(A,Rt,0,Q); // triangulate first scene point
	
	const double c1=Q[0][2], c2=Rt[6]*Q[0][0]+Rt[7]*Q[0][1]+Rt[8]*Q[0][2]+Rt[11];
	if (c1>0 && c2>0); // Rt2 is correct
	else if (c1<0 && c2<0) // R2 is correct, t2=-t2
	{
		Rt[9]=-Rt[9];
		Rt[10]=-Rt[10];
		Rt[11]=-Rt[11];
	}
	else // in this case R2 = Ht R2, where Ht = -I + 2(t2 t2')/(t2't2)
	{
		const double b=2./(Rt[9]*Rt[9]+Rt[10]*Rt[10]+Rt[11]*Rt[11]);
		for (int j=0; j<3; ++j)
		{
			const int j3=j+3, j6=j+6;
			const double w=b*(Rt[9]*Rt[j]+Rt[10]*Rt[j3]+Rt[11]*Rt[j6]);
			Rt[j]=Rt[9]*w-Rt[j];
			Rt[j3]=Rt[10]*w-Rt[j3];
			Rt[j6]=Rt[11]*w-Rt[j6];
		}
		triang2v(A,Rt,0,Q);
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
		triang2v(A,Rt,i,Q); // triangulate all scene points
	
	// check that the rest of scene points are in front of the first two cameras
	// if any of them is not, then the solution is dropped
	for (int i=1; i<NPOINTS; ++i)
		if (Q[i][2]<0 || Rt[6]*Q[i][0]+Rt[7]*Q[i][1]+Rt[8]*Q[i][2]+Rt[11]<0) return 0;

	return 1;
}



// output is either 1 or 0 (no solution found)
bool getCameras(const double A[NVIEWS+1][3][NPOINTS], const double Ht[9], const double B1[6], const double e[2][3], Camera &cam)
{
	const double l0[3]={B1[0]*e[0][0]+B1[1]*e[0][1]+B1[2]*e[0][2], B1[1]*e[0][0]+B1[3]*e[0][1]+B1[4]*e[0][2], B1[2]*e[0][0]+B1[4]*e[0][1]+B1[5]*e[0][2]};
	const double l1[3]={B1[0]*e[1][0]+B1[1]*e[1][1]+B1[2]*e[1][2], B1[1]*e[1][0]+B1[3]*e[1][1]+B1[4]*e[1][2], B1[2]*e[1][0]+B1[4]*e[1][1]+B1[5]*e[1][2]};
	const double xa[3]={l0[1]*l1[2]-l0[2]*l1[1], l0[2]*l1[0]-l0[0]*l1[2], l0[0]*l1[1]-l0[1]*l1[0]};
	const double u[3]={0.,xa[2],-xa[1]}; // u is arbitrary nonzero vector
	const double a=(u[1]*xa[2]-u[2]*xa[1])*e[0][0]+(u[2]*xa[0]-u[0]*xa[2])*e[0][1]+(u[0]*xa[1]-u[1]*xa[0])*e[0][2];
	const double b=l0[0]*u[0]+l0[1]*u[1]+l0[2]*u[2];
	const double F[9]={a*B1[0], a*B1[1]+b*xa[2], a*B1[2]-b*xa[1], a*B1[1]-b*xa[2], a*B1[3], a*B1[4]+b*xa[0], a*B1[2]+b*xa[1], a*B1[4]-b*xa[0], a*B1[5]};
	double E[9];
	mult(Ht,F,E); // essential matrix E=H'*F

	double U[9], Vt[9];
	svdE(E,U,Vt);

	cam.Rt[0][0]=U[0]*Vt[3]-U[1]*Vt[0]+U[2]*Vt[6];
	cam.Rt[0][1]=U[0]*Vt[4]-U[1]*Vt[1]+U[2]*Vt[7];
	cam.Rt[0][2]=U[0]*Vt[5]-U[1]*Vt[2]+U[2]*Vt[8];
	cam.Rt[0][3]=U[3]*Vt[3]-U[4]*Vt[0]+U[5]*Vt[6];
	cam.Rt[0][4]=U[3]*Vt[4]-U[4]*Vt[1]+U[5]*Vt[7];
	cam.Rt[0][5]=U[3]*Vt[5]-U[4]*Vt[2]+U[5]*Vt[8];
	cam.Rt[0][6]=U[6]*Vt[3]-U[7]*Vt[0]+U[8]*Vt[6];
	cam.Rt[0][7]=U[6]*Vt[4]-U[7]*Vt[1]+U[8]*Vt[7];
	cam.Rt[0][8]=U[6]*Vt[5]-U[7]*Vt[2]+U[8]*Vt[8];
	cam.Rt[0][9]=U[2];
	cam.Rt[0][10]=U[5];
	cam.Rt[0][11]=U[8];

	double Q[NPOINTS][3];
	if (!cheirality(A,cam.Rt[0],Q)) return 0;

	return p4p(A,Q,cam); // use p4p algorithm to get 3rd camera matrix
}



// evaluate cost function at point theta
void costFunction(const double &theta, const auxArrays &S, Camera &cam)
{
	const double u[2]={sin(theta), cos(theta)};
	double B1[6], e[4][2][3];
	const int nv=getEpipoles(S,u,B1,e);
	
	int kmin;
	Camera cam1[4];
	cam.Err=1.;
	for (int k=0; k<nv; ++k)
	{
		if (!getCameras(S.A,S.Ht,B1,e[k],cam1[k])) continue;
		if (cam1[k].Err<cam.Err)
		{
			cam.Err=cam1[k].Err;
			kmin=k;
		}
	}
	if (cam.Err==1.) return;

	cam=cam1[kmin];
}



// method of golden section
void golden(const double lm[3], const auxArrays &S, Camera &cam)
{
	const double tol=1.e-12;
	double x0=lm[0], x1, x2, x3=lm[2];
	Camera cam1, cam2;
	
	if (x3-lm[1]>lm[1]-x0)
	{
		x1=lm[1];
		x2=x3-RATIO*(x3-x1);
		cam1=cam;
		costFunction(x2,S,cam2);
	}
	else
	{
		x2=lm[1];
		x1=x0+RATIO*(x2-x0);
		cam2=cam;
		costFunction(x1,S,cam1);
	}
	for (int k=0; k<MAXIT && fabs(x3-x0)>tol*(fabs(x1)+fabs(x2)); ++k)
	{
		if (cam1.Err<cam2.Err)
		{
			x3=x2;
			x2=x1;
			x1=x0+RATIO*(x2-x0);
			cam2=cam1;
			costFunction(x1,S,cam1);
		}
		else
		{
			x0=x1;
			x1=x2;
			x2=x3-RATIO*(x3-x1);
			cam1=cam2;
			costFunction(x2,S,cam2);
		}
	}
	if (cam1.Err<cam2.Err) cam=cam1; else cam=cam2;
}



// find local minima of the cost function
int localMinima(const auxArrays &S, Camera cam[MAXLM])
{
	double theta[2*MAXLM+1], lm[MAXLM][3];
	const double N=PI/(double)(2*MAXLM);
	Camera cam1[2*MAXLM+1];
	
	theta[0]=0.;
	for (int i=1; i<2*MAXLM+1; ++i) theta[i]=theta[i-1]+N;
	
	#pragma omp parallel for shared(theta,S,cam1) schedule(dynamic) num_threads(NUM_THREADS)
	for (int i=0; i<2*MAXLM; ++i) costFunction(theta[i],S,cam1[i]);
	
	cam1[2*MAXLM]=cam1[0];
	
	int n=0;
	for (int i=1; i<2*MAXLM; ++i) // find local minima
	{
		const int im1=i-1, ip1=i+1;
		if (cam1[i].Err<cam1[im1].Err && cam1[i].Err<cam1[ip1].Err)
		{
			lm[n][0]=theta[im1];
			lm[n][1]=theta[i];
			lm[n][2]=theta[ip1];
			cam[n++]=cam1[i];
		}
	}

	#pragma omp parallel for shared(n,lm,S,cam) schedule(dynamic) num_threads(NUM_THREADS)
	for (int i=0; i<n; ++i) golden(lm[i],S,cam[i]); // polish each local minimum

	return n;
}



// main function
// output is either 1 or 0 (no solution found)
bool relative4p3v(const double data[NVIEWS][3][NPOINTS], Camera &cam_est)
{
	auxArrays S;
	S.getA(data);
	S.getBG();

	Camera cam[MAXLM];
	const int n=localMinima(S,cam);
	if (!n) return 0; // no local minimum found
	
	int imin=0;
	double minErr=cam[0].Err;
	for (int i=1; i<n; ++i) // find the global minimum of the cost function
		if (cam[i].Err<minErr)
		{
			minErr=cam[i].Err;
			imin=i;
		}
	cam_est=cam[imin];

	return 1;
}