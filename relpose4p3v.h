#include "auxArrays.h"



// evaluate polynomial p at s, Horner scheme
// output is sextic polynomial p1 in one variable w
inline void get_p1(const double p[22], const double &s, double p1[7])
{
	p1[0]=s*s*((((p[4]*s+p[3])*s+p[2])*s+p[1])*s+p[0]);
	p1[1]=s*(((((p[10]*s+p[9])*s+p[8])*s+p[7])*s+p[6])*s+p[5]);
	p1[2]=(((((p[17]*s+p[16])*s+p[15])*s+p[14])*s+p[13])*s+p[12])*s+p[11];
	p1[3]=(((((-p[18]*s+p[19])*s-p[20])*s+p[21])*s+p[20])*s+p[19])*s+p[18];
	p1[4]=(((((p[11]*s-p[12])*s+p[13])*s-p[14])*s+p[15])*s-p[16])*s+p[17];
	p1[5]=((((p[5]*s-p[6])*s+p[7])*s-p[8])*s+p[9])*s-p[10];
	p1[6]=(((p[0]*s-p[1])*s+p[2])*s-p[3])*s+p[4];
}



// get value of u as a null-vector of matrix F1
double get_u(const double F[9][9], const double &s, const double &w)
{
	double F1[6];
	F1[0]=((F[0][1]*s+F[0][2]*w+F[0][3])*s+F[0][4]*w+F[0][5])*s+F[0][6]*w;
	F1[1]=((F[1][1]*w+F[1][3])*s+(F[1][2]*w+F[1][4])*w+F[1][6])*s+(F[1][5]*w+F[1][7])*w;
	F1[2]=((F[2][2]*w+F[2][1]*s+F[2][4])*w+F[2][3]*s+F[2][6])*w+F[2][5]*s;
	F1[3]=((F[3][1]*s+F[3][2]*w+F[3][3])*s+F[3][4]*w+F[3][5])*s+F[3][6]*w;
	F1[4]=((F[4][1]*w+F[4][3])*s+(F[4][2]*w+F[4][4])*w+F[4][6])*s+(F[4][5]*w+F[4][7])*w;
	F1[5]=((F[5][2]*w+F[5][1]*s+F[5][4])*w+F[5][3]*s+F[5][6])*w+F[5][5]*s;

	return (F1[0]*F1[5]-F1[3]*F1[2])/(F1[3]*F1[1]-F1[0]*F1[4]);
}



// output is either 1 or 0 (no solution found)
bool getCameras(const double A[NVIEWS+1][3][NPOINTS], const double B[9][10], const double &u, const double &s, const double &w, Camera &cam)
{
	// get rotation matrix R by Cayley's formula
	const double u2=u*u, v=u*s, v2=v*v, w2=w*w, uv=u*v, uw=u*w, vw=v*w;
	const double fac=1./(1.+u2+v2+w2), fac2=2.*fac, tu=1.-u2, tvw=v2-w2;
	cam.Rt[0][0]=(1.+u2-v2-w2)*fac;
	cam.Rt[0][1]=(uv+w)*fac2;
	cam.Rt[0][2]=(uw-v)*fac2;
	cam.Rt[0][3]=(uv-w)*fac2;
	cam.Rt[0][4]=(1.-u2+v2-w2)*fac;
	cam.Rt[0][5]=(vw+u)*fac2;
	cam.Rt[0][6]=(uw+v)*fac2;
	cam.Rt[0][7]=(vw-u)*fac2;
	cam.Rt[0][8]=(1.-u2-v2+w2)*fac;

	const double B1[3]={B[0][0]*u2+B[0][3]*v2+B[0][4]*vw+B[0][5]*w2+B[0][6]*u+B[0][9], B[1][1]*uv+B[1][2]*uw+B[1][7]*v+B[1][8]*w, B[2][1]*uv+B[2][2]*uw+B[2][7]*v+B[2][8]*w};

	// compute translation vector t2, formulas follow from the epipolar constraint for Q1
	cam.Rt[0][9]=B1[2]*cam.Rt[0][2];
	cam.Rt[0][10]=B1[2]*cam.Rt[0][5];
	cam.Rt[0][11]=-B1[0]*cam.Rt[0][2]-B1[1]*cam.Rt[0][5];
	
	// normalize translation vector so that ||t2|| = 1
	const double fac1=1./sqrt(cam.Rt[0][9]*cam.Rt[0][9]+cam.Rt[0][10]*cam.Rt[0][10]+cam.Rt[0][11]*cam.Rt[0][11]);
	cam.Rt[0][9]*=fac1;
	cam.Rt[0][10]*=fac1;
	cam.Rt[0][11]*=fac1;

	double Q[NPOINTS][3];
	if (!cheirality(A,cam.Rt[0],Q)) return 0;

	return p4p(A,Q,cam); // use p4p algorithm to get 3rd camera matrix
}



// evaluate cost function at point s
void costFunction(const double &s, const auxArrays &S, Camera &cam)
{
	double p1[7], w[6];
	get_p1(S.p,s,p1);
	const int nw=solvePoly(p1,w); // find all real roots of the sextic polynomial
	
	int kmin;
	Camera cam1[6];
	cam.Err=1.;
	for (int k=0; k<nw; ++k)
	{
		const double u=get_u(S.F,s,w[k]);
		if (!getCameras(S.A,S.B,u,s,w[k],cam1[k])) continue;
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
// this function adapted from "Numerical recipes in C" by Press et al.
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
	double s[2*MAXLM+1], lm[MAXLM][3];
	const double N=1./(double)MAXLM;
	Camera cam1[2*MAXLM+1];
	
	s[0]=-1.;
	for (int i=1; i<MAXLM; ++i)
	{
		s[i]=s[i-1]+N;
		s[2*MAXLM-1-i]=-s[i];
	}
	s[2*MAXLM-1]=1.;
	s[2*MAXLM]=-1./s[1];
	
	#pragma omp parallel for shared(s,S,cam1) schedule(dynamic) num_threads(NUM_THREADS)
	for (int i=0; i<2*MAXLM-1; ++i) costFunction(s[i],S,cam1[i]);

	cam1[2*MAXLM-1]=cam1[0];
	cam1[2*MAXLM]=cam1[1];
	
	int n=0;
	for (int i=1; i<2*MAXLM; ++i)
	{ // find local minima
		const int im1=i-1, ip1=i+1;
		if (cam1[i].Err<cam1[im1].Err && cam1[i].Err<cam1[ip1].Err)
		{
			lm[n][0]=s[im1];
			lm[n][1]=s[i];
			lm[n][2]=s[ip1];
			cam[n++]=cam1[i];
		}
	}
	
	#pragma omp parallel for shared(n,lm,S,cam) schedule(dynamic) num_threads(NUM_THREADS)
	for (int i=0; i<n; ++i) golden(lm[i],S,cam[i]); // polish each local minimum by golden section method

	return n;
}



// main function
// output is either 1 or 0 (no solution found)
bool relpose4p3v(const double data[NVIEWS][3][NPOINTS], Camera &cam_est)
{
	auxArrays S;
	S.trans(data);
	getA(S.A);
	S.getB();
	S.getF();
	S.decicPoly();

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

	S.itrans(cam_est.Rt);
	return 1;
}