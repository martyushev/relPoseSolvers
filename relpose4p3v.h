#include "auxArrays.h"
#include "p4p.h"



// construct two rotation matrices
void getRot(const double ax[3], const double cs[2], double R[2][9])
{
	const double ax2[6]={ax[0]*ax[1], ax[1]*ax[2], ax[2]*ax[0], ax[0]*ax[0], ax[1]*ax[1], ax[2]*ax[2]};
	const double b[6]={cs[0]*ax2[0], cs[0]*ax2[1], cs[0]*ax2[2], cs[0]*ax2[3], cs[0]*ax2[4], cs[0]*ax2[5]};
	const double a[3]={cs[1]*ax[0], cs[1]*ax[1], cs[1]*ax[2]}, c[3]={cs[0]-b[3], cs[0]-b[4], cs[0]-b[5]};
	const double t[3]={ax2[0]-b[0], ax2[1]-b[1], ax2[2]-b[2]}, r[3]={ax2[0]+b[0], ax2[1]+b[1], ax2[2]+b[2]};
	R[0][0]=c[0]+ax2[3];
	R[0][1]=-a[2]+t[0];
	R[0][2]=a[1]+t[2];
	R[0][3]=a[2]+t[0];
	R[0][4]=c[1]+ax2[4];
	R[0][5]=-a[0]+t[1];
	R[0][6]=-a[1]+t[2];
	R[0][7]=a[0]+t[1];
	R[0][8]=c[2]+ax2[5];
	
	R[1][0]=-c[0]+ax2[3];
	R[1][1]=a[2]+r[0];
	R[1][2]=-a[1]+r[2];
	R[1][3]=-a[2]+r[0];
	R[1][4]=-c[1]+ax2[4];
	R[1][5]=a[0]+r[1];
	R[1][6]=a[1]+r[2];
	R[1][7]=-a[0]+r[1];
	R[1][8]=-c[2]+ax2[5];
}



// construct the quartic polynomial
void getPol(const auxArrays &S, const double u[2], double p[5], double d[5])
{
	double E1[6], G1[9];
	if (fabs(u[1])>fabs(u[0]))
	{
		const double f=u[0]/u[1];
		for (int i=0; i<6; ++i)
			E1[i]=S.E[0][i]*f+S.E[1][i];
		for (int i=0; i<9; ++i)
			G1[i]=(((S.G[0][i]*f+S.G[1][i])*f+S.G[2][i])*f+S.G[3][i])*f+S.G[4][i];
	}
	else
	{
		const double f=u[1]/u[0];
		for (int i=0; i<6; ++i)
			E1[i]=S.E[1][i]*f+S.E[0][i];
		for (int i=0; i<9; ++i)
			G1[i]=(((S.G[4][i]*f+S.G[3][i])*f+S.G[2][i])*f+S.G[1][i])*f+S.G[0][i];
	}
	
	// normalize conic G1
	double fac=0;
	for (int i=0; i<9; ++i) fac+=fabs(G1[i]);
	fac=1./fac;
	for (int i=0; i<9; ++i) G1[i]=G1[i]*fac;
	
	const double E02=2.0*E1[0], E12=2.0*E1[1], G02=2.0*G1[0], G13=G1[1]+G1[3], G26=G1[2]+G1[6], G57=G1[5]+G1[7], a=E02*(G13*E1[2]-E1[1]*G26);
	const double c[3]={E1[0]*(G13*E1[5]-E12*G1[8]), E02*(G13*E1[4]-E1[1]*G57), E1[0]*(G13*E1[3]-E12*G1[4])};
	
	d[0]=E1[0]*G1[8]-G1[0]*E1[5];
	d[1]=E1[0]*G57-G02*E1[4];
	d[2]=E1[0]*G1[4]-G1[0]*E1[3];
	d[3]=E1[0]*G26-G02*E1[2];
	d[4]=E1[0]*G13-E12*G1[0];
	
	const double phi=E1[0]*(E1[0]*G13-E12*G1[0]), ab=a*d[3], bb=d[3]*d[3], cc=d[4]*d[4], phi2=2.0*phi, t=phi2*d[0]-ab, d32=2.0*d[3];
	p[0]=(phi*d[0]-ab)*d[0]+bb*c[0];
	p[1]=d[4]*(d32*c[0]-a*d[0])+d[1]*t+bb*c[1];
	p[2]=d[2]*t+d[1]*(phi*d[1]-a*d[4])+bb*c[2]+d32*d[4]*c[1]+cc*c[0];
	p[3]=d[4]*(d32*c[2]-a*d[2])+phi2*d[1]*d[2]+cc*c[1];
	p[4]=phi*d[2]*d[2]+cc*c[2];
}



// check the cheirality constraint
bool cheirality(const double A[4][3][_P], double Rt[12], double XQ[_P][3])
{
	triang2v(A, Rt, 0, XQ[0]); // triangulate first scene point
	
	const double c1=XQ[0][2]*A[0][2][0], c2=(Rt[6]*XQ[0][0]+Rt[7]*XQ[0][1]+Rt[8]*XQ[0][2]+Rt[11])*A[1][2][0];
	if (c1>0 && c2>0); // rt2 is correct
	else if (c1<0 && c2<0) // r2 is correct, t2=-t2
	{
		Rt[9]=-Rt[9];
		Rt[10]=-Rt[10];
		Rt[11]=-Rt[11];
	}
	else // in this case r2 = H_t2 r2, where H_t2 = -I + 2 (t2 t2')/(t2' t2)
	{
		const double beta=2./(Rt[9]*Rt[9]+Rt[10]*Rt[10]+Rt[11]*Rt[11]);
		for (int j=0; j<3; ++j)
		{
			const int j3=j+3, j6=j+6;
			const double w=beta*(Rt[9]*Rt[j]+Rt[10]*Rt[j3]+Rt[11]*Rt[j6]);
			Rt[j]=Rt[9]*w-Rt[j];
			Rt[j3]=Rt[10]*w-Rt[j3];
			Rt[j6]=Rt[11]*w-Rt[j6];
		}
		triang2v(A, Rt, 0, XQ[0]);
		const double c1=XQ[0][2]*A[0][2][0], c2=(Rt[6]*XQ[0][0]+Rt[7]*XQ[0][1]+Rt[8]*XQ[0][2]+Rt[11])*A[1][2][0];
		if (c1>0 && c2>0);
		else
		{
			Rt[9]=-Rt[9];
			Rt[10]=-Rt[10];
			Rt[11]=-Rt[11];
		}
	}
	
    for (int i=0; i<_P; ++i)
		triang2v(A, Rt, i, XQ[i]); // triangulate all scene points
	
	// check that the rest of scene points are in front of the first two cameras
	// if any of them is not, then the solution is dropped
	for (int i=1; i<_P; ++i)
		if (XQ[i][2]*A[0][2][i]<0 || (Rt[6]*XQ[i][0]+Rt[7]*XQ[i][1]+Rt[8]*XQ[i][2]+Rt[11])*A[1][2][i]<0)
			return 0;
	return 1;
}



// output is either 1 or 0 (no solution found)
bool evalCameras(const double A[4][3][_P], const double e[2][3], Camera &cam)
{
	cam.Rt[0][9]=e[1][0];
	cam.Rt[0][10]=e[1][1];
	cam.Rt[0][11]=e[1][2];
	    
	const double ee=e[0][0]*e[1][0]+e[0][1]*e[1][1]+e[0][2]*e[1][2];
	const double cs[2]={ee, sqrt(1.0-ee*ee)};
	double ax[3]={-e[0][2]*e[1][1]+e[0][1]*e[1][2], e[0][2]*e[1][0]-e[0][0]*e[1][2], -e[0][1]*e[1][0]+e[0][0]*e[1][1]};
	const double fac=1./cs[1];
	ax[0]*=fac;
	ax[1]*=fac;
	ax[2]*=fac;
	
	double R[2][9], p[2][3], q[2][2][3], Y[2][9], eps[2];
	getRot(ax, cs, R);
	
	for (int i=0; i<2; ++i)
	{
		p[i][0]=A[1][1][i]*e[1][2]-A[1][2][i]*e[1][1];
		p[i][1]=A[1][2][i]*e[1][0]-A[1][0][i]*e[1][2];
		p[i][2]=A[1][0][i]*e[1][1]-A[1][1][i]*e[1][0];
		for (int k=0; k<2; ++k)
		{
			q[k][i][0]=R[k][0]*A[0][0][i]+R[k][1]*A[0][1][i]+R[k][2]*A[0][2][i];
			q[k][i][1]=R[k][3]*A[0][0][i]+R[k][4]*A[0][1][i]+R[k][5]*A[0][2][i];
			q[k][i][2]=R[k][6]*A[0][0][i]+R[k][7]*A[0][1][i]+R[k][8]*A[0][2][i];
		}
	}
	
	const double e12=e[1][2]*e[1][2], e10=e[1][0]*e[1][0], e11=e[1][1]*e[1][1];
	const double X[6]={1.-2.*(e12+e11), 2.*e[1][1]*e[1][0], 2.*e[1][2]*e[1][0], 1.-2.*(e12+e10), 2.*e[1][2]*e[1][1], 1.-2.*(e11+e10)};
	const double pX[3]={p[0][0]*X[0]+p[0][1]*X[1]+p[0][2]*X[2], p[0][0]*X[1]+p[0][1]*X[3]+p[0][2]*X[4], p[0][0]*X[2]+p[0][1]*X[4]+p[0][2]*X[5]};
	const double pe[3]={p[0][1]*e[1][2]-p[0][2]*e[1][1], p[0][2]*e[1][0]-p[0][0]*e[1][2], p[0][0]*e[1][1]-p[0][1]*e[1][0]};
	for (int i=0; i<2; ++i)
	{
		const double a=q[i][0][0]*pX[0]+q[i][0][1]*pX[1]+q[i][0][2]*pX[2], 
			b=q[i][0][0]*pe[0]+q[i][0][1]*pe[1]+q[i][0][2]*pe[2], 
			c=q[i][0][0]*p[0][0]+q[i][0][1]*p[0][1]+q[i][0][2]*p[0][2];
		
		if (fabs(a)>fabs(c))
		{
			const double a1=1.0/a, ba=b*a1, ca=c*a1, x=(ba>0.0)? -ba-sqrt(ba*ba-ca):-ba+sqrt(ba*ba-ca), x2=x*x, fac=1.0/(1.0+x2);
			const double s1=X[1]*x2, s2=2.0*e[1][2]*x, s3=X[2]*x2, s4=2.0*e[1][1]*x, s5=X[4]*x2, s6=2.0*e[1][0]*x;
			Y[i][0]=(X[0]*x2+1.0)*fac;
			Y[i][1]=(s1-s2)*fac;
			Y[i][2]=(s3+s4)*fac;
			Y[i][3]=(s1+s2)*fac;
			Y[i][4]=(X[3]*x2+1.0)*fac;
			Y[i][5]=(s5-s6)*fac;
			Y[i][6]=(s3-s4)*fac;
			Y[i][7]=(s5+s6)*fac;
			Y[i][8]=(X[5]*x2+1.0)*fac;
		}
		else
		{
			const double c1=1.0/c, bc=b*c1, ac=a*c1, x=(bc>0.0)? -bc-sqrt(bc*bc-ac):-bc+sqrt(bc*bc-ac), x2=x*x, fac=1./(1.+x2);
			Y[i][0]=(X[0]+x2)*fac;
			Y[i][1]=(X[1]-2.*e[1][2]*x)*fac;
			Y[i][2]=(X[2]+2.*e[1][1]*x)*fac;
			Y[i][3]=(X[1]+2.*e[1][2]*x)*fac;
			Y[i][4]=(X[3]+x2)*fac;
			Y[i][5]=(X[4]-2.*e[1][0]*x)*fac;
			Y[i][6]=(X[2]-2.*e[1][1]*x)*fac;
			Y[i][7]=(X[4]+2.*e[1][0]*x)*fac;
			Y[i][8]=(X[5]+x2)*fac;
		}
		eps[i]=fabs(q[i][1][0]*(p[1][0]*Y[i][0]+p[1][1]*Y[i][3]+p[1][2]*Y[i][6])+
		q[i][1][1]*(p[1][0]*Y[i][1]+p[1][1]*Y[i][4]+p[1][2]*Y[i][7])+
		q[i][1][2]*(p[1][0]*Y[i][2]+p[1][1]*Y[i][5]+p[1][2]*Y[i][8]));
	}
	
	const int ii=(eps[1]<eps[0]);
	for (int i=0; i<3; ++i)
	{
		const int i3=i+3, i6=i+6;
		cam.Rt[0][i]=Y[ii][0]*R[ii][i]+Y[ii][1]*R[ii][i3]+Y[ii][2]*R[ii][i6];
		cam.Rt[0][i3]=Y[ii][3]*R[ii][i]+Y[ii][4]*R[ii][i3]+Y[ii][5]*R[ii][i6];
		cam.Rt[0][i6]=Y[ii][6]*R[ii][i]+Y[ii][7]*R[ii][i3]+Y[ii][8]*R[ii][i6];
	}
	
	double XQ[_P][3];
	if (!cheirality(A, cam.Rt[0], XQ)) return 0;

	return p4p(A, XQ, cam); // use p4p algorithm to get 3rd camera matrix
}



// evaluate cost function at point theta
void costFunction(const double &theta, const auxArrays &S, Camera &cam)
{
	const double u[2]={cos(theta), sin(theta)};
	double del[5], p1[5], w[4], V[2][9];
	
	if (fabs(u[1])>fabs(u[0]))
	{
		const double f=u[0]/u[1];
		for (int i=0; i<9; ++i)
		{
			V[0][i]=S.V0[i];
			V[1][i]=((S.V1[0][i]*f+S.V1[1][i])*f+S.V1[2][i])*f+S.V1[3][i];
		}
	}
	else
	{
		const double f=u[1]/u[0];
		for (int i=0; i<9; ++i)
		{
			V[0][i]=S.V0[i];
			V[1][i]=((S.V1[3][i]*f+S.V1[2][i])*f+S.V1[1][i])*f+S.V1[0][i];
		}
	}

	getPol(S, u, p1, del);
	
	// find all real roots of degree-4 polynomial p1
	const int nw=solveQuartic(p1, w);
	int kmin;
	
	Camera cam1[4];
	cam.Err=1.;
	for (int k=0; k<nw; ++k)
	{
		const double v=-((del[2]*w[k]+del[1])*w[k]+del[0])/(del[4]*w[k]+del[3]);
		double e[2][3]; // epipoles e[1]=t2, e[1]=R2.e[0]
		for (int j=0; j<2; ++j)
		{
			e[j][0]=V[j][0]*v+V[j][1]*w[k]+V[j][2];
			e[j][1]=V[j][3]*v+V[j][4]*w[k]+V[j][5];
			e[j][2]=V[j][6]*v+V[j][7]*w[k]+V[j][8];
			const double norm=1./sqrt(e[j][0]*e[j][0]+e[j][1]*e[j][1]+e[j][2]*e[j][2]);
			e[j][0]*=norm;
			e[j][1]*=norm;
			e[j][2]*=norm;
		}
		
		if (!evalCameras(S.A, e, cam1[k])) continue;
		
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
		costFunction(x2, S, cam2);
	}
	else
	{
		x2=lm[1];
		x1=x0+RATIO*(x2-x0);
		cam2=cam;
		costFunction(x1, S, cam1);
	}
	for (int k=0; k<MAXIT && fabs(x3-x0)>tol*(fabs(x1)+fabs(x2)); ++k)
	{
		if (cam1.Err<cam2.Err)
		{
			x3=x2;
			x2=x1;
			x1=x0+RATIO*(x2-x0);
			cam2=cam1;
			costFunction(x1, S, cam1);
		}
		else
		{
			x0=x1;
			x1=x2;
			x2=x3-RATIO*(x3-x1);
			cam1=cam2;
			costFunction(x2, S, cam2);
		}
	}
	if (cam1.Err<cam2.Err) cam=cam1; else cam=cam2;
}



// find local minima of the cost function
int localMinima(const auxArrays &S, Camera cam[MAXLM])
{
	double theta[2*MAXLM+1], lm[MAXLM][3];
	const double N1=PI/(double)(2*MAXLM);
	Camera cam1[2*MAXLM+1];
	
	theta[0]=0.0;
	for (int i=1; i<2*MAXLM+1; ++i)
		theta[i]=theta[i-1]+N1;
	
	#pragma omp parallel for shared(theta,S,cam1) schedule(dynamic) num_threads(NUM_THREADS)
	for (int i=0; i<2*MAXLM-1; ++i) // evaluate cost function at each theta[i]
		costFunction(theta[i], S, cam1[i]);
	
	cam1[2*MAXLM-1]=cam1[0];
	cam1[2*MAXLM]=cam1[1];
	
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
	for (int i=0; i<n; ++i)
		golden(lm[i], S, cam[i]); // polish each local minimum by golden section method

	return n;
}



// main function
// output is either 1 or 0 (no solution found)
bool relative4p3v(const double A0[_V][3][_P], Camera &cam_est)
{
	auxArrays S;
	S.iniPerm(A0);
	S.getA();
	S.getEG();

	Camera cam[MAXLM];
	const int n=localMinima(S, cam);
	if (!n) return 0; // no local minimum found
	
	int imin=0;
	double Epsmin=cam[0].Err;
	for (int i=1; i<n; ++i) // find the global minimum of the cost function
		if (cam[i].Err<Epsmin)
		{
			Epsmin=cam[i].Err;
			imin=i;
		}
	cam_est=cam[imin];

	return 1;
}