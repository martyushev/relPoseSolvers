#include "quartic.h"



// triangulate scene points XQ_i using a pair of cameras (I 0) and Rt=(R t)
void triang2v(const double A[4][3][_P], const double Rt[12], const int &i, double XQ[3])
{
	double XQ3;
	const double t1=Rt[11]*A[1][1][i]-Rt[10]*A[1][2][i], t2=Rt[11]*A[1][0][i]-Rt[9]*A[1][2][i];
	if (fabs(t1)>fabs(t2))
		XQ3=t1/((Rt[3]*A[0][0][i]+Rt[4]*A[0][1][i]+Rt[5]*A[0][2][i])*A[1][2][i]-
			(Rt[6]*A[0][0][i]+Rt[7]*A[0][1][i]+Rt[8]*A[0][2][i])*A[1][1][i]);
	else
		XQ3=t2/((Rt[0]*A[0][0][i]+Rt[1]*A[0][1][i]+Rt[2]*A[0][2][i])*A[1][2][i]-
			(Rt[6]*A[0][0][i]+Rt[7]*A[0][1][i]+Rt[8]*A[0][2][i])*A[1][0][i]);

	XQ[0]=A[0][0][i]*XQ3;
	XQ[1]=A[0][1][i]*XQ3;
	XQ[2]=A[0][2][i]*XQ3;
}



// average rotation matrix
void aveRot(double Rt[12])
{
	// normalize the first column of Rt
	double fac=1./sqrt(Rt[0]*Rt[0]+Rt[3]*Rt[3]+Rt[6]*Rt[6]);
	Rt[0]*=fac;
	Rt[3]*=fac;
	Rt[6]*=fac;

	// make first two columns orthogonal
	fac=Rt[0]*Rt[1]+Rt[3]*Rt[4]+Rt[6]*Rt[7];
	Rt[1]-=Rt[0]*fac;
	Rt[4]-=Rt[3]*fac;
	Rt[7]-=Rt[6]*fac;

	// normalize the second column of Rt
	fac=1./sqrt(Rt[1]*Rt[1]+Rt[4]*Rt[4]+Rt[7]*Rt[7]);
	Rt[1]*=fac;
	Rt[4]*=fac;
	Rt[7]*=fac;

	// the last column is the cross product of the first two columns 
	Rt[2]=Rt[3]*Rt[7]-Rt[4]*Rt[6];
	Rt[5]=Rt[1]*Rt[6]-Rt[0]*Rt[7];
	Rt[8]=Rt[0]*Rt[4]-Rt[1]*Rt[3];
}



// output is either 1 or 0 (no solution found)
bool p4p(const double A[4][3][_P], const double XQ[_P][3], Camera &cam)
{
	double p3[4][12], d10[3], d20[3], d21[3];
	const double cosA[3]={2.*(A[2][0][1]*A[2][0][2]+A[2][1][1]*A[2][1][2]+A[2][2][1]*A[2][2][2]),
						  2.*(A[2][0][0]*A[2][0][2]+A[2][1][0]*A[2][1][2]+A[2][2][0]*A[2][2][2]),
						  2.*(A[2][0][0]*A[2][0][1]+A[2][1][0]*A[2][1][1]+A[2][2][0]*A[2][2][1])};

    for (int i=0; i<3; ++i)
	{
		d10[i]=XQ[1][i]-XQ[0][i];
		d20[i]=XQ[2][i]-XQ[0][i];
		d21[i]=XQ[2][i]-XQ[1][i];
	}
    double den=d10[0]*d20[1]-d10[1]*d20[0];
    if (!den) return 0;
	den=1./den;
    
	const double D[3]={d21[0]*d21[0]+d21[1]*d21[1]+d21[2]*d21[2], d20[0]*d20[0]+d20[1]*d20[1]+d20[2]*d20[2], d10[0]*d10[0]+d10[1]*d10[1]+d10[2]*d10[2]};
    const double c[4]={(d10[1]*d20[2]-d10[2]*d20[1])*den, (d10[2]*d20[0]-d10[0]*d20[2])*den, 1./(c[0]*c[0]+c[1]*c[1]+1.), 0.5*den};

    const double D1=1.0/D[1], t1=D[0]*D1, t2=D[2]*D1, t3=t1+t2;
	const double cosA2[4]={cosA[0]*cosA[0], cosA[1]*cosA[1], cosA[2]*cosA[2], cosA[0]*cosA[2]};
	const double t[3]={t1-t2, 1.-t[0], 1.+t[0]}, t6=t2*cosA2[0], t7=t1*cosA2[2], t8=t3*cosA2[3], t9=t8-cosA2[3], t02=2.*t[0];

	// coefficients of quartic polynomial pol
	double pol[5], v[4];
	pol[4]=t[1]*t[1]-t6;
	pol[3]=(t6+t02*t[1])*cosA[1]+t9;
	pol[2]=cosA2[0]-t6+cosA2[2]-t7-t8*cosA[1]+t[0]*t[0]*(2.0+cosA2[1])-2.0;
	pol[1]=(t7-t02*t[2])*cosA[1]+t9;
	pol[0]=t[2]*t[2]-t7;

	// real roots of quartic polynomial pol
	const int nv=solveQuartic(pol, v);
	int imin;

	cam.Err=1.;
	for (int i=0; i<nv; ++i)
	{
		if(v[i]<=0) continue;
		const double vcosA=v[i]*cosA[0];
        if (cosA[2]==vcosA) continue;
		const double u=(t[2]-(t[1]*v[i]+t[0]*cosA[1])*v[i])/(cosA[2]-vcosA);
		if (u<=0) continue;
        const double uv=v[i]*v[i]+u*(u-vcosA);
        if (uv<=0) continue;

		double Duv=D[0]/uv, fL0=sqrt(Duv), fL1=u*fL0, fL2=v[i]*fL0;
		const double s1=Duv+D[2]-fL1*fL1, s2=Duv+D[1]-fL2*fL2;
        const double dx=(s1*d20[1]-s2*d10[1])*c[3], dy=(s2*d10[0]-s1*d20[0])*c[3];
		const double t1=(c[0]*dx+c[1]*dy)*c[2], t0=(dx*dx+dy*dy-Duv)*c[2], t2=t1*t1-t0;
		if (t2<0) continue;

		// coordinates of 3rd camera center
		const double u1=(t1>0)? -t1-sqrt(t2):-t1+sqrt(t2);
		double O3[3]={c[0]*u1+dx+XQ[0][0], c[1]*u1+dy+XQ[0][1], u1+XQ[0][2]}, w0[3], w1[3], w2[3];

		// triplets (XP[0]-O3, XP[1]-O3, XP[2]-O3) and (data[][0], data[][1], data[][2]) must define the same orientation
		for (int k=0; k<3; ++k)
		{
			w0[k]=XQ[0][k]-O3[k];
			w1[k]=XQ[1][k]-O3[k];
			w2[k]=XQ[2][k]-O3[k];
		}

		const double det2=w2[2]*(w0[0]*w1[1]-w0[1]*w1[0])-w2[1]*(w0[0]*w1[2]-w0[2]*w1[0])+w2[0]*(w0[1]*w1[2]-w0[2]*w1[1]);
		if ((det2>0)-(A[3][0][3]>0))
		{
			const double u2=t0/u1;
			O3[0]=c[0]*u2+dx+XQ[0][0];
			O3[1]=c[1]*u2+dy+XQ[0][1];
			O3[2]=u2+XQ[0][2];
		}

		fL0=1./fL0;
		fL1=1./fL1;
		fL2=1./fL2;
        const double W[9]={(XQ[0][0]-O3[0])*fL0, (XQ[1][0]-O3[0])*fL1, (XQ[2][0]-O3[0])*fL2, (XQ[0][1]-O3[1])*fL0, (XQ[1][1]-O3[1])*fL1,
            (XQ[2][1]-O3[1])*fL2, (XQ[0][2]-O3[2])*fL0, (XQ[1][2]-O3[2])*fL1, (XQ[2][2]-O3[2])*fL2};

		p3[i][0]=A[3][0][0]*W[0]+A[3][0][1]*W[1]+A[3][0][2]*W[2];
		p3[i][3]=A[3][1][0]*W[0]+A[3][1][1]*W[1]+A[3][1][2]*W[2];
		p3[i][6]=A[3][2][0]*W[0]+A[3][2][1]*W[1]+A[3][2][2]*W[2];
		p3[i][1]=A[3][0][0]*W[3]+A[3][0][1]*W[4]+A[3][0][2]*W[5];
		p3[i][4]=A[3][1][0]*W[3]+A[3][1][1]*W[4]+A[3][1][2]*W[5];
		p3[i][7]=A[3][2][0]*W[3]+A[3][2][1]*W[4]+A[3][2][2]*W[5];
		p3[i][2]=A[3][0][0]*W[6]+A[3][0][1]*W[7]+A[3][0][2]*W[8];
		p3[i][5]=A[3][1][0]*W[6]+A[3][1][1]*W[7]+A[3][1][2]*W[8];
		p3[i][8]=A[3][2][0]*W[6]+A[3][2][1]*W[7]+A[3][2][2]*W[8];

		aveRot(p3[i]);

		p3[i][9]=-p3[i][0]*O3[0]-p3[i][1]*O3[1]-p3[i][2]*O3[2];
		p3[i][10]=-p3[i][3]*O3[0]-p3[i][4]*O3[1]-p3[i][5]*O3[2];
		p3[i][11]=-p3[i][6]*O3[0]-p3[i][7]*O3[1]-p3[i][8]*O3[2];
		
		// check that all scene points are in front of the third camera
		// if any of them is not, then solution is dropped
		for (int k=0; k<_P; ++k)
			if ((p3[i][6]*XQ[k][0]+p3[i][7]*XQ[k][1]+p3[i][8]*XQ[k][2]+p3[i][11])*A[2][2][k]<0) continue;
		
		// compute the squared reprojection error for the remaining scene point
		const double b2=1./(p3[i][6]*XQ[3][0]+p3[i][7]*XQ[3][1]+p3[i][8]*XQ[3][2]+p3[i][11]);
		const double b0=(p3[i][0]*XQ[3][0]+p3[i][1]*XQ[3][1]+p3[i][2]*XQ[3][2]+p3[i][9])*b2-A[2][0][3];
		const double b1=(p3[i][3]*XQ[3][0]+p3[i][4]*XQ[3][1]+p3[i][5]*XQ[3][2]+p3[i][10])*b2-A[2][1][3];
		const double err1=b0*b0+b1*b1;

		if (err1<cam.Err)
		{
			cam.Err=err1;
			imin=i;
		}
	}

	if (cam.Err==1.) return 0;
	for (int k=0; k<12; ++k) cam.Rt[1][k]=p3[imin][k];
	return 1;
}