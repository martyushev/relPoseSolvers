// p4p algorithm
// output is either 1 or 0 (no solution found)
bool p4p(const double q[3], const double K[18], const double Q[NPOINTS][3], Camera &cam)
{
	double d10[3], d20[3], d21[3];
	for (int i=0; i<3; ++i)
	{
		d10[i]=Q[1][i]-Q[0][i];
		d20[i]=Q[2][i]-Q[0][i];
		d21[i]=Q[2][i]-Q[1][i];
	}
	const double fac=1./(d10[0]*d20[1]-d10[1]*d20[0]);
	const double D[3]={d21[0]*d21[0]+d21[1]*d21[1]+d21[2]*d21[2], d20[0]*d20[0]+d20[1]*d20[1]+d20[2]*d20[2], d10[0]*d10[0]+d10[1]*d10[1]+d10[2]*d10[2]};
	const double c0=(d10[1]*d20[2]-d10[2]*d20[1])*fac, c1=(d10[2]*d20[0]-d10[0]*d20[2])*fac, c2=1./(c0*c0+c1*c1+1.), c3=0.5*fac;

	const double t0=1./D[1], t1=D[0]*t0, t2=D[2]*t0, t3=t1+t2;
	const double r0=t1-t2, r1=1.-r0, r2=1.+r0, t5=2.*r0, t6=t2*K[12], t7=t1*K[14], t8=t3*K[15], t9=t8-K[15];

	// coefficients of quartic polynomial
	double p[5], v[4];
	p[0]=r2*r2-t7;
	p[1]=(t7-t5*r2)*K[10]+t9;
	p[2]=K[16]-t6-t7-t8*K[10]+r0*r0*(2.+K[13])-2.;
	p[3]=(t6+t5*r1)*K[10]+t9;
	p[4]=r1*r1-t6;
	
	const int nv=solveQuartic(p,v); // find all real roots of the polynomial
	
	int imin;
	double P3[4][12];
	cam.Err=1;
	for (int i=0; i<nv; ++i)
	{
		if (v[i]<=0) continue;
		const double vcosA=v[i]*K[9];
		if (K[11]==vcosA) continue;
		const double u=(r2-(r1*v[i]+r0*K[10])*v[i])/(K[11]-vcosA);
		if (u<=0) continue;
		const double uv=v[i]*v[i]+u*(u-vcosA);
		if (uv<=0) continue;
		
		double Duv=D[0]/uv, fL0=sqrt(Duv), fL1=u*fL0, fL2=v[i]*fL0;
		const double s1=Duv+D[2]-fL1*fL1, s2=Duv+D[1]-fL2*fL2;
		const double dx=(s1*d20[1]-s2*d10[1])*c3, dy=(s2*d10[0]-s1*d20[0])*c3;
		const double t0=(dx*dx+dy*dy-Duv)*c2, t1=(c0*dx+c1*dy)*c2, t2=t1*t1-t0;
		if (t2<0) continue;
		
		// coordinates of the 3rd camera center
		const double u1=(t1>0)? -t1-sqrt(t2):-t1+sqrt(t2);
		double O3[3]={c0*u1+dx+Q[0][0], c1*u1+dy+Q[0][1], u1+Q[0][2]}, w[3][3];
		for (int j=0; j<3; ++j) for (int k=0; k<3; ++k) w[j][k]=Q[j][k]-O3[k];

		const double det2=w[2][2]*(w[0][0]*w[1][1]-w[0][1]*w[1][0])-w[2][1]*(w[0][0]*w[1][2]-w[0][2]*w[1][0])+w[2][0]*(w[0][1]*w[1][2]-w[0][2]*w[1][1]);
		if ((det2>0)-(K[17]>0))
		{ // triplets (Q[0]-O3, Q[1]-O3, Q[2]-O3) and (q[2][0], q[2][1], q[2][2]) must define the same orientation
			const double u2=t0/u1;
			O3[0]=c0*u2+dx+Q[0][0];
			O3[1]=c1*u2+dy+Q[0][1];
			O3[2]=u2+Q[0][2];
		}

		fL0=1./fL0;
		fL1=1./fL1;
		fL2=1./fL2;
		const double W[9]={(Q[0][0]-O3[0])*fL0, (Q[0][1]-O3[1])*fL0, (Q[0][2]-O3[2])*fL0, (Q[1][0]-O3[0])*fL1, (Q[1][1]-O3[1])*fL1, (Q[1][2]-O3[2])*fL1, (Q[2][0]-O3[0])*fL2, (Q[2][1]-O3[1])*fL2, (Q[2][2]-O3[2])*fL2};
		mult(K,W,P3[i]);
		P3[i][9]=-P3[i][0]*O3[0]-P3[i][1]*O3[1]-P3[i][2]*O3[2];
		P3[i][10]=-P3[i][3]*O3[0]-P3[i][4]*O3[1]-P3[i][5]*O3[2];
		P3[i][11]=-P3[i][6]*O3[0]-P3[i][7]*O3[1]-P3[i][8]*O3[2];
		
		// compute the squared reprojection error for the 4th scene point
		const double b2=1./(P3[i][6]*Q[3][0]+P3[i][7]*Q[3][1]+P3[i][8]*Q[3][2]+P3[i][11]);
		const double b0=(P3[i][0]*Q[3][0]+P3[i][1]*Q[3][1]+P3[i][2]*Q[3][2]+P3[i][9])*b2-q[0];
		const double b1=(P3[i][3]*Q[3][0]+P3[i][4]*Q[3][1]+P3[i][5]*Q[3][2]+P3[i][10])*b2-q[1];
		
		const double err=b0*b0+b1*b1;
		if (err<cam.Err)
		{
			cam.Err=err;
			imin=i;
		}
	}
	if (cam.Err==1) return 0;

	for (int k=0; k<12; ++k) cam.Rt[2][k]=P3[imin][k];

	return 1;
}