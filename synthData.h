struct Camera
{
	double Rt[NVIEWS][12];
	double Err;
};



// generate camera centers
void cameraCenters(double O[NVIEWS][3])
{
	O[0][0]=O[0][1]=O[0][2]=0; // 1st camera center is at the origin

	double theta=rnd(0,PI), phi=rnd(-PI,PI), b=1; // b is baseline length

	if (MOTION==1)
	{ // sideway motion
		O[1][0]=sin(phi)*b;
		O[1][1]=cos(phi)*b;
		O[1][2]=O[2][2]=0;
		for (int k=0; k<2; ++k) O[2][k]=0.5*O[1][k]+b*rnd(-0.1,0.1);
	}
	else if (MOTION==2)
	{ // forward motion
		O[1][0]=0;
		O[1][1]=0;
		O[1][2]=b;
		for (int k=0; k<3; ++k) O[2][k]=0.5*O[1][k];
	}
	else // (MOTION==0 || MOTION==3)
	{
		const double t=sin(theta)*b;
		O[1][0]=sin(phi)*t;
		O[1][1]=cos(phi)*t;
		O[1][2]=cos(theta)*b;
		for (int k=0; k<3; ++k) O[2][k]=0.5*O[1][k]+b*rnd(-0.1,0.1);
	}
}



// generate camera matrices
void cameraMatrices(const double O[NVIEWS][3], double Rt[NVIEWS][12])
{
	// the 1st camera matrix is [I 0]
	Rt[0][0]=Rt[0][4]=Rt[0][8]=1;
	Rt[0][1]=Rt[0][2]=Rt[0][3]=Rt[0][5]=Rt[0][6]=Rt[0][7]=Rt[0][9]=Rt[0][10]=Rt[0][11]=0;

	const double max_ang=(MOTION==3)? 0:PI/8.; // maximum value for Euler angles
	for (int j=1; j<NVIEWS; ++j)
	{
		// uniformly distributed Euler angles
		const double phi=rnd(-max_ang,max_ang), theta=rnd(0,2.*max_ang), psi=rnd(-max_ang,max_ang);
		const double cphi=cos(phi), sphi=sin(phi), ctheta=cos(theta), stheta=sin(theta), cpsi=cos(psi), spsi=sin(psi);

		// jth camera matrix
		const double R1[9]={cphi, sphi, 0, -sphi, cphi, 0, 0, 0, 1};
		const double R2[9]={1, 0, 0, 0, ctheta, stheta, 0, -stheta, ctheta};
		const double R3[9]={cpsi, spsi, 0, -spsi, cpsi, 0, 0, 0, 1};
		double R12[9];
		mult(R1,R2,R12);
		mult(R12,R3,Rt[j]);
		
		Rt[j][9]=-Rt[j][0]*O[j][0]-Rt[j][1]*O[j][1]-Rt[j][2]*O[j][2];
		Rt[j][10]=-Rt[j][3]*O[j][0]-Rt[j][4]*O[j][1]-Rt[j][5]*O[j][2];
		Rt[j][11]=-Rt[j][6]*O[j][0]-Rt[j][7]*O[j][1]-Rt[j][8]*O[j][2];
	}

	// normalize translation vectors so that ||t2||=1
	const double fac=1./sqrt(Rt[1][9]*Rt[1][9]+Rt[1][10]*Rt[1][10]+Rt[1][11]*Rt[1][11]);
	for (int j=1; j<NVIEWS; ++j)
	{
		Rt[j][9]*=fac;
		Rt[j][10]*=fac;
		Rt[j][11]*=fac;
	}
}



// generate scene points
void scenePoints(double Q[NPOINTS][3], const int &k)
{
	const double d=10, fov=PI/4., w=2.*d*tan(0.5*fov), h=w*0.8, depth=SCENE? 0:0.5*d;
	for (int i=0; i<NPOINTS; ++i)
	{
		Q[i][0]=0.5*rnd(-w,w);
		Q[i][1]=0.5*rnd(-h,h);
		Q[i][2]=rnd(d,d+depth);
	}

	if (k==0);
	// degenerate point configurations
	else if (k==1)
		{ // 3 points on a line
		const double a=rnd(0.2,0.8);
		Q[3][0]=Q[1][0]+a*(Q[2][0]-Q[1][0]);
		Q[3][1]=Q[1][1]+a*(Q[2][1]-Q[1][1]);
		Q[3][2]=Q[1][2]+a*(Q[2][2]-Q[1][2]);
	}
	else if (k==2)
	{ // points on a circle
		const double r=8;
		for (int i=0; i<NPOINTS; ++i)
		{
			double phi=rnd(0,2.*PI);
			Q[i][0]=r*cos(phi);
			Q[i][1]=r*sin(phi);
			Q[i][2]=d;
		}
	}
	else if (k==3)
	{ // points at vertices of a rectangle
		Q[0][0]=-w;
		Q[0][1]=-h;
		Q[0][2]=d;
		Q[1][0]=-w;
		Q[1][1]=h;
		Q[1][2]=d;
		Q[2][0]=w;
		Q[2][1]=-h;
		Q[2][2]=d;
		Q[3][0]=w;
		Q[3][1]=h;
		Q[3][2]=d;
	}
}



// map scene points onto image planes
void projectPoints(const double P[NVIEWS][12], const double Q[NPOINTS][3], double q[NVIEWS][NPOINTS][3])
{
	for (int i=0; i<NPOINTS; ++i)
	{
		const double fac=1./Q[i][2];
		q[0][i][0]=Q[i][0]*fac;
		q[0][i][1]=Q[i][1]*fac;
		q[0][i][2]=1;
		for (int j=1; j<NVIEWS; ++j)
		{
			double PjQ[3];
			PjQ[0]=P[j][0]*Q[i][0]+P[j][1]*Q[i][1]+P[j][2]*Q[i][2]+P[j][9];
			PjQ[1]=P[j][3]*Q[i][0]+P[j][4]*Q[i][1]+P[j][5]*Q[i][2]+P[j][10];
			PjQ[2]=P[j][6]*Q[i][0]+P[j][7]*Q[i][1]+P[j][8]*Q[i][2]+P[j][11];
			const double fac=1./PjQ[2];
			q[j][i][0]=PjQ[0]*fac;
			q[j][i][1]=PjQ[1]*fac;
			q[j][i][2]=1;
		}
	}
}



// add normally distributed image noise
void imageNoise(const double sigma, double q[NVIEWS][NPOINTS][3])
{
	for (int i=0; i<NPOINTS; ++i)
		for (int j=1; j<NVIEWS; ++j)
		{
			double noise[2];
			rndn(0.,sigma,noise); // generate normally distributed image noise
			q[j][i][0]+=noise[0];
			q[j][i][1]+=noise[1];
		}
}



// generate synthetic image points and ground truth camera matrices
void synthData(double q[NVIEWS][NPOINTS][3], Camera &cam)
{
	double O[NVIEWS][3], Q[NPOINTS][3], pixel=1e-4;
	cameraCenters(O);
	cameraMatrices(O,cam.Rt);
	scenePoints(Q,0);
	projectPoints(cam.Rt,Q,q);
	imageNoise(pixel*NOISE,q);
}