struct Camera
{
	double Rt[NVIEWS1][12];
	double Err;
};



// generate camera centers
void cameraCenters(double O[NVIEWS1][3])
{
	// cameras' centers, 1st camera center is at the origin
	double theta=rnd(0,PI), phi=rnd(-PI,PI), b=1.; // b is baseline length

	if (!MOTION)
	{
		const double t=sin(theta)*b;
		O[0][0]=sin(phi)*t;
		O[0][1]=cos(phi)*t;
		O[0][2]=cos(theta)*b;
	}
	else if (MOTION==1) // sideway motion
	{
		O[0][0]=sin(phi)*b;
		O[0][1]=cos(phi)*b;
		O[0][2]=0.;
	}
	else // forward motion
	{
		O[0][0]=0.;
		O[0][1]=0.;
		O[0][2]=b;
	}

	for (int k=0; k<3; ++k) O[1][k]=0.5*O[0][k]+(MOTION? 0.:b*rnd(-0.1,0.1));
}



// generate camera matrices
void cameraMatrices(const double O[NVIEWS1][3], double Rt[NVIEWS1][12])
{
	for (int j=0; j<NVIEWS1; ++j)
	{
		// uniformly distributed Euler angles
		const double max_ang=PI/8.; // maximum value for Euler angles
		const double phi=rnd(-max_ang,max_ang), theta=rnd(0.,2.*max_ang), psi=rnd(-max_ang,max_ang);
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
	const double scale=1./sqrt(Rt[0][9]*Rt[0][9]+Rt[0][10]*Rt[0][10]+Rt[0][11]*Rt[0][11]);
	for (int j=0; j<NVIEWS1; ++j)
	{
		Rt[j][9]*=scale;
		Rt[j][10]*=scale;
		Rt[j][11]*=scale;
	}
}



// generate a scene point
void scenePoints(double Q[NPOINTS][3])
{
	const double d=10, fov=PI/4., w=2.*d*tan(0.5*fov), h=w*0.8, depth=SCENE? 0.:d*rnd(0.,1.);
	for (int i=0; i<NPOINTS; ++i)
	{
		Q[i][0]=0.5*rnd(-w,w);
		Q[i][1]=0.5*rnd(-h,h);
		Q[i][2]=rnd(d,d+depth);
	}
}



// map scene points onto image planes
void projectPoints(const double P[NVIEWS1][12], const double Q[NPOINTS][3], double data[NVIEWS][3][NPOINTS])
{
	for (int i=0; i<NPOINTS; ++i)
	{
		const double fac=1./Q[i][2];
		data[0][0][i]=Q[i][0]*fac;
		data[0][1][i]=Q[i][1]*fac;
		data[0][2][i]=1.;
		for (int j=0; j<NVIEWS1; ++j)
		{
			double PjQ[3];
			PjQ[0]=P[j][0]*Q[i][0]+P[j][1]*Q[i][1]+P[j][2]*Q[i][2]+P[j][9];
			PjQ[1]=P[j][3]*Q[i][0]+P[j][4]*Q[i][1]+P[j][5]*Q[i][2]+P[j][10];
			PjQ[2]=P[j][6]*Q[i][0]+P[j][7]*Q[i][1]+P[j][8]*Q[i][2]+P[j][11];
			const double fac=1./PjQ[2];
			const int j1=j+1;
			data[j1][0][i]=PjQ[0]*fac;
			data[j1][1][i]=PjQ[1]*fac;
			data[j1][2][i]=1.;
		}
	}
}



// add normally distributed image noise
void imageNoise(const double sigma, double data[NVIEWS][3][NPOINTS])
{
	for (int i=0; i<NPOINTS; ++i)
		for (int j=1; j<NVIEWS; ++j)
		{
			double noise[2];
			rndn(0.,sigma,noise); // generate normally distributed image noise
			data[j][0][i]+=noise[0];
			data[j][1][i]+=noise[1];
		}
}



// generate synthetic data and ground truth camera matrices
void synthData(double data[NVIEWS][3][NPOINTS], Camera &cam)
{
	double O[NVIEWS1][3], Q[NPOINTS][3], pixel=1.e-4;
	cameraCenters(O);
	cameraMatrices(O,cam.Rt);
	scenePoints(Q);
	projectPoints(cam.Rt,Q,data);
	imageNoise(pixel*NOISE, data);
}