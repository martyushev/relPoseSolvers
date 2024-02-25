// uniformly distributed random numbers from [a, b]
double rnd(const double a, const double b)
{
	return (b-a)*rand()/(double)RAND_MAX+a;
}



// multiply two 3x3 matrices C=A.B
void matrix_mult(const double A[9], const double B[9], double C[9])
{
	for (int i=0; i<3; ++i)
	{
		const int i3=3*i;
		for (int k=0; k<3; ++k)
			C[i3+k]=A[i3]*B[k]+A[i3+1]*B[k+3]+A[i3+2]*B[k+6];
	}
}



// generate camera centers
void cameraCenters(double XO[_V1][3])
{
	// cameras' centers, 1st camera center is at the origin
	double theta=rnd(0,PI), phi=rnd(-PI,PI), b=0.1; // b is baseline length

	const double t=sin(theta)*b;
	XO[0][0]=sin(phi)*t;
	XO[0][1]=cos(phi)*t;
	XO[0][2]=cos(theta)*b;

	for (int k=0; k<3; ++k) XO[1][k]=0.5*XO[0][k]+rnd(-0.01,0.01);
}



// generate camera matrices
void cameraMatrices(const double XO[_V1][3], double Rt[_V1][12])
{
	for (int j=0; j<_V1; ++j)
	{
		// uniformly distributed Euler angles
		const double max_ang=PI/8.; // maximal value for Euler angles
		const double phi=rnd(-max_ang, max_ang), theta=rnd(0, 2.*max_ang), psi=rnd(-max_ang, max_ang);
		const double cphi=cos(phi), sphi=sin(phi), ctheta=cos(theta), stheta=sin(theta), cpsi=cos(psi), spsi=sin(psi);

		// compute jth camera matrix
		const double R1[9]={cphi, sphi, 0, -sphi, cphi, 0, 0, 0, 1};
		const double R2[9]={1, 0, 0, 0, ctheta, stheta, 0, -stheta, ctheta};
		const double R3[9]={cpsi, spsi, 0, -spsi, cpsi, 0, 0, 0, 1};
		double R12[9];
		matrix_mult(R1, R2, R12);
		matrix_mult(R12, R3, Rt[j]);
		
		Rt[j][9]=-Rt[j][0]*XO[j][0]-Rt[j][1]*XO[j][1]-Rt[j][2]*XO[j][2];
		Rt[j][10]=-Rt[j][3]*XO[j][0]-Rt[j][4]*XO[j][1]-Rt[j][5]*XO[j][2];
		Rt[j][11]=-Rt[j][6]*XO[j][0]-Rt[j][7]*XO[j][1]-Rt[j][8]*XO[j][2];
	}
}



// generate a scene point
void scenePoint(double XP[3])
{
	const double d=1, fov=PI/4, w=2*d*tan(0.5*fov), h=w*0.8, depth=rnd(0,1);
	XP[0]=0.5*rnd(-w, w);
	XP[1]=0.5*rnd(-h, h);
	XP[2]=rnd(d, d+depth);
}



// map scene points onto image planes
void project(const double P[_V1][12], double data[_V][3][_P])
{
	double XQ[_P][3];
	for (int i=0; i<_P; ++i)
	{
		scenePoint(XQ[i]);
		const double fac=1./XQ[i][2];
		data[0][0][i]=XQ[i][0]*fac;
		data[0][1][i]=XQ[i][1]*fac;
		data[0][2][i]=1.;
		for (int j=0; j<_V1; ++j)
		{
			double PjXQ[3];
			PjXQ[0]=P[j][0]*XQ[i][0]+P[j][1]*XQ[i][1]+P[j][2]*XQ[i][2]+P[j][9];
			PjXQ[1]=P[j][3]*XQ[i][0]+P[j][4]*XQ[i][1]+P[j][5]*XQ[i][2]+P[j][10];
			PjXQ[2]=P[j][6]*XQ[i][0]+P[j][7]*XQ[i][1]+P[j][8]*XQ[i][2]+P[j][11];
			const double fac=1./PjXQ[2];
			data[j+1][0][i]=PjXQ[0]*fac;
			data[j+1][1][i]=PjXQ[1]*fac;
			data[j+1][2][i]=1.;
		}
	}
}



// generate synthetic data and ground truth camera matrices
void synthData(double data[_V][3][_P], Camera &cam)
{
	double XO[_V1][3];

	cameraCenters(XO);
	cameraMatrices(XO, cam.Rt);
	project(cam.Rt, data);

	// normalize translation vectors so that ||t2||=1
	const double scale=1./sqrt(cam.Rt[0][9]*cam.Rt[0][9]+cam.Rt[0][10]*cam.Rt[0][10]+cam.Rt[0][11]*cam.Rt[0][11]);
	for (int j=0; j<_V1; ++j)
	{
		cam.Rt[j][9]*=scale;
		cam.Rt[j][10]*=scale;
		cam.Rt[j][11]*=scale;
	}

	cam.Err=0;
}