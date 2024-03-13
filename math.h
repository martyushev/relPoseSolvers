#define PI 3.1415926535897932385
#define RATIO 0.6180339887498948482
#define THIRD 0.3333333333333333333
#define SQRT3 1.7320508075688772935



// uniformly distributed random numbers from [a, b]
inline double rnd(const double &a, const double &b)
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



// multiply two 3x3 matrices C=A*B
void mult(const double A[9], const double B[9], double C[9])
{
	for (int i=0; i<3; ++i)
	{
		const int j=3*i;
		for (int k=0; k<3; ++k)
			C[j+k]=A[j]*B[k]+A[j+1]*B[k+3]+A[j+2]*B[k+6];
	}
}



// basis of null space of mxn matrix A using specifically tailored QR factorization
// result is n-m n-vectors Q
// note that matrix A changes during the computation
void nullQR(double A[9][2*NPOINTS], double Q[9][2], const int n, const int m)
{
	const int nm=n-m, m1=m-1;

	// Hauseholder vectors
	for (int j=0; j<m; ++j)
	{
		double t=0;
		for (int i=j+1; i<n; ++i) t+=A[i][j]*A[i][j];
		double mu=sqrt(A[j][j]*A[j][j]+t);
		mu=(A[j][j]<0)? 1./(A[j][j]-mu):1./(A[j][j]+mu);

		for (int i=j+1; i<n; ++i) A[i][j]*=mu;

		double beta=-2./(1.+t*(mu*mu));
		for (int k=j+1; k<m; ++k)
		{
			double w=A[j][k];
			for (int i=j+1; i<n; ++i) w+=A[i][k]*A[i][j];
			w*=beta;
			for (int i=j+1; i<n; ++i) A[i][k]+=A[i][j]*w;
		}
	}

	// multiply m Householder matrices, inverse order is more efficient
	// we only need last n-m columns of the resulting matrix
	for (int k=0; k<m1; ++k) memset(Q[k],0,sizeof(double)*m1);
	
	// start from the mth matrix
	double beta=1.;
	for (int i=m; i<n; ++i) beta+=A[i][m1]*A[i][m1];

	beta=-2./beta;
	for (int k=0; k<nm; ++k)
	{
		Q[m1][k]=A[k+m][m1]*beta;
		for (int i=m; i<n; ++i)
			Q[i][k]=(k==i-m)? 1.+A[i][m1]*Q[m1][k]:A[i][m1]*Q[m1][k];
	}
	
	// multiply by the remaining m-1 matrices
	for (int j=m1-1; j>=0; --j)
	{
		double beta=1.;
		for (int i=j+1; i<n; ++i) beta+=A[i][j]*A[i][j];
		beta=-2./beta;
		for (int k=0; k<nm; ++k)
		{
			double w=Q[j][k];
			for (int i=j+1; i<n; ++i) w+=Q[i][k]*A[i][j];
			w*=beta;
			Q[j][k]+=w;
			for (int i=j+1; i<n; ++i) Q[i][k]+=A[i][j]*w;
		}
	}
}



// compute adjoint to 3x3 matrix A, i.e. det(A) A^{-1}
void adjoint(const double A[9], double B[9], bool sym)
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



// transpose of 3x3 matrix A
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



// solve for real roots of cubic polynomial x^3 + a[2]*x^2 + a[1]*x + a[0]
int solveCubic(const double a[3], double x[3])
{
	const double p=-THIRD*a[2], p2=p*p, f=p2-THIRD*a[1], f3=f*f*f;
	const double g=(p2-0.5*a[1])*p-0.5*a[0], discr=g*g-f3; // (minus) discriminant
	if (discr<0)
	{ // 3 simple real roots
		const double t=sqrt(f), k=g/(f*t), phi=THIRD*acos(k); // phi in [0, pi/3]
		const double s=cos(phi)*t, t1=p-s, t2=SQRT3*sqrt(f-s*s);
		x[0]=2.*s+p;
		x[1]=t1-t2;
		x[2]=t1+t2;
		return 3;
	}
	else if (discr>0)
	{ // 1 real and 2 complex conjugate roots
		const double t=sqrt(discr), s1=((g>0)-(t>0))? g-t:g+t, s2=f3/s1;
		const double t1=(s1>0)? pow(s1,THIRD):-pow(-s1,THIRD), t2=f/t1;
		x[0]=t1+t2+p;
		//x[1]=-0.5*(t1+t2)+p;  // real part of the complex conjugate pair, imaginary part is 0.5*(t1-t2)*SQRT3
		return 1;
	}
	else // discr==0
	{
		if (g)
		{ // 2 multiple real roots
			const double s=(g>0)? pow(g,THIRD):-pow(-g,THIRD);
			x[0]=2.*s+p;
			x[1]=x[2]=p-s;
			return 2;
		}
		else
		{ // 3 multiple real roots
			x[0]=x[1]=x[2]=p;
			return 1;
		}
	}
}



// solve for real roots of quartic polynomial a[4]*x^4 + a[3]*x^3 + a[2]*x^2 + a[1]*x + a[0]
int solveQuartic(const double a[5], double x[4])
{
	const double a4=0.25/a[4], b[4]={a[0]*a4, a[1]*a4, a[2]*a4, a[3]*a4};
	const double e=b[3]*b[3], q=b[2]-1.5*e, r=0.5*b[1]+(e-b[2])*b[3], s=b[0]-b[1]*b[3]+(b[2]-0.75*e)*e, q2=2.*q;

	if (!q && !r && !s)
	{ // 4 multiple real roots
		x[0]=-b[3];
		return 1;
	}
	
	const double c[3]={-r*r, q*q-s, q2};
	double y[3];

	const int nr3=solveCubic(c,y);

	double y_pos=0; // find positive root of the cubic equation
	for (int i=0; i<nr3; ++i)
		if (y[i]>0)
		{
			y_pos=y[i];
			break;
		}
	if (!y_pos) return 0;

	const double k=sqrt(y_pos), l=y_pos+q2, m=2.*r/k, f1=m-l, f2=-m-l;
	int nr4;
	if (f1>0)
	{ // 2 simple real roots
		const double d=sqrt(f1), t1=-k-b[3];
		x[0]=t1-d;
		x[1]=t1+d;
		nr4=2;
	}
	else if (f1<0)
		nr4=0;
	else // f1==0
	{ // 2 multiple real roots
		x[0]=-k-b[3];
		nr4=1;
	}

	if (f2>0)
	{ // 2 more simple real roots
		const double d=sqrt(f2), t2=k-b[3];
		x[nr4++]=t2-d;
		x[nr4++]=t2+d;
	}
	else if (f2<0)
		return nr4;
	else // f2==0
		x[nr4++]=k-b[3];
	
	return nr4;
}