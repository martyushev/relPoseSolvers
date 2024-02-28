// solve for real roots of cubic polynomial x^3 + a[2]*x^2 + a[1]*x + a[0]
int solveCubic(const double a[3], double x[3])
{
	const double third=0.3333333333333333, sqrt3=1.732050807568877;
	const double p=-third*a[2], p2=p*p, f=p2-third*a[1], f3=f*f*f;
	const double g=(p2-0.5*a[1])*p-0.5*a[0], discr=g*g-f3; // (minus) discriminant
	if (discr<0)
	{ // 3 simple real roots
		const double t=sqrt(f), k=g/(f*t), phi=third*acos(k); // phi in [0, pi/3]
		const double s=cos(phi)*t, t1=p-s, t2=sqrt3*sqrt(f-s*s);
		x[0]=2.*s+p;
		x[1]=t1-t2;
		x[2]=t1+t2;
		return 3;
	}
	else if (discr>0)
	{ // 1 real and 2 complex conjugate roots
		const double t=sqrt(discr), s1=((g>0)-(t>0))? g-t:g+t, s2=f3/s1;
		const double t1=(s1>0)? pow(s1,third):-pow(-s1,third), t2=f/t1;
		x[0]=t1+t2+p;
		//x[1]=-0.5*(t1+t2)+p;  // real part of complex conjugate pair, imaginary part is 0.5*(t1-t2)*sqrt3
		return 1;
	}
	else // discr==0
	{
		if (g)
		{ // 2 multiple real roots
			const double s=(g>0)? pow(g,third):-pow(-g,third);
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

	const int nr3=solveCubic(c, y);

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