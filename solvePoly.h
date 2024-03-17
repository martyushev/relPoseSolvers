#define THIRD 0.3333333333333333333
#define SQRT3 1.7320508075688772935
#define MAXIT1 50 // number of iterations for root isolating and polishing
#define ACC 1.e-14 // root is polished to an approximate accuracy ACC
#define DEG 6 // degree of polynomial



// all real roots of a polynomial p(x) lie in (bound[0], bound[1])
// Kioustelidis' method
void bounds(const double p[DEG], double bound[2])
{
	double M=fabs(p[0]);
	bound[0]=bound[1]=0.;
	for (int i=0; i<DEG; ++i)
	{
		const double c=fabs(p[i]);
		if (c>M) M=c;
		if (p[i]<0)
		{
			const double t=pow(-p[i], 1./(double)(DEG-i));
			if (t>bound[1]) bound[1]=t;
			if (t>bound[0] && !(i%2)) bound[0]=t;
		}
		else if (i%2)
		{
			const double t=pow(p[i], 1./(double)(DEG-i));
			if (t>bound[0]) bound[0]=t;
		}
	}
	M+=1.; // Cauchy bound
	bound[0]*=-2.;
	if (bound[0]<-M) bound[0]=-M;
	bound[1]*=2.;
	if (bound[1]>M) bound[1]=M;
}



// get data to compute the Sturm sequence for a polynomial p(x) at any point
void quotients(const double p[DEG], const double dp[DEG-1], double q[DEG+1][2])
{
	double r_[DEG], r[DEG-1];
	
	r_[DEG-1]=1.;
	q[DEG][1]=p[DEG-1]-dp[DEG-2];
	r_[0]=dp[0];
	r[0]=dp[0]*q[DEG][1]-p[0];
	for (int j=1; j<DEG-1; ++j)
	{
		r_[j]=dp[j];
		r[j]=dp[j]*q[DEG][1]+r_[j-1]-p[j];
	}

	for (int i=DEG-1; i>=2; --i)
	{
		const int i1=i-1;
		const double ri1=1./r[i1];
		q[i][0]=r_[i]*ri1;
		q[i][1]=(r_[i1]-q[i][0]*r[i-2])*ri1;
		const double t=r_[0];
		r_[0]=r[0];
		r[0]=r[0]*q[i][1]-t;
		for (int j=1; j<i1; ++j)
		{
			const double t=r_[j];
			r_[j]=r[j];
			r[j]=r[j]*q[i][1]+q[i][0]*r_[j-1]-t;
		}
		r_[i1]=r[i1];
	}

	q[1][0]=r_[1];
	q[1][1]=r_[0];
	q[0][1]=r[0];
}



// evaluate polynomial p(x) at x0
double evalPoly(const double p[DEG], const double &x0)
{
	double s=x0+p[DEG-1];
	for (int i=DEG-2; i>=0; --i) s=s*x0+p[i];
	return s;
}



// number of sign changes in a sequence seq[]
int nchanges(const double seq[DEG+1])
{
	int s1, s2=(seq[0]>0.)? 1:((seq[0]<0.)? -1:0), s=0;
	for (int i=1; i<DEG+1; ++i)
	{
		if (!seq[i]) continue;
		s1=s2;
		s2=(seq[i]>0.)? 1:-1;
		if (!(s1+s2)) ++s;
	}
	return s;
}



// evaluate Sturm sequence at a
int evalSturmSeq(const double q[DEG+1][2], const double &a)
{
	double sa[DEG+1];
	// initialize sa
	sa[0]=q[0][1];
	sa[1]=q[1][0]*a+q[1][1];
	// compute sa recursively
	for (int i=2; i<DEG; ++i)
		sa[i]=(q[i][0]*a+q[i][1])*sa[i-1]-sa[i-2];
    sa[DEG]=(a+q[DEG][1])*sa[DEG-1]-sa[DEG-2]; // since q[DEG][0]=1.
	return nchanges(sa);
}



// isolate all real roots of polynomial p(x)
int isolateRoots(const double p[DEG], const double dp[DEG-1], double Isol[DEG][2])
{
	int nIsol=0, nTree=1, nIters=1, min=0, sTree[2*MAXIT1+1][2];
	double Tree[2*MAXIT1+1][2], q[DEG+1][2];

	// initialize the tree
	bounds(p,Tree[0]); // all real roots of polynomial p(x) lie in (Tree[0][0], Tree[0][1])

	quotients(p,dp,q);
	sTree[0][0]=evalSturmSeq(q,Tree[0][0]);
	sTree[0][1]=evalSturmSeq(q,Tree[0][1]);

	while (nTree>min)
	{
		const double a=Tree[min][0], b=Tree[min][1];
		const int sa=sTree[min][0], sb=sTree[min][1];
		const int s=sa-sb; // counts the number of real roots in (a, b)
		++min;

		if (s==1)
		{ // an isolated root found
			Isol[nIsol][0]=a;
			Isol[nIsol++][1]=b;
		}
		else if (s>1)
		{ // proceed to make subdivision
			const int nTree1=nTree+1;
			const double mid=0.5*(a+b);
			// add intervals (a, mid] and (mid, b] to Isol
			Tree[nTree][1]=Tree[nTree1][0]=mid;
			sTree[nTree][1]=sTree[nTree1][0]=evalSturmSeq(q,mid);
			Tree[nTree][0]=a;
			Tree[nTree1][1]=b;
			sTree[nTree][0]=sa;
			sTree[nTree1][1]=sb;
			nTree+=2;
			++nIters;
			if (nIters>MAXIT1)
			{ // perhaps some roots are too close
				Isol[nIsol][0]=a;
				Isol[nIsol++][1]=mid;
				Isol[nIsol][0]=mid;
				Isol[nIsol++][1]=b;
				return nIsol;
			}
		}
	}
	return nIsol;
}



// polish the root of a polynomial p(x) known to lie between xl and x2, Ridders' method
// this function adapted from "Numerical recipes in C" by Press et al.
// the output is either 1 or 0 (no solution found)
bool polishRoots(const double p[DEG], const double &x1, const double &x2, double &ans)
{
	double fl=evalPoly(p,x1), fh=evalPoly(p,x2);
	if (!fh)
    {
        ans=x2;
        return 1;
    }

	if ((fl>0)? (fh<0):(fh>0))
	{
		double xl=x1, xh=x2;
		ans=0.5*(x1+x2);
		for (int j=1; j<=MAXIT1; ++j)
		{
			const double xm=0.5*(xl+xh), fm=evalPoly(p,xm), s=sqrt(fm*fm-fl*fh);
			if (!s) return 1;
			const double xnew=(fl<fh)? xm+(xl-xm)*fm/s:xm+(xm-xl)*fm/s;
			if (fabs(xnew-ans)<=ACC) return 1;
			ans=xnew;
			const double fnew=evalPoly(p,ans);
			if (!fnew) return 1;
			if (fnew>=0? (fm<0):(fm>0))
			{
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			}
			else if (fnew>=0? (fl<0):(fl>0))
			{
				xh=ans;
				fh=fnew;
			}
			else
			{
				xl=ans;
				fl=fnew;
			}
			if (fabs(xh-xl)<=ACC) return 1;
		}

		return 0;
	}
	else return 0;
}



// solve for real roots of a square-free polynomial p(x) of degree DEG
// output is the number of roots found
// for polynomials of degree 2,3,4, use closed-form solvers for efficiency
int solvePoly(const double p[DEG+1], double x[DEG])
{
	// copy and normalize the input polynomial p(x) and its derivative dp(x)
	double p1[DEG], dp1[DEG-1];
	const double pdeg=1./p[DEG], dpdeg=1./(double)DEG;
	p1[0]=p[0]*pdeg;
	for (int i=1; i<DEG; ++i)
	{
		p1[i]=p[i]*pdeg;
		dp1[i-1]=(double)i*p1[i]*dpdeg;
	}	

	double Isol[DEG][2];
	const int nIsol=isolateRoots(p1,dp1,Isol); // isolate all real roots of p(x)

	int nr=0;
	for (int i=0; i<nIsol; ++i)
	{ // find isolated real root of p(x)
		if (!polishRoots(p1,Isol[i][0],Isol[i][1],x[nr])) continue;
		++nr;
	}

	return nr;
}



// solve for real roots of quadratic polynomial p[2]*x^2 + p[1]*x + p[0]
// output is the number of roots found
int solveQuadratic(const double p[3], double x[2])
{
	const double discr=p[1]*p[1]-4.*p[0]*p[2];
	if (discr>0)
	{ // 2 simple real roots
		const double q=-0.5*(p[1]+((p[1]>0)? p[1]:-p[1])*sqrt(discr));
		x[0]=q/p[2];
		x[1]=p[0]/q;
		return 2;
	}
	else if (discr<0) return 0;
	else // discr==0
	{ // 2 multiple real roots
		x[0]=-0.5*p[1]/p[2];
		return 1;
	}

}


// solve for real roots of cubic polynomial x^3 + a[2]*x^2 + a[1]*x + a[0]
// output is the number of roots found
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
		{ // 1 simple and 2 multiple real roots
			const double s=(g>0)? pow(g,THIRD):-pow(-g,THIRD);
			x[0]=2.*s+p; // simple root
			x[1]=p-s; // multiple root
			return 2;
		}
		else
		{ // 3 multiple real roots
			x[0]=p;
			return 1;
		}
	}
}



// solve for real roots of quartic polynomial a[4]*x^4 + a[3]*x^3 + a[2]*x^2 + a[1]*x + a[0]
// output is the number of roots found
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

	double y_pos=0; // find a positive root of the cubic equation
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
	else if (f1<0) nr4=0;
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
	else if (f2<0) return nr4;
	else // f2==0
		x[nr4++]=k-b[3]; // 2 more multiple real roots
	
	return nr4;
}



#undef THIRD
#undef SQRT3
#undef MAXIT1
#undef ACC
#undef DEG