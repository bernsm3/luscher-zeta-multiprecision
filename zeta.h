#ifndef __ZETA_H__
#define __ZETA_H__

#include <cstdlib>
#include <iostream>

#include "myfloat.h"

using std::pow;
using std::cout;
using std::endl;

#define bigbound Float(10)
#define q Float(sqrt(q2))
#define FloatType template <class Float> const Float
#define _pi pi<Float>()

//SPHERICAL HARMONIC POLYNOMIALS

template <class T> 
class Point {
public:
    Point(const T& x, const T& y, const T& z) {
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z/r);
        phi = atan(y/x);
    }

    Point(const Point& p, const T& a): r(a * p.r), theta(p.theta), phi(p.phi) {}
    
    const Point& times(const T& a) const { 
        Point* p = new Point(*this,a); 
        return *p; 
    }
    
    const T r2() const { return r*r; }
    
    T r;
    T theta;
    T phi;
};


FloatType inline y(int l, int m, const Point<Float>& n) {
    return pow(n.r,l)*spherical_harmonic_r(l,m,n.theta,n.phi);
}

//SUM ARGUMENTS

FloatType inline small_kernel(int l, int m, const Point<Float>& n, const Float& t) {
    return y(l, m, n.times(2*_pi))*exp(-1*_pi*_pi*n.r2() /t);
}

FloatType inline big_kernel(int l, int m, const Point<Float>& n, const Float& t) {
    return y(l,m,n)*exp(-1*t*n.r2());
}

FloatType inline inner_summand(int l, int m, const Point<Float>& n, const Float& q2) {
	return y(l,m,n)/(n.r2() - q2);
}

//SUMMATION

//This is where the extended precision of myfloat is needed to evaluate a sum of
// many small contributions (values of small_kernel, big_kernel, or inner_summand). 
//More discussion by romberg_integrate().
FloatType sum(const Float (*func)(int, int, const Point<Float>&, const Float&),
const Float& low, const Float& high, int l, int m, const Float& param) {
    Float result = 0;
    //Here is the only explicit dependence on a non-templated type in this file.
    //If you are not using myfloat.h, replace it with the floating-point type
    // of your choice. The reason is that Boost number types are objects, not primitives,
    // so the usual (int)foo cast doesn't work with them.
    myfloat lambdafloor = trunc(high*high);
    int ubound = lambdafloor.convert_to<int>();
    myfloat low_ceil = ceil(low*low);
    int lbound = low_ceil.convert_to<int>();
    //This is the slowest part of the calculation, esp. evaluated many times
    // in romberg_integrate()
    int magnitude,i,j,k;
    for (i=-ubound; i<=ubound; i++) {
    for (j=-ubound; j<=ubound; j++) {
    for (k=-ubound; k<=ubound; k++) {
        magnitude = i*i + j*j + k*k;
    	if (magnitude <= ubound && magnitude >= lbound) {
            Point<Float> p = Point<Float>( Float(i), Float(j), Float(k) );
    	    result += func(l,m,p,param);
        }
    }}}
    return result;
}
        
//INTEGRANDS

//Let lambda (cutoff) = floor(q) to sum over all n with q^2 - n^2 > 0
//This includes as many terms as possible in the 

//These represent the integral of e^(tq^2) K_lm^lambda(t,0) - 1/((4pi)^2 t^3/2)
// on t=[0,+INF] (eq. B.7, B.9). The t^-3/2 term is missing from one_to_inf() because it is
// evaluated analytically and added back in luscher_zeta() to simplify these functions.
// What would have been a non-integrable t^-3/2 singularity at t=0 is subtracted off,
// but there remains a t^-1/2 singularity at t=0 which is handled through a change of
// variables in midpoint_*(). The integration range [1,INF] is made finite through
// a different change of variables.

//R&G eq. 110
FloatType zero_to_one(int l, int m, const Float& q2, const Float& t) {
    Float constterm;
    if (!l && !m) 
        constterm = pow(4*_pi,-2)*pow(t,-3./2.);
    else constterm = 0;
    int i;
    if (l%4 == 0) i=1;
    else if (l%4 == 1 || l%4 == 3) i=0;
    else i=-1;
	return ( (Float(i)/pow(2*t,l))*pow(4*_pi*t,-3./2.)*sum(&small_kernel,Float(0),bigbound,l,m,t) - 
		pow(2*_pi,-3)*sum(&big_kernel,Float(0),q,l,m,t) )*
		exp(t*q2) - constterm;
}

//R&G eq. 109
FloatType one_to_inf(int l, int m, const Float& q2, const Float& t) {
	return pow(2*_pi,-3)*sum(&big_kernel,q,bigbound,l,m,t)*exp(t*q2);
}


//INTEGRATION

// The midpoint and Romberg algorithms are from Numerical Recipes in C, 2nd. ed.
// The changes of variables are from Numerical Recipes, 3rd. ed.

//Identity function (no change of variable)
FloatType inline identity(const Float (*func)(int, int, const Float&, const Float&), 
int l, int m, const Float& q2, const Float& t) {
    return func(l,m,q2,t);
}

//Change of variable - for integrals on [x,+INF] that decay exponentially
FloatType inline exponential(const Float (*func)(int, int, const Float&, const Float&), 
int l, int m, const Float& q2, const Float& t) {
    return ( func(l,m,q2,-log(t)) )/t;
}


//Change of variable - for integrals on [0,x] with a 1/sqrt singularity at 0
FloatType inline invsqrt(const Float (*func)(int, int, const Float&, const Float&), 
int l, int m, const Float& q2, const Float& t) {
    return 2*t*( (*func)(l,m,q2,t*t) );
}

//A new COV can be implemented easily by defining it as above, adding its name to this
// enum type, and adding its case in romberg_integrate()
enum midpt {normal,cov_exponential,cov_invsqrt};

//The finite difference used is an n-point midpoint method because it does not evaluate
// f at the endpoints, so even if the function blows up at the endpoints a decent 
// answer can be gotten without knowing an appropriate COV to suppress it.
//Remembering the previous function call in a static variable avoids unnecessarily 
// duplicating work, see discussion in Numerical Recipes
FloatType midpoint(const Float (*func)(int, int, const Float&, const Float&), 
const Float (*transform)(const Float (*)(int, int, const Float&, const Float&), int, int, const Float&, const Float&),
int l, int m, const Float& q2, const Float& a, const Float& b, const int n) {
    static Float s; //result of midpoint(n-1)
    if (n==1) {
        Float pt = 0.5*(a+b);
        return (s=(b-a)*transform(func,l,m,q2,pt));
    } else {
        int it, i;
        Float t, num_points, sum, spacing1, spacing2;
        for (it=1,i=1; i<n-1; i++)
            it *= 3;
        num_points = it;
        spacing1 = (b-a)/(3*num_points);
        spacing2 = spacing1 + spacing1;
        t = a + 0.5*spacing1;
        sum = 0;
        for (i=1; i<=it; i++) {
            sum += transform(func,l,m,q2,t);
            t += spacing2;
            sum += transform(func,l,m,q2,t);
            t += spacing1;
        }
        s=(s+(b-a)*sum/num_points)/Float(3);
        return s;
    }
}

//The driver routine - evaluate a finite difference several times for increasing n
// (decreasing stepsize h) and extrapolate to the h=0 limit.
//Note that the precision of the answer, FRAC_ACC ("fractional accuracy") can be much 
// less than the capacity of the myfloat type. myfloat is needed to evaluate func
// (zero_to_one or one_to_inf) in a numerically stable way, but after that the integral
// can be computed to whatever precision is desired. This means that assuming FRAC_ACC
// is less than myfloat's machine precision, it will bound the final accuracy of
// luscher_zeta(). If myfloat gives 20 decimal digits of precision but FRAC_ACC is 10^-4,
// luscher_zeta() will be accurate to 4 digits, not 20!
FloatType romberg_integrate(int NUM_PTS, const Float (*func)(int, int, const Float&, const Float&),
int l, int m, const Float& q2, const Float& aa, const Float& bb, midpt rule) {
	//int NUM_PTS = 6;
	int MAX_STEPS = 14;
	Float FRAC_ACC = pow(10,-4);
    Float s[MAX_STEPS+2], h[MAX_STEPS+2]; //stores calls to midpoint() and their relative stepsizes, 1-indexed
	h[1] = 1.0;
    Float a,b;
    switch (rule) {
        case cov_exponential : a = 0; b = exp(-aa); break;
        case cov_invsqrt     : a = 0; b = sqrt(bb); break;
        case normal          : a = aa; b = bb;      break;
    }

	for (int j=1; j<=MAX_STEPS; j++) {

		switch (rule) {
			case cov_exponential : s[j] = midpoint(func,exponential,l,m,q2,a,b,j); break;
			case cov_invsqrt     : s[j] = midpoint(func,invsqrt,l,m,q2,a,b,j);     break;
			case normal          : s[j] = midpoint(func,identity,l,m,q2,a,b,j);    break;
		}
		if (j>=NUM_PTS) {
			//Neville's algorithm for polynomial extrapolation
			int table_idx = 1;
			int zero = j-NUM_PTS; //start of the used part of the arrays
			Float diff = abs(h[zero+1]);
			Float c[NUM_PTS+1];
			Float d[NUM_PTS+1];
			Float dift;
			for (int i=1; i<=NUM_PTS; i++) {
				if ( (dift=abs(h[zero+i])) < diff ) {
            		table_idx = i;
            		diff = dift;
            	}
				c[i] = s[zero+i];
				d[i] = s[zero+i];
			}
			Float ans = s[zero + table_idx--];
			Float ho, hp, w, error;
			for (int m=1; m<NUM_PTS; m++) {
				for (int i=1; i<=NUM_PTS-m; i++) {
					ho = h[zero+i];
					hp = h[zero+i+m];
					w = c[i+1] - d[i];
					if ((ho-hp) <= epsilon<Float>()) //floating point equality
                        std::cerr << "Roundoff error in romberg_integrate" << endl;
					Float den = w/(ho-hp);
					d[i] = hp*den;
					c[i] = ho*den;
				}
				ans += (error=(2*table_idx < (NUM_PTS-m) ? c[table_idx+1] : d[table_idx--]));
			}
			if (abs(error) < FRAC_ACC*abs(ans)) {
				return ans;
			}
		}
		s[j+1] = s[j];
		h[j+1] = h[j]/Float(9);
	}
    std::cerr << "romberg_integrate failed to converge- too many steps" << endl;
    return 0;
}

//ZETA

// All sum terms with n^2>q^2, converted to the integral over K.
FloatType luscher_zeta_integralterm(int l, int m, const Float& q2) {
	Float int1 = romberg_integrate(5,zero_to_one,l,m,q2,Float(0),Float(1),cov_invsqrt);
	Float int2 = romberg_integrate(5,one_to_inf,l,m,q2,Float(1),Float(50),cov_exponential);
    Float constterm;
    if (!l && !m) 
        constterm = -pow(_pi,-2)/Float(8);
    else constterm = 0;
    return pow(2*_pi,3)*( constterm + int1 + int2 );
 }
 
//printing (optional)
template <class Float> void lzprint(const Float& q2, const Float& sum1,
const Float& intl, const Float& ans) {  
    if (abs(ans)>5*pow(10,1)) cout<<'#'; //y-value cutoff for plotting
    cout<<q2<<' '<<sum1<<' '<<intl<<' '<<ans<<endl;
}

//driver routine
FloatType _luscher_zeta(int l, int m, const Float& q2, bool reduced) {
    if (q2<0)
        std::cerr<<"luscher_zeta: nonnegative arguments only."<<endl;

    Float sumbound;
    if (reduced) sumbound = (q2<1) ? Float(0) : sqrt(q2-1); //max(sqrt(q2-1)),0)
    else sumbound = q;
    Float sum1 = sum(inner_summand,Float(0),sumbound,l,m,q2);
    Float intl = luscher_zeta_integralterm(l,m,q2);
    Float ans = sum1 + intl;

    lzprint(q2,sum1,intl,ans);
    
    return ans;
}

//The final answer. Should be self-explanatory.
//This function looks much like tan((q2+0.5)pi), with poles at each integer q2. It also
// gets steeper with increasing q2 and, interestingly, sometimes skips its integer poles.
//The poles are caused by the sum1 term for |n|=q, and for some q2 there is no n that
// satisfies this (that is, there is no n=(a,b,c) s.t. q2=a^2+b^2+c^2). The q2 that
// do not have a pole are given by the complement of OEIS:A000378.
FloatType luscher_zeta(int l, int m, const Float& q2) {
    return _luscher_zeta(l,m,q2,false);
}

FloatType reduced_luscher_zeta(int l, int m, const Float& q2) {
    return _luscher_zeta(l,m,q2,true);
}

#endif
