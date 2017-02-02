#ifndef __ZETA_H__
#define __ZETA_H__

#include <cstdlib>
#include <iostream>

#include "myfloat.h"

using std::pow;
using std::cout;
using std::endl;

// Templated function type
#define FloatType template <class Float> const Float
// Max range to carry out sums
#define bigbound Float(10)
// Math
#define q Float(sqrt(q2))
#define _pi pi<Float>()
// Change this if needed- final answer will be answer to 1 part in FRAC_ACC .
// Note that this can be much less than the capacity of the myfloat type.
// myfloat is needed to evaluate summands in a numerically stable way, but after
// that the integral can be computed to whatever precision is desired. Assuming
// FRAC_ACC is less than myfloat's machine precision, it will bound the accuracy
// of luscher_zeta().
#define FRAC_ACC Float(pow(10,-4))

//SUMMANDS

FloatType inline small_kernel(const Float& n2, const Float& t) {
    return exp(-1*_pi*_pi*n2/t);
}

FloatType inline big_kernel(const Float& n2, const Float& t) {
    return exp(-1*t*n2);
}

FloatType inline inner_summand(const Float& n2, const Float& q2) {
	return Float(1)/(n2 - q2);
}

//SUMMATION

//Only the first octant in 3-space is summed over. All the summands are
//rotationally symmetric (depend only on magnitude of n), so this extends the
//sum to all  3-space by symmetry.
FloatType inline extend(int n1, int n2, int n3, const Float& val) {
	int zeros = 0;
	if (!n1) zeros++;
    if (!n2) zeros++;
    if (!n3) zeros++;
	
	if (zeros==0) return 8*val;
	if (zeros==1) return 4*val;
	if (zeros==2) return 2*val;
	return val;
}

//This is where the extended precision of myfloat is needed to evaluate a sum of
//many small contributions (values of small_kernel, big_kernel, or
//inner_summand).
FloatType sum(const Float (*func)(const Float&, const Float&), const Float& low, 
const Float& high, const Float& param) {
    Float result = 0;
    //Here is the only explicit dependence on a non-templated type in this file.
    //If you are not using myfloat.h, replace it with the floating-point type
    //of your choice. The reason is that Boost number types are objects, not
    //primitives,  so the usual (int)foo cast doesn't work with them.
    myfloat lambdafloor = trunc(high*high);
    int l = lambdafloor.convert_to<int>();
    myfloat low_ceil = ceil(low*low);
    int l2 = low_ceil.convert_to<int>();
    //This is the slowest part of the calculation, esp. evaluated many times
    // in romberg_integrate(), and using extend() to start these loops at 0 
    // instead of -l gives an 8x speedup.
    int m,i,j,k;
    for (i=0; i<=l; i++) {
    for (j=0; j<=l; j++) {
    for (k=0; k<=l; k++) {
        m = i*i + j*j + k*k;
    	if (m <= l && m >= l2) 
    	    result += extend(i,j,k, func((Float)m,param) );
    }}}
    return result;
}
        
//INTEGRANDS

//Let lambda (cutoff) = floor(q) to sum over all n with q^2 - n^2 > 0
//This includes as many terms as possible in the first sum.

//These represent the integral of e^(tq^2) K_lm^lambda(t,0) - 1/((4pi)^2 t^3/2)
//on t=[0,+INF] (eq. B.7, B.9). The t^-3/2 term is missing from one_to_inf()
//because it is  evaluated analytically and added back in luscher_zeta() to
//simplify these functions.  What would have been a non-integrable t^-3/2
//singularity at t=0 is subtracted off,  but there remains a t^-1/2 singularity
//at t=0 which is handled through a change of  variables in midpoint_*(). The
//integration range [1,INF] is made finite through  a different change of
//variables.

FloatType zero_to_one(const Float& q2, const Float& t) {
	return ( pow(4*_pi*t,-3./2.)*sum(&small_kernel,Float(0),bigbound,t) - 
		pow(2*_pi,-3)*sum(&big_kernel,Float(0),q,t) )*
		pow(4*_pi,-0.5)*exp(t*q2) - pow(4*_pi,-2)*pow(t,-3./2.);
}
FloatType one_to_inf(const Float& q2, const Float& t) {
	return pow(2*_pi,-3)*sum(&big_kernel,q,bigbound,t)*exp(t*q2)*pow(4*_pi,-0.5);
}


//INTEGRATION

// The midpoint and Romberg algorithms are from Numerical Recipes in C, 2nd. ed.
// The changes of variables are from Numerical Recipes, 3rd. ed.

//Identity function (no change of variable)
FloatType inline identity(const Float (*func)(const Float&, const Float&), 
const Float& q2, const Float& t) {
    return func(q2,t);
}

//Change of variable - for integrals on [x,+INF] that decay exponentially
FloatType inline exponential(const Float (*func)(const Float&, const Float&), 
const Float& q2, const Float& t) {
    return ( (*func)(q2,-log(t)) )/t;
}


//Change of variable - for integrals on [0,x] with a 1/sqrt singularity at 0
FloatType inline invsqrt(const Float (*func)(const Float&, const Float&), 
const Float& q2, const Float& t) {
    return 2*t*( (*func)(q2,t*t) );
}

//A new COV can be implemented easily by defining it as above, adding its name
//to this enum type, and adding its case in romberg_integrate()
enum midpt {normal,cov_exponential,cov_invsqrt};

//The finite difference used is an n-point midpoint method because it does not
//evaluate f at the endpoints, so even if the function blows up at the
//endpoints a decent answer can be gotten without knowing an appropriate COV
//to suppress it. Remembering the previous function call in a static variable
//avoids unnecessarily duplicating work, see discussion in Numerical Recipes
FloatType midpoint(const Float (*func)(const Float&, const Float&), 
const Float (*transform)(const Float (*)(const Float&, const Float&), const Float&, const Float&),
const Float& q2, const Float& a, const Float& b, const int n) {
    static Float s; //result of midpoint(n-1)
    if (n==1) {
        Float pt = 0.5*(a+b);
        return (s=(b-a)*transform(func,q2,pt));
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
            sum += transform(func,q2,t);
            t += spacing2;
            sum += transform(func,q2,t);
            t += spacing1;
        }
        s=(s+(b-a)*sum/num_points)/Float(3);
        return s;
    }
}

//The driver routine - evaluate a finite difference several times for increasing
//n  (decreasing stepsize h) and extrapolate to the h=0 limit.
FloatType romberg_integrate(int NUM_PTS, const Float (*func)(const Float&, const Float&),
const Float& q2, const Float& aa, const Float& bb, midpt rule) {
	int MAX_STEPS = 14;
    //stores calls to midpoint() and their relative stepsizes, 1-indexed
    Float s[MAX_STEPS+2], h[MAX_STEPS+2];
	h[1] = 1.0;
    Float a,b;
    switch (rule) {
        case cov_exponential : a = 0; b = exp(-aa); break;
        case cov_invsqrt     : a = 0; b = sqrt(bb); break;
        case normal          : a = aa; b = bb;      break;
    }

	for (int j=1; j<=MAX_STEPS; j++) {

		switch (rule) {
			case cov_exponential : s[j] = midpoint(func,exponential,q2,a,b,j); break;
			case cov_invsqrt     : s[j] = midpoint(func,invsqrt,q2,a,b,j);     break;
			case normal          : s[j] = midpoint(func,identity,q2,a,b,j);    break;
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
FloatType luscher_zeta_integralterm(const Float& q2) {
    Float int1 = romberg_integrate(5,zero_to_one,q2,Float(0),Float(1),cov_invsqrt);
    Float int2 = romberg_integrate(5,one_to_inf,q2,Float(1),Float(50),cov_exponential);
    Float constterm = -pow(_pi,-2)/Float(8);
    return pow(2*_pi,3)*( constterm + int1 + int2 );
 }
 
//printing (optional)
template <class Float> void lzprint(const Float& q2, const Float& sum1,
const Float& intl, const Float& ans) {  
    if (abs(ans)>5*pow(10,1)) cout<<'#';
    cout<<q2<<' '<<sum1<<' '<<intl<<' '<<ans<<endl;
}

//driver routine
FloatType _luscher_zeta00(const Float& q2, bool reduced) {
    if (q2<0)
        std::cerr<<"luscher_zeta: nonnegative arguments only."<<endl;

    Float sumbound;
    if (reduced) sumbound = (q2<1) ? Float(0) : sqrt(q2-1); //max(sqrt(q2-1)),0)
    else sumbound = q;
    Float sum1 = sum(inner_summand,Float(0),sumbound,q2) / sqrt(4*_pi);
    Float intl = luscher_zeta_integralterm(q2);
    Float ans = sum1 + intl;

    lzprint(q2,sum1,intl,ans);
    
    return ans;
}

//The final answer. Should be self-explanatory. This function looks much like
//tan((q2+0.5)pi), with poles at each integer q2. It also  gets steeper with
//increasing q2 and, interestingly, sometimes skips its integer poles. The poles
//are caused by the sum1 term for |n|=q, and for some q2 there is no n that
//satisfies this (that is, there is no n=(a,b,c) s.t. q2=a^2+b^2+c^2). The q2
//that  have a pole are given by OEIS:A000378. In Luscher's original paper
//(1990) this is called the "singular set" (eq. 3.3)
FloatType luscher_zeta00(const Float& q2) {
    return _luscher_zeta00(q2,false);
}

FloatType reduced_luscher_zeta00(const Float& q2) {
    return _luscher_zeta00(q2,true);
}

#endif
