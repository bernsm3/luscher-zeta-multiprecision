#include <cmath>
#include <ctime>
#include <vector>

#include "myfloat.h"
// --- pick one of these ---
// #include "zeta00.h"
#include "zeta.h"
// -------------------------

using namespace std;

int main(int argc, char* argv[]) {

	cout << std::setprecision(std::numeric_limits<myfloat>::digits10);

    if (argc>4) {
            std::cerr << "Spurious arguments. Usage: ./zeta -log(stepsize) OR "
		  << "./zeta start end -log(stepsize) OR ./zeta start end" << endl;
    }
    int SCALE = 3;
    if (argc==2) SCALE = atoi(argv[1]);
    if (argc==4) SCALE = atoi(argv[3]);
    myfloat step = pow(10,-SCALE);
    myfloat i, end;
    i = 0.005; end=1;
    if (argc>2) { i=atof(argv[1]); end=atof(argv[2]); }
    for (; i<end; i+=step) {
    	luscher_zeta<myfloat>(0,0,i);
    	// luscher_zeta00<myfloat>(i);
    }

	return 0;
}


// ------------
//An implementation of the example in sec. 1 of "Two-particle states on a torus".

// condition (a), eqn. 1.3
complex<myfloat> zeta_ratio(myfloat q2) {
	myfloat re = luscher_zeta<myfloat>(0,0,q2);
	// myfloat re = luscher_zeta00<myfloat>(q2);
	myfloat im = pow(pi<myfloat>(),3.0/2.0)*sqrt(q2);
	return complex<myfloat>(re,im)/complex<myfloat>(re,-1*im);
}

// condition (b)
//Return true if a squared magnitude is in A^{++} i.e. if it contains 2 
//integer-valued 3-vectors independent under cubic rotation, reflection. 
//In other words, checks r_squared for membership in OEIS:A124966
bool check_sqr_mag(int r_squared) {
	vector<int> v1(3,0), v2(3,0);
	int r = (int)trunc(sqrt(r_squared));
	for (int i=0; i<r; i++) {
		for (int j=i; j<r; j++) {
			int running_sum = i*i + j*j;
			//if (running_sum < r_squared) continue;
			float k = sqrt(r_squared - running_sum);
			if (k!=trunc(k)) continue;
			//v1[0] = i; v1[1] = j; v1[2] = (int)k;
			for (int a=0; a<r; a++) {
				if (a==i) continue;
				for (int b=a; b<r; b++) {
					if (b==j) continue;
					int sum2 = a*a + b*b;
					float c = sqrt(r_squared-sum2);
					if (c==trunc(c)	&& c!=k)
						return true;
				}
			}	
		}
	}
	return false;
}
// ------------