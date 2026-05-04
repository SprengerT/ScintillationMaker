////////////////////////////////////////////////////////////////////////

// g++ -Wall -O2 -fopenmp -shared -Wl,-soname,ScintillationMaker -o ScintillationMaker.so -fPIC ScintillationMaker.cpp

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;



extern "C"
void  simulate_SimpleScreen (double *E_real, double *E_im, int N_t, int N_nu, int N_th, double *nu, double *t, double *thx, double *thy, double *mu, double *ph )
{
	#define  REAL(i_t,i_nu)  E_real[(i_t)*N_nu + (i_nu)]
	#define  IMAG(i_t,i_nu)  E_im[(i_t)*N_nu + (i_nu)]
	
	#pragma omp parallel for
	for (int i_nu = 0; i_nu < N_nu; i_nu++){
		double Phase;
        const double v_nu = nu[i_nu];
        for (int i_t = 0; i_t < N_t; i_t++){
            double realsum = 0.;
            double imagsum = 0.;
            const double v_t = t[i_t];
            #pragma omp simd reduction(+:realsum,imagsum)
			for (int i_th = 0; i_th < N_th; i_th++){
                Phase = ph[i_th] + v_nu*(thy[i_th]*thy[i_th] + thy[i_th]*thy[i_th]) + v_t*thx[i_th];
                const double v_mu = mu[i_th];
				realsum += v_mu*cos(Phase);
				imagsum += v_mu*sin(Phase);
			}
			REAL(i_t,i_nu) = realsum;
            IMAG(i_t,i_nu) = imagsum;
		}
	}
	
	#undef  REAL
	#undef  IMAG
	#undef  TAU
} 
