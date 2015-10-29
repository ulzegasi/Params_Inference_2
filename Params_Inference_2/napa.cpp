#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <float.h>
#include <array>
#include <random>
#include <algorithm>
#include <functional>
#include <omp.h>
#include "napa.h"
#include "vecfun.h"
#include "Vfun.h"
#include "aVfun.h"
#include "adept.h"
using namespace std;


void napa(vector<double> & theta, vector<double> & u, vector<double> & p, vector<double> & mp, vector<double> & bq, vector<double> & lnr_der,
	int counter, int n, int j, int N, int nparams, double sigma, double T, double dt, double dtau, double m_stg, double m_bdy, 
	vector<double> & time_respa_f, vector<double> & time_respa_s, int nsample_burnin, vector<adouble> & x,
	adept::Stack & stack)
{
	// vector<double> & force_old, vector<double> & force_new, vector<double> & dV_vec, 
	
	static vector<double> force_old(nparams+N);
	static vector<double> force_new(nparams+N);

	// Fast outer propagator (V_N), step dtau/2 
	clock_t timef = clock();
    for (int s = 1; s <= n; ++s)
	{
        for (int k = 2; k <=j; ++k)
		{
            double u_old = u[(s-1)*j+k-1];
            u[(s-1)*j+k-1] = u[(s-1)*j+k-1]*cos(w_stg(k, T, dt, m_stg)*dtau/2.0) 
				+ p[nparams+(s-1)*j+k-1]*sin(w_stg(k, T, dt, m_stg)*dtau/2.0) / (m_stg*w_stg(k, T, dt, m_stg));
            p[nparams+(s-1)*j+k-1] = p[nparams+(s-1)*j+k-1]*cos(w_stg(k, T, dt, m_stg)*dtau/2.0) 
				- m_stg*w_stg(k, T, dt, m_stg)*u_old*sin(w_stg(k, T, dt, m_stg)*dtau/2.0);
		}
	}
	#pragma omp master
	{
		time_respa_f[counter - nsample_burnin - 1] += ((float)(clock() - timef) / CLOCKS_PER_SEC);
	}

	// Slow inner propagator (V_n, V_1), step dtau
    clock_t time1 = clock();

	dV_fun(stack, n, j, N, sigma, T, dt, bq, lnr_der, theta, u, x, force_old);
	// force_old = vtimes(-1.0,force_old);  // NOT NECESSARY: sign changed already in the definition of aV_n_1 

	for (int s = 1; s <= (n+1); ++s)
		u[(s-1)*j] += dtau * ( p[nparams+(s-1)*j]  + (dtau/2.0) * force_old[nparams+(s-1)*j] ) / m_bdy;
	for (int ix = 1; ix <= nparams; ++ix)
		theta[ix-1] += dtau * ( p[ix-1]  + (dtau/2.0) * force_old[ix-1] ) / mp[ix-1];
	
	dV_fun(stack, n, j, N, sigma, T, dt, bq, lnr_der, theta, u, x, force_new);
	// force_new = vtimes(-1.0,force_new);  // NOT NECESSARY: sign changed already in the definition of aV_n_1 

	for (int ix = 1; ix <= nparams+N; ++ix)
		p[ix - 1] += (dtau / 2)*(force_old[ix - 1] + force_new[ix - 1]);
	#pragma omp master
	{
		time_respa_s[counter - nsample_burnin - 1] += ((float)(clock() - time1) / CLOCKS_PER_SEC);
	}
	
	// Again, fast outer propagator (V_N), step dtau/2 
	timef = clock();
    for (int s = 1; s <= n; ++s)
	{
        for (int k = 2; k <=j; ++k)
		{
            double u_old = u[(s-1)*j+k-1];
            u[(s-1)*j+k-1] = u[(s-1)*j+k-1]*cos(w_stg(k, T, dt, m_stg)*dtau/2.0) 
				+ p[nparams+(s-1)*j+k-1]*sin(w_stg(k, T, dt, m_stg)*dtau/2.0) / (m_stg*w_stg(k, T, dt, m_stg));
            p[nparams+(s-1)*j+k-1] = p[nparams+(s-1)*j+k-1]*cos(w_stg(k, T, dt, m_stg)*dtau/2.0) 
				- m_stg*w_stg(k, T, dt, m_stg)*u_old*sin(w_stg(k, T, dt, m_stg)*dtau/2.0);
		}
	}
	#pragma omp master
	{
		time_respa_f[counter - nsample_burnin - 1] += ((float)(clock() - timef) / CLOCKS_PER_SEC);
	}
    
}