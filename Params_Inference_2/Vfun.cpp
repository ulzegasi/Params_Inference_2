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
using namespace std;

// ***************************************************************** ///
// ***************************************************************** //
double V_N(int n, int j, double T, double dt, const vector<double> & u)
{
	double out = 0.0;
	for (int s = 1; s <= n; ++s)
	{
		for (int k = 2; k <= j; ++k)
		{
			out += 0.5*(T/dt)*(double)k/(k-1)*pow(u[(s-1)*j+k-1],2);
		}
	}
	return out;
};
// ***************************************************************** ///
// ***************************************************************** //
double V_n(int n, int j, double sigma, double T, double dt,
	const vector<double> & bq, const vector<double> & theta, const vector<double> & u)
{
	double beta  = theta[0]; // double gamma = theta[1];
	double out = pow((bq[n]-beta*u[n*j]),2) / (2.0*pow(sigma,2));

	for (int s = 1; s <= n; ++s)
        out += 0.5*(T/dt)*(1.0/(double)j)*pow((u[(s-1)*j]-u[s*j]),2) + 
		(pow((bq[s-1]-beta*u[(s-1)*j]),2)/(2.0*pow(sigma,2))); 

    return out;

};
// ***************************************************************** ///
// ***************************************************************** //
double V_1(int n, int j, int N, double T, double dt, 
	const vector<double> & lnr_der, const vector<double> & theta, const vector<double> & u)
{
	double beta  = theta[0]; 
	double gamma = theta[1];
	double rho   = (2.0+gamma)*beta/(2.0*gamma);

	double out = (1.0/gamma)*exp(-beta*u[N-1]) + u[N-1]*((T/beta)*lnr_der[N-2]+rho) - 
		(1.0/gamma)*exp(-beta*u[0]) - u[0]* ((T/beta)*lnr_der[0]+rho);

	for (int s = 1; s <= (n-1); ++s)
		out += (dt/T)*0.5*pow( (((T/beta) * lnr_der[s*j] + rho )-(beta/gamma)*exp(-beta*u[s*j])), 2 ) -
			(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*u[s*j]) -
			u[s*j] *(T/beta) * (lnr_der[s*j] - lnr_der[s*j-1]);

	for (int s = 1; s <= n; ++s)
	{
		for (int k = 2; k <= j; ++k)
		{
			double tmp = ((double)(j-k+1))*u[(s-1)*j]/(double)j;
			for (int l = k; l<=(j+1); ++l)
                tmp += ((double)(k-1))*u[(s-1)*j+l-1]/((double)(l-1));
			out += (dt/T)*0.5* pow( (((T/beta)*lnr_der[(s-1)*j+k-1]+rho )-(beta/gamma)*exp(-beta*tmp)) , 2) -
				(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*tmp) -
				tmp * (T/beta) * (lnr_der[(s-1)*j+k-1] - lnr_der[(s-1)*j+k-2]);
		}
	}

	return out;
}