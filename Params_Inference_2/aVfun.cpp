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
#include "adept.h"
using namespace std;
using adept::adouble;

// The sign of out in each function aV_* must be changed !
//  => F = dV_vec

// ***************************************************************** //
// ***************************Function***************************** //
//adouble aV_N(int n, int j, double T, double dt, const vector<adouble> & u)
//{
//	adouble out = 0.0;
//	for (int s = 1; s <= n; ++s)
//	{
//		for (int k = 2; k <= j; ++k)
//		{
//			out += 0.5*(T/dt)*(double)k/(k-1)*pow(u[(s-1)*j+k-1],2);
//		}
//	}
//	return out;
//};
// ***************************Jacobian***************************** //
//void dV_N(adept::Stack & stack, int n, int j, double T, double dt, const vector<double> & u, vector<double> & dy_dx)
//{
//	int nx = u.size();
//	vector<adouble> x(nx);
//	adept::set_values(&x[0],nx,&u[0]);
//	stack.new_recording();
//	adouble y = aV_N(n, j, T, dt, x);
//	y.set_gradient(1.0);
//	stack.compute_adjoint();
//	adept::get_gradients(&x[0], nx, &dy_dx[0]);
//	/*stack.independent(&x[0],nx);
//	stack.dependent(y);
//	stack.jacobian(&dy_dx[0]);*/
//}

// ***************************************************************** //
// ***************************Function***************************** //
//double V_n(int n, int j, double sigma, double T, double dt,
//	const vector<double> & bq, const vector<double> & theta, const vector<double> & u)
//{
//	double beta  = theta[0]; // double gamma = theta[1];
//	double out = pow((bq[n]-beta*u[n*j]),2) / (2.0*pow(sigma,2));
//
//	for (int s = 1; s <= n; ++s)
//        out += 0.5*(T/dt)*(1.0/(double)j)*pow((u[(s-1)*j]-u[s*j]),2) + 
//		(pow((bq[s-1]-beta*u[(s-1)*j]),2)/(2.0*pow(sigma,2))); 
//
//    return out;
//};

// Above is the original function V_n. Below, vectors theta and u are merged into a vector x
// All elements of x corresponding to elements of u have their original indexes plus ntheta = theta.size()
// In this way, the first ntheta elements of x correspond to theta, the remaining N elements correspond to u

//adouble aV_n(int n, int j, double sigma, double T, double dt,
//	const vector<double> & bq, const vector<adouble> & x, int ntheta)
//{
//	adouble out  = pow((bq[n]-x[0]*x[n*j+ntheta]),2) / (2.0*pow(sigma,2));
//	for (int s = 1; s <= n; ++s)
//        out += 0.5*(T/dt)*(1.0/(double)j)*pow((x[(s-1)*j+ntheta]-x[s*j+ntheta]),2) + 
//		(pow((bq[s-1]-x[0]*x[(s-1)*j+ntheta]),2)/(2.0*pow(sigma,2))); 
//    return out;
//};
// ***************************Jacobian***************************** //
//void dV_n(adept::Stack & stack, int n, int j, double sigma, double T, double dt, 
//	const vector<double> & bq, const vector<double> & theta, const vector<double> & u, vector<adouble> & x, vector<double> & dy_dx)
//{
//	int ntheta = theta.size();
//	int nu = u.size();
//	int nx = ntheta + nu;
//	adept::set_values(&x[0],ntheta,&theta[0]);  // Construct a vector of adoubles where the first ntheta elements are theta
//	adept::set_values(&x[ntheta],nu,&u[0]);     // and the remaining nu elements are u
//	stack.new_recording();
//	adouble y = aV_n(n, j, sigma, T, dt, bq, x, ntheta);
//	y.set_gradient(1.0);
//	stack.compute_adjoint();
//	adept::get_gradients(&x[0], nx, &dy_dx[0]);
//	/*stack.independent(&x[0],nx);
//	stack.dependent(y);
//	stack.jacobian(&dy_dx[0]);*/
//}

// ***************************************************************** //
// ***************************Function***************************** //
//double V_1(int n, int j, int N, double T, double dt, 
//	const vector<double> & lnr_der, const vector<double> & theta, const vector<double> & u)
//{
//	double beta  = theta[0]; 
//	double gamma = theta[1];
//	double rho   = (2.0+gamma)*beta/(2.0*gamma);
//
//	double out = (1.0/gamma)*exp(-beta*u[N-1]) + u[N-1]*((T/beta)*lnr_der[N-2]+rho) - 
//		(1.0/gamma)*exp(-beta*u[0]) - u[0]* ((T/beta)*lnr_der[0]+rho);
//
//	for (int s = 1; s <= (n-1); ++s)
//		out += (dt/T)*0.5*pow( (((T/beta) * lnr_der[s*j] + rho )-(beta/gamma)*exp(-beta*u[s*j])), 2 ) -
//			(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*u[s*j]) -
//			u[s*j] *(T/beta) * (lnr_der[s*j] - lnr_der[s*j-1]);
//
//	for (int s = 1; s <= n; ++s)
//	{
//		for (int k = 2; k <= j; ++k)
//		{
//			double tmp = ((double)(j-k+1))*u[(s-1)*j]/(double)j;
//			for (int l = k; l<=(j+1); ++l)
//                tmp += ((double)(k-1))*u[(s-1)*j+l-1]/((double)(l-1));
//			out += (dt/T)*0.5* pow( (((T/beta)*lnr_der[(s-1)*j+k-1]+rho )-(beta/gamma)*exp(-beta*tmp)) , 2) -
//				(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*tmp) -
//				tmp * (T/beta) * (lnr_der[(s-1)*j+k-1] - lnr_der[(s-1)*j+k-2]);
//		}
//	}
//
//	return out;
//}

// Above is the original function V_1. Below, vectors theta and u are merged into a vector x
// All elements of x corresponding to elements of u have their original indexes plus ntheta = theta.size()
// In this way, the first ntheta elements of x correspond to theta, the remaining N elements correspond to u

//adouble aV_1(int n, int j, int N, double T, double dt,
//	const vector<double> & lnr_der, const vector<adouble> & x, int ntheta)
//{
//	adouble beta  = x[0]; 
//	adouble gamma = x[1];
//	adouble rho   = (2.0+x[1])*x[0]/(2.0*x[1]);
//
//	adouble out = (1.0/gamma)*exp(-beta*x[N-1+ntheta]) + x[N-1+ntheta]*((T/beta)*lnr_der[N-2]+rho) - 
//		(1.0/gamma)*exp(-beta*x[ntheta]) - x[ntheta]* ((T/beta)*lnr_der[0]+rho);
//
//	for (int s = 1; s <= (n-1); ++s)
//		out += (dt/T)*0.5*pow( (((T/beta) * lnr_der[s*j] + rho )-(beta/gamma)*exp(-beta*x[s*j+ntheta])), 2 ) -
//			(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*x[s*j+ntheta]) -
//			x[s*j+ntheta] *(T/beta) * (lnr_der[s*j] - lnr_der[s*j-1]);
//
//	for (int s = 1; s <= n; ++s)
//	{
//		for (int k = 2; k <= j; ++k)
//		{
//			adouble tmp = ((double)(j-k+1))*x[(s-1)*j+ntheta]/(double)j;
//			for (int l = k; l<=(j+1); ++l)
//                tmp += ((double)(k-1))*x[(s-1)*j+l-1+ntheta]/((double)(l-1));
//			out += (dt/T)*0.5* pow( (((T/beta)*lnr_der[(s-1)*j+k-1]+rho )-(beta/gamma)*exp(-beta*tmp)) , 2) -
//				(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*tmp) -
//				tmp * (T/beta) * (lnr_der[(s-1)*j+k-1] - lnr_der[(s-1)*j+k-2]);
//		}
//	}
//
//	return out;
//};
// ***************************Jacobian***************************** //
//void dV_1(adept::Stack & stack, int n, int j, int N, double T, double dt, 
//	const vector<double> & lnr_der, const vector<double> & theta, const vector<double> & u, vector<adouble> & x, vector<double> & dy_dx)
//{
//	int ntheta = theta.size();
//	int nu = u.size();
//	int nx = ntheta + nu;
//	adept::set_values(&x[0],ntheta,&theta[0]);  // Construct a vector of adoubles where the first ntheta elements are theta
//	adept::set_values(&x[ntheta],nu,&u[0]);     // and the remaining nu elements are u
//	stack.new_recording();
//	adouble y = aV_1(n, j, N, T, dt, lnr_der, x, ntheta);
//	y.set_gradient(1.0);
//	stack.compute_adjoint();
//	adept::get_gradients(&x[0], nx, &dy_dx[0]);
//	/*stack.independent(&x[0],nx);
//	stack.dependent(y);
//	stack.jacobian(&dy_dx[0]);*/
//}

// ********************************************************************************************* //
// ********************************************************************************************* //
// *************************** Merge V_n and V_1 ***************************** //
// ********************************************************************************************* //
// ********************************************************************************************* //
adouble aV_n_1(int n, int j, int N, double sigma, double T, double dt,
	const vector<double> & bq, const vector<double> & lnr_der, const vector<adouble> & x, int ntheta)
{
	adouble beta  = x[0]; 
	adouble gamma = x[1];
	adouble rho   = (2.0+x[1])*x[0]/(2.0*x[1]);

	adouble out = (1.0/gamma)*exp(-beta*x[N-1+ntheta]) + x[N-1+ntheta]*((T/beta)*lnr_der[N-2]+rho) - 
		(1.0/gamma)*exp(-beta*x[ntheta]) - x[ntheta]* ((T/beta)*lnr_der[0]+rho);

	for (int s = 1; s <= (n-1); ++s)
		out += (dt/T)*0.5*pow( (((T/beta) * lnr_der[s*j] + rho )-(beta/gamma)*exp(-beta*x[s*j+ntheta])), 2 ) -
			(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*x[s*j+ntheta]) -
			x[s*j+ntheta] *(T/beta) * (lnr_der[s*j] - lnr_der[s*j-1]);

	for (int s = 1; s <= n; ++s)
	{
		for (int k = 2; k <= j; ++k)
		{
			adouble tmp = ((double)(j-k+1))*x[(s-1)*j+ntheta]/(double)j;
			for (int l = k; l<=(j+1); ++l)
                tmp += ((double)(k-1))*x[(s-1)*j+l-1+ntheta]/((double)(l-1));
			out += (dt/T)*0.5* pow( (((T/beta)*lnr_der[(s-1)*j+k-1]+rho )-(beta/gamma)*exp(-beta*tmp)) , 2) -
				(dt/T)*0.5*(pow(beta,2)/gamma)*exp(-beta*tmp) -
				tmp * (T/beta) * (lnr_der[(s-1)*j+k-1] - lnr_der[(s-1)*j+k-2]);
		}
	}

	out += pow((bq[n]-beta*x[n*j+ntheta]),2) / (2.0*pow(sigma,2));
	for (int s = 1; s <= n; ++s)
        out += 0.5*(T/dt)*(1.0/(double)j)*pow((x[(s-1)*j+ntheta]-x[s*j+ntheta]),2) + 
		(pow((bq[s-1]-beta*x[(s-1)*j+ntheta]),2)/(2.0*pow(sigma,2)));
	
	out *= (-1); // Change sign of the output => F = (-1)*dV = d(-1*V)

	return out;
};

// ***************************Full Jacobian***************************** //
void dV_fun(adept::Stack & stack, int n, int j, int N, double sigma, double T, double dt, 
	const vector<double> & bq, const vector<double> & lnr_der, const vector<double> & theta,
	const vector<double> & u, vector<adouble> & x, vector<double> & dy_dx)
{
	int ntheta = theta.size();
	int nu = u.size();
	int nx = ntheta + nu;
	adept::set_values(&x[0],ntheta,&theta[0]);  // Construct a vector of adoubles where the first ntheta elements are theta
	adept::set_values(&x[ntheta],nu,&u[0]);     // and the remaining nu elements are u
	stack.new_recording();
	(aV_n_1(n, j, N, sigma, T, dt, bq, lnr_der, x, ntheta)).set_gradient(1.0);
	stack.compute_adjoint();
	adept::get_gradients(&x[0], nx, &dy_dx[0]);
	/*stack.independent(&x[0],nx); stack.dependent(y); stack.jacobian(&dy_dx[0]);*/
}