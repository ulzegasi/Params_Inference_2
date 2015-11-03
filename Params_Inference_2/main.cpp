//
// Parameter inference with SDE model 
//  
// Created by Simone Ulzega, Carlo Albert (October, 2015)
//
// This is the first parallel version of the code
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include "math.h"
#include <float.h>
#include <array>
#include <random>
#include <algorithm>
#include <functional>
#include <iomanip>
#include "arraysize.h"
#include "napa.h"
#include "vecfun.h"
#include "Vfun.h"
#include "adept_source.h"
#include "aVfun.h"
#include <omp.h>
//#include "timer.hpp" // On UNIX
using namespace std;
using adept::adouble;

// Useful directories
// const string dir  = "C:/Users/ulzegasi/Cpp/Params_Inference_2/";
const string dir2 = "C:/Users/ulzegasi/Cpp/Params_Inference_2/temp_data/";
const string dir3 = "C:/Users/ulzegasi/Cpp/Params_Inference_2/input_data/";

// Read time interval from file
bool tdat( vector<double> & t_limits )
{
	ifstream ifs(dir3+"t.dat");
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir3+"t.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		int range = 5002;
		string line;
		int pos = 0, ix = 0;
		string col1; double col2;
		
		while (getline(ifs, line) && ix < range) // read one line from ifs
		{
			if (ix == 1 || ix == (range-1))
			{
				istringstream iss(line); // access line as a stream
				iss >> col1 >> col2;
				t_limits[pos++]=col2;
			}
			++ix;
		}
		ifs.close();
		return true;
	}
}

// Read input data from file
bool indat( vector<double> & y )
{
	ifstream ifs( dir3 + "y_n10_K50_g02_s10_sinr.dat" );
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir3 + "y_n10_K50_g02_s10_sinr.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		string line;
		double yvalue;
		while (getline(ifs, line))   // read one line from ifs
		{
			istringstream iss(line); // access line as a stream
			while (iss >> yvalue)
				y.push_back(yvalue);
		}
		ifs.close();
		return true;
	}
}

// Read "true" system realization from file
bool Sdat( vector<double> & S )
{
	ifstream ifs( dir3 + "St_n10_K50_g02_s10_sinr.dat" );
	if (! ifs)
	{
		cerr << "\nInvalid file name or location "
			<< dir3 + "St_n10_K50_g02_s10_sinr.dat" << " ... Aborting...\n\n";
		return false;
	}
	else
	{	
		string line;
		double Svalue;
		while (getline(ifs, line)) // read one line from ifs
		{
			istringstream iss(line); // access line as a stream
			iss >> Svalue;
			S.push_back(Svalue);
		}
		ifs.close();
		return true;
	}
}

//// Generate a vector of normally distributed random numbers with mean mu (default 0) and standard deviation std (default 1)
//void nrand(vector<double> & vec, size_t seed, double mu = 0.0, double std = 1.0)
//{
//	static mt19937_64 engine(seed);
//	// static mt19937_64 engine(time(NULL));
//	static normal_distribution<double> norm_dist(mu,std);
//	for (int ix = 0; ix < vec.size(); ++ix)
//		vec[ix] = norm_dist(engine);
//}
//
//// Generate a random number from uniform ditribution between a and b (default 0 and 1) 
//double urand(size_t seed, double a = 0.0, double b = 1.0)
//{
//	static mt19937_64 engine(seed);
//	// static mt19937_64 engine(time(NULL));
//	static uniform_real_distribution<double> unif_dist(a,b);
//	double res = unif_dist(engine);
//	return res;
//}

int main()
{

	// ******************************************************* //
	// *******************  Preliminaries ******************** //
	// ******************************************************* //

	// **************  System parameters *************** //
	const int nparams = 2;
	const int n = 10;
	const int j = 30;
	const int N = int (n*j+1);

	// **************  Time parameters *************** //
	double T, dt = 0;
	vector<double> t_limits(2);
	vector<double> t(N);
	vector<int> ty(n+1);

	if (!tdat(t_limits))   // Read time points from file
		return -1;

	T  = t_limits.back()-t_limits.front();  // Total time
	dt = T/(double)(N-1);					// Time step

	t[0] = t_limits.front();				// Time points (N)
	for (int ix = 1; ix < N; ++ix)
		t[ix] = t[ix-1] + dt;

	for (int ix = 0; ix < n+1; ++ix)	   // Indexes of measurement points (= n+1 boundary beads)
		ty[ix] = ix*j+1;

	// **************  HMC parameters *************** //
	const int nsample_burnin = 0;         // Number of points in the MCMC
	int nsample_eff = 25000;
	size_t num_threads = 2;
	int nsample = nsample_eff + nsample_burnin;

	const double dtau = 0.25;  // MD time step
	const int n_napa = 3;      // Number of NAPA steps

	const double true_K = 50.0;   // Retention time
	const double true_gam = 0.2;  // Dimensionless noise parameter
	const double true_bet = sqrt(T*true_gam/true_K);
	const double sigma = 0.10;     // Measurement noise
	vector<double> true_theta(nparams);  // Parameters to be inferred (beta, tau)
	true_theta[0] = true_bet;
	true_theta[1] = true_gam;


	double K = 200.0;    // Initial state
	double gam = 0.5;
	double bet = sqrt(T*gam/K);                                   
	vector<double> theta(nparams);
	vector<double> * theta_pt = &theta;
	theta[0] = bet;
	theta[1] = gam;

	double K_min   = 0.0;   // Parameter limits
	double gam_min = 0.0;
	double K_max   = 1000.0;
	double gam_max = 5.0;

	// **************  Generate data  *************** // 
	vector<double> r(N);         // Input signal (rain)
	vector<double> lnr_der(N,0); // Log-derivative of the rain
	vector<double> y;
	vector<double> S;
	vector<double> long_t;
	vector<double> q_init(n+1);
	vector<double> bq(n+1);
	vector<double> q(N,0);
	
	for (int ix = 0; ix < N; ++ix)
		r[ix] = pow(sin(t[ix]/100),2.0) + 0.1;

	for (int ix = 0; ix < N-1; ++ix)
		lnr_der[ix] = (log(r[ix+1])-log(r[ix]))/dt;

	if (!indat(y))   // Read data points from file
		return -1;

	if (!Sdat(S))   // Read "true" system realization from file
		return -1;

	double tiny_dt = T/((S.size())-1);	// Time step for the "true" system realization 
	long_t.push_back(t[0]);				    // Time points (N)
	for (int ix = 1; ix < S.size(); ++ix)
		long_t.push_back(long_t[ix-1] + tiny_dt);

	for (int ix = 0; ix < n+1; ++ix)
	{
		q_init[ix] = (1/true_bet)*log(y[ix]/r[ty[ix]-1]);
		bq[ix] = true_bet*q_init[ix];
	}

	for (int s = 0; s < n; ++s)
	{
		double step = (q_init[s+1]-q_init[s])/j;
		q[ty[s]-1] = q_init[s];
		q[ty[s+1]-1] = q_init[s+1];
		for (int ix = 1; ix < j; ++ix)
			q[ty[s]-1+ix] = q[ty[s]-1] + ix*step;
	}

	// **************************************************************** //
	// ******************  Hamiltonian Monte Carlo  ******************* //
	// **************************************************************** //
	
	// **************  Transformations q -> u  *************** // 
    // (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

	vector<double> u(N);

	for (int s = 0; s <= (n-1); ++s)
	{
		u[s*j] = q[s*j];
		for (int k = 2; k <= j; ++k)
			// IMPORTANT: the transformation that should be applied is :
			// u[s*j+k-1] = q[s*j+k-1] - ( (k-1)*q[s*j+k] + q[s*j] )/k;
			// BUT: it is always = 0 when the{ q } between data points are linearly distributed
			// For staging beads as linear interpolation of data points, therefore
			u[s*j + k - 1] = 0.0;
	}
	u[N-1] = q[N-1];

	// **************  Init chains storage  *************** //
	vector< vector<double> > theta_sample;
	vector< vector<double> > u_sample;
	vector<double> energies;

	theta_sample.push_back(theta);
	u_sample.push_back(u);

	// **************  Init container for masses  *************** //
	vector<double> mp(nparams+N);
	vector<double> sqrt_mp(nparams+N);
	vector<double> * mp_pt = &mp;
	vector<double> * sqrt_mp_pt = &sqrt_mp;

	// **************  Containers for AD  *************** //
	// vector<double> force_old(nparams+N);
	// vector<double> force_new(nparams+N);
	// vector<double> dVN(N);
	/*vector<double> dVn(nparams+N);
	vector<double> dV1(nparams+N);*/
	// vector<double> dV_vec(nparams+N);
	// adept::Stack stack;
	// vector<adouble> x(nparams+N);

	// **************  Masses (burn-in)  *************** //
	double m_bdy = 10.0;     // m = m_q / dt
	double m_stg = 1.0;      // we assume m_q prop. to dt ==> m = costant
	double m_theta_burnin_value = 1.0;
	vector<double> m_theta(nparams,m_theta_burnin_value);

	for (int ix = 0; ix < nparams; ++ix)
		(*mp_pt)[ix] = m_theta[ix];
	for (int s = 1; s < n+1; ++s)
	{
		(*mp_pt)[nparams+(s-1)*j] = m_bdy;                       
		for (int k = 2; k < j+1; ++k) 
			(*mp_pt)[nparams+(s-1)*j+k-1] = m_stg;                   
	}
	(*mp_pt)[nparams+N-1] = m_bdy;  

	*sqrt_mp_pt = vsqrt(*mp_pt);
	
	// ************** HMC loops *************** //
	int reject_counter_burnin = 0;

	std::cout << "\nStarting HMC loops (burn-in)...\n---------------------------------\n";
	clock_t tinit = clock();
	// On UNIX:
	// timer tinit;
	// tinit.start();

	for (int counter = 1; counter <= nsample_burnin; ++counter)
		std::cout << "\nHMC loops here ... to be implemented ...\n";

	std::cout << "\nHMC loops (burn-in) terminated ...\n---------------------------------\n";

	// ************** Redefinition of masses (effective HMC loops) *************** //
	double m_bdy_burnin = m_bdy;
	double m_stg_burnin = m_stg;
	vector<double> m_theta_burnin = m_theta;

	m_bdy = 720.0;       // m = m_q / dt
	m_stg = 320.0;       // we assume m_q prop. to dt ==> m = costant     
	m_theta[0] = 150.0;  // refers to beta
	m_theta[1] = 130.0;  // refers to gamma

	for (int ix = 0; ix < nparams; ++ix)
		(*mp_pt)[ix] = m_theta[ix];
	for (int s = 1; s < n+1; ++s)
	{
		(*mp_pt)[nparams+(s-1)*j] = m_bdy;                       
		for (int k = 2; k < j+1; ++k) 
			(*mp_pt)[nparams+(s-1)*j+k-1] = m_stg;                   
	}
	(*mp_pt)[nparams+N-1] = m_bdy;

	*sqrt_mp_pt = vsqrt(*mp_pt);

	// ************** Effective HMC loops *************** //
	int reject_counter = 0;
	std::cout << "\nStarting effective HMC loops...\n---------------------------------\n";
	
	if ((nsample_eff%num_threads) != 0)
	{
		
		std::cout << "\nMC iterations " << nsample_eff
			<< " not compatible with number of threads " << num_threads << endl;

		nsample_eff = num_threads * round((double)nsample_eff / num_threads);

		std::cout << "Parameter nsample_eff will be set to " << nsample_eff << endl;

		nsample = nsample_eff + nsample_burnin;

	}

	// **************  Containers for execution times  *************** //
	vector<double> time_respa(nsample_eff);
	vector<double> time_respa_s(nsample_eff);
	vector<double> time_respa_f(nsample_eff);
	
	// **************  To jump out of the parallel region in case of NaN energy  *************** //
	bool abort = false;
	int abort_counter;
	int abort_thread;

	size_t seed = 18711221013;  // Can be used instead of time(NULL) for reproducibility

	#pragma omp parallel num_threads(num_threads)
	{
		#pragma omp master
		{
			std::cout << "Activated threads: " << omp_get_num_threads() << endl;
		}
		// **************  Init thread-local vector container for normal random numbers  *************** //
		vector<double> randvec(nparams + N);

		// **************  Init thread-local vector containers to store chains  *************** //
		vector< vector<double> > local_theta_sample;
		vector< vector<double> > local_u_sample;
		vector<double> local_energies;

		// **************  Init local vector container for momenta  *************** //
		vector<double> p; //(nparams+N);
		// vector<double> * p_pt = &p;

		// **************  Thread-local energies  *************** //
		double H_old, H_new;
		double accept_prob;

		// **************  Init thread-local containers for temp values  *************** //
		vector<double> theta_save(nparams);
		vector<double> u_save(N);

		// **************  Containers for AD  *************** //
		adept::Stack stack;
		vector<adouble> x(nparams + N);

		// **************  Theta and u must be thread-local!  *************** //
		vector<double> local_theta = theta;
		vector<double> local_u = u;

		// **************  Thread-local RNG  *************** //
		// mt19937_64 engine(time(NULL)*(omp_get_thread_num() + 1));
		mt19937_64 engine(seed*(omp_get_thread_num() + 1));
		std::function <double()> nrand = std::bind(normal_distribution<>(0, 1), std::ref(engine));
		std::function <double()> urand = std::bind(uniform_real_distribution<>(0, 1), std::ref(engine));

		int local_reject_counter = 0;

		#pragma omp flush(abort)
		// timer t0; t0.start();
		for (int counter = (nsample_burnin + 1); counter <= nsample; ++counter)
		{
			if (!abort)
			{
				// clock_t t0 = clock();
				// Sample momenta
				std::generate(randvec.begin(), randvec.end(), nrand);
				// nrand(randvec, seed*(omp_get_thread_num() + 1));
				// *p_pt = vtimes(*sqrt_mp_pt, randvec);
				p = vtimes(*sqrt_mp_pt, randvec);

				// Calculate energy
				// H_old = vsum(vdiv(vsquare(*p_pt), vtimes(2.0, *mp_pt))) + V_N(n, j, T, dt, u) + V_n(n, j, sigma, T, dt, bq, local_theta, u) + V_1(n, j, N, T, dt, lnr_der, local_theta, u);
				H_old = vsum(vdiv(vsquare(p), vtimes(2.0, *mp_pt))) + V_N(n, j, T, dt, local_u) + V_n(n, j, sigma, T, dt, bq, local_theta, local_u) + V_1(n, j, N, T, dt, lnr_der, local_theta, local_u);
				/*std::cout << H_old << endl;*/
				local_energies.push_back(H_old);

				if (isnan(H_old) != 0) {
					abort = true;
					#pragma omp flush (abort)
					abort_counter = counter;
					abort_thread = omp_get_thread_num();
				}

				// Save current state
				// (*theta_save_pt) = local_theta;
				// (*u_save_pt) = u;
				theta_save = local_theta;
				u_save = local_u;

				// MD Integration:
				clock_t t1 = clock();
				// timer t1; t1.start();

				for (int counter_napa = 1; counter_napa <= n_napa; ++counter_napa)
				{
					napa(local_theta, local_u, p, mp, bq, lnr_der, counter, n, j, N, nparams, sigma, T, dt, dtau, m_stg, m_bdy,
						time_respa_f, time_respa_s, nsample_burnin, x, stack);
					// force_old, force_new, dV_vec, 
				}
				#pragma omp master
				{
					time_respa[counter - nsample_burnin - 1] = ((float)(clock() - t1) / CLOCKS_PER_SEC);
					// t1.stop();
					// time_respa[counter - nsample_burnin - 1] = t1.get_timing();
				}

				// Calculate energy of proposal state => Metropolis accept/reject

				/*H_new = vsum(vdiv(vsquare(*p_pt), vtimes(2.0, *mp_pt))) + V_N(n, j, T, dt, u) + V_n(n, j, sigma, T, dt, bq, local_theta, u) + V_1(n, j, N, T, dt, lnr_der, local_theta, u);
				std::cout << setprecision(16) << H_new << endl;*/

				if (local_theta[1] <= gam_min || local_theta[1] >= gam_max || (T*local_theta[1] / pow(local_theta[0], 2) <= K_min) || (T*local_theta[1] / pow(local_theta[0], 2) >= K_max))
				{
					// local_theta = (*theta_save_pt); u = (*u_save_pt);
					local_theta = theta_save; local_u = u_save;
					local_reject_counter += 1;
				}
				else if (find_if(local_u.begin(), local_u.end(), [&](double el){return (el > 10.0); }) != local_u.end() ||
					find_if(local_u.begin(), local_u.end(), [&](double el){return (el < -10.0); }) != local_u.end())
				{
					// local_theta = (*theta_save_pt); u = (*u_save_pt);
					local_theta = theta_save; local_u = u_save;
					local_reject_counter += 1;
				}
				else
				{
					// H_new = vsum(vdiv(vsquare(*p_pt), vtimes(2.0, *mp_pt))) + V_N(n, j, T, dt, u) + V_n(n, j, sigma, T, dt, bq, theta, u) + V_1(n, j, N, T, dt, lnr_der, local_theta, u);
					H_new = vsum(vdiv(vsquare(p), vtimes(2.0, *mp_pt))) + V_N(n, j, T, dt, local_u) + V_n(n, j, sigma, T, dt, bq, local_theta, local_u) + V_1(n, j, N, T, dt, lnr_der, local_theta, local_u);
					accept_prob = min(1.0, exp(H_old - H_new));
					if (urand() > accept_prob)
					{
						// local_theta = (*theta_save_pt); u = (*u_save_pt);
						local_theta = theta_save; local_u = u_save;
						local_reject_counter += 1;
					}
				}

				local_theta_sample.push_back(local_theta);
				local_u_sample.push_back(local_u);

				#pragma omp critical (output)
				{
					if (counter % 1000 == 0)
					{
						// t0.stop();
						std::cout << "\n" << counter << " loops in thread " << omp_get_thread_num() << " completed in "
							<< ((float)(clock() - tinit) / CLOCKS_PER_SEC) << " seconds\n";
						//  << t0.get_timing() << " seconds\n";
					}
						
				}
			}  // end of if(!abort)
		}  // end of for loop

		#pragma omp barrier  // wait for all threads

		#pragma omp critical (writing_chains)
		{
			for (vector< vector<double> >::iterator row_it = local_theta_sample.begin(); 
				row_it != local_theta_sample.end(); ++row_it)
			{
				theta_sample.push_back(*row_it);
			}
				
			for (vector< vector<double> >::iterator row_it = local_u_sample.begin();
				row_it != local_u_sample.end(); ++row_it)
			{
				u_sample.push_back(*row_it);
			}
			
			for (vector<double>::iterator it = local_energies.begin();
				it != local_energies.end(); ++it)
			{
				energies.push_back(*it);
			}

			reject_counter += local_reject_counter;

		} // end of critical section (writing chains)

	}
	// END OF PARALLEL REGION
	if (abort)
	{
		cerr << "\n\nIteration " << abort_counter << ", thread " << abort_thread
			<< " --> energy values diverged...\n\n";
		return -1;
	}

	// ***************************************************************** //
	// **********************  End of HMC loops  *********************** //
	// ***************************************************************** //
	std::cout << "\n\n";
	std::cout << "Number of active threads: " << num_threads << endl;
	std::cout << "Run completed in " << ((float)(clock()-tinit)/CLOCKS_PER_SEC) << " seconds\n";
	// tinit.stop();
	// std::cout << "Run completed in " << tinit.get_timing() << " seconds\n";
	std::cout << "NAPA cycles in " << vsum(time_respa) << " seconds\n";
	std::cout << "Slow NAPA in " << vsum(time_respa_s) << " seconds\n";
	std::cout << "Fast NAPA in " << vsum(time_respa_f) << " seconds\n"; 
	
	// **************  Back transformations u -> q , y  *************** // 
    // (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)
	
	clock_t back_transform = clock();

	vector< vector<double> > qs(u_sample.size(), vector<double>(u_sample[0].size()));
	vector< vector<double> > predy(u_sample.size(), vector<double>(u_sample[0].size()));

	int sample_size = num_threads*nsample_eff;

	/*for (int sample_ind = 0; sample_ind < (nsample+1); ++sample_ind)
	{
		for (int s = 1; s <= n; ++s)
		{
			qs[sample_ind][(s-1)*j] = u_sample[sample_ind][(s-1)*j];
			for (int k = 2; k <= j; ++k)
			{
				qs[sample_ind][(s-1)*j+k-1] = ((double)(j-k+1)/j)*u_sample[sample_ind][(s-1)*j];
				for (int l = k; l <= j+1; ++l)
				{
					qs[sample_ind][(s-1)*j+k-1] += ((double)(k-1)/(l-1))*u_sample[sample_ind][(s-1)*j+l-1];
				}
			}
		}
		qs[sample_ind][N - 1] = u_sample[sample_ind][N - 1];
	}*/

	// ABOVE, slow method (cache thrashing?)
	// BELOW, much faster
	for (int sample_ind = 0; sample_ind < (sample_size + 1); ++sample_ind)
	{
		for (int s = 1; s <= n; ++s)
		{
			qs[sample_ind][(s - 1)*j] = u_sample[sample_ind][(s - 1)*j];
		}
		qs[sample_ind][N - 1] = u_sample[sample_ind][N - 1];
	}

	for (int sample_ind = 0; sample_ind < (sample_size + 1); ++sample_ind)
	{
		for (int s = n; s > 0; --s)
		{
			for (int k = j; k > 1; --k)
			{
				qs[sample_ind][(s - 1)*j + k - 1] = u_sample[sample_ind][(s - 1)*j + k - 1] +
					((double)(k - 1.0) / k)*qs[sample_ind][(s - 1)*j + k] +
					((double)(1.0 / k))*u_sample[sample_ind][(s - 1)*j];
			}
		}
	}

	for (int sample_ind = 0; sample_ind < (sample_size + 1); ++sample_ind)
	{
		for (int ind = 0; ind < N; ++ind)
		{
			predy[sample_ind][ind] = r[ind]*exp(theta_sample[sample_ind][0]*qs[sample_ind][ind]); // y = r*exp(beta*q)
		}
	}

	std::cout << "\nBack transformations u -> q,y completed in " 
		<< ((float)(clock()-back_transform)/CLOCKS_PER_SEC) << " seconds\n";


	// ******************************************************* //
	// *******************  Save results ******************** //
	// ******************************************************* //

	const string fname = "_omptest";

	// *******************  Parameters ******************** //
	string param_names[] = {"N", "j", "n", "t[1]", "dt", "nparams", "true_K", "true_gam", "sigma", "K", "gam", 
	"nsample_burnin", "nsample_eff", "m_bdy_burnin", "m_bdy", "m_theta_burnin", "m_theta_bet", 
	"m_theta_gam", "m_stg_burnin", "m_stg", "dtau", "n_napa", "num_threads"};

	double param_values[] = {N, j, n, t[0], dt, nparams, true_K, true_gam, sigma, K, gam, 
	nsample_burnin, nsample_eff, m_bdy_burnin, m_bdy, m_theta_burnin[0], m_theta[0], m_theta[1], 
	m_stg_burnin, m_stg, dtau, n_napa, num_threads};

	ofstream ofs_params(dir2 + "params" + fname + ".dat");
	if (! ofs_params)
	{
		cerr << "\nInvalid output file name or location " << dir2+"params_"+fname+".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int ix = 0; ix < array_size(param_names); ++ix)
			ofs_params << param_names[ix] << "\t" << param_values[ix] << endl;
	}
	ofs_params.close();

	// *******************  Last state ******************** //
	ofstream ofs_lastq(dir2 + "last_qs" + fname + ".dat");
	if (! ofs_lastq)
	{
		cerr << "\nInvalid output file name or location " << dir2+"last_qs"+fname+".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int ix = 0; ix < qs.size(); ++ix)
			ofs_lastq << qs[ix].back() << endl;
	}
	ofs_lastq.close();

	// *******************  u chains ******************** //
	// vector< vector<double> > u_chains(nsample+2, vector<double>(4));
	int ind_vec[4];
	ind_vec[0] = ty[floor((double)ty.size()/2) - 1];
	ind_vec[1] = ty[floor((double)ty.size()/2) - 1] + floor((double)(j-1)/2);
	ind_vec[2] = ty[floor((double)ty.size()/2) - 1] + j;
	ind_vec[3] = ty[floor((double)ty.size()/2) -1 ] + j + floor((double)(j-1)/2);

	ofstream ofs_uchains(dir2 + "u_chains" + fname + ".dat");
	if (! ofs_uchains)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "u_chains" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int col = 0; col < 4; ++col)
		{
			ofs_uchains << ind_vec[col] << " ";
		}
		ofs_uchains << endl;
		for (int row = 0; row < sample_size + 1; ++row)
		{
			for (int col = 0; col < 4; ++col)
			{
				ofs_uchains << u_sample[row][ind_vec[col] - 1] << " ";
			}
			ofs_uchains << endl;
		}
	}
	ofs_uchains.close();

	// *******************  y chains for spaghetti ******************** //
	int red_range = min(200,sample_size+1);
	double step = (double)sample_size / (red_range - 1);
	vector<int> red_ind;
	for (int ix = 0; ix < red_range; ++ix)
		red_ind.push_back((int)(1 + ix*step));

	ofstream ofs_predy(dir2 + "predy" + fname + ".dat");
	if (! ofs_predy)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "predy" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int row = 0; row < red_ind.size(); ++row)
		{
			for (int col = 0; col < N; ++col)
			{
				ofs_predy << predy[red_ind[row]-1][col] << " ";
			}
			ofs_predy << endl;
		}
	}
	ofs_predy.close();

	// *******************  q chains for spaghetti ******************** //
	ofstream ofs_predq(dir2 + "predq" + fname + ".dat");
	if (! ofs_predq)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "predq" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int row = 0; row < red_ind.size(); ++row)
		{
			for (int col = 0; col < N; ++col)
			{
				ofs_predq << qs[red_ind[row]-1][col] << " ";
			}
			ofs_predq << endl;
		}
	}
	ofs_predq.close();

	// *******************  theta ******************** //
	ofstream ofs_thetas(dir2 + "thetas" + fname + ".dat");
	if (! ofs_thetas)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "thetas" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int row = 0; row < sample_size + 1; ++row)
		{
			ofs_thetas << T * theta_sample[row][1] / pow(theta_sample[row][0],2) 
				<< "\t" << theta_sample[row][1] << endl;
		}
	}
	ofs_thetas.close();

	// *******************  energies ******************** //
	ofstream ofs_energies(dir2 + "energies" + fname + ".dat");
	if (! ofs_energies)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "energies" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		for (int row = 0; row < sample_size; ++row)
			ofs_energies << energies[row] << endl;
	}
	ofs_energies.close();

	// *******************  reject rates ******************** //
	ofstream ofs_reject(dir2 + "reject" + fname + ".dat");
	if (! ofs_reject)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "reject" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		ofs_reject << (double)reject_counter_burnin/nsample_burnin*100 << " "
			<< (double)reject_counter/sample_size*100 << endl;
	}
	ofs_reject.close();
	
	// *******************  I/O data ******************** //
	step = (double)(S.size()-1) / (N - 1);
	vector<double> red_S;
	for (int ix = 0; ix < N; ++ix)
		red_S.push_back(S[(int)(1 + ix*step)-1]);

	ofstream ofs_iodata(dir2 + "iodata" + fname + ".dat");
	if (! ofs_iodata)
	{
		cerr << "\nInvalid output file name or location " << dir2 + "iodata" + fname + ".dat" << " ... Aborting...\n\n";
	}
	else
	{
		ofs_iodata << "System_realization" << " ";
		for (int ix = 0; ix < red_S.size(); ++ix)
			ofs_iodata << red_S[ix] << " ";
		ofs_iodata << endl;

		ofs_iodata << "Rain_input" << " ";
		for (int ix = 0; ix < r.size(); ++ix)
			ofs_iodata << r[ix] << " ";
		ofs_iodata << endl;

		ofs_iodata << "Output_flow" << " ";
		for (int ix = 0; ix < y.size(); ++ix)
			ofs_iodata << y[ix] << " ";
		ofs_iodata << endl;
	}
	ofs_iodata.close();
	
	std::cout << "\n\n";
	return 0;
}