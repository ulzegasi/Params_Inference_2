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

inline double w_stg(int k, double T, double dt, double m_stg) {
	double res = sqrt(T*k/(dt*(k-1)*m_stg));
	return res;
};

void napa(vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &,
	int, int, int, int, int, double, double, double, double, double, double,
	vector<double> &, vector<double> &, int, vector<adouble> &,
	adept::Stack &);

// vector<double> &, vector<double> &, vector<double> &,