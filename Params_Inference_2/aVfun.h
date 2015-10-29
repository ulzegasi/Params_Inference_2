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

adouble aV_N(int, int, double, double, const vector<adouble> &);
void dV_N(adept::Stack &, int, int, double, double, const vector<double> &, vector<double> &);

adouble aV_n(int, int, double, double, double, const vector<double> &, const vector<adouble> &, int);
void dV_n(adept::Stack & , int, int, double, double, double,
	const vector<double> &, const vector<double> &, const vector<double> &, vector<adouble> &, vector<double> &);

adouble aV_1(int, int, int, double, double, const vector<double> &, const vector<adouble> &, int);
void dV_1(adept::Stack &, int, int, int, double, double, 
	const vector<double> &, const vector<double> &, const vector<double> &, vector<adouble> &, vector<double> &);

adouble aV_n_1(int, int, int, double, double, double,
	const vector<double> &, const vector<double> &, const vector<adouble> &, int);
void dV_fun(adept::Stack &, int, int, int, double, double, double, 
	const vector<double> &, const vector<double> &, const vector<double> &,
	const vector<double> &, vector<adouble> &, vector<double> &);