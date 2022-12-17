/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file EpiIndicators.h
 * @brief functions for the EpiIndicators method to compute the
 *        functional relationship between 2 epidemiological indicators
 *
 * @author Luis Alvarez <lalvarez.mat@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <string.h>
#include <time.h>
#include <vector>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/// OPTIMIZATION OF THE RATIO AND SHIFT BETWEEN INDICATORS USING DIFFERENT INITIALIZATIONS OF s(t)
double shift_and_ratio_optimization(
    vector<string> &dates,vector<double> &c,vector<double> &d,vector<double> &r,vector<double> &s,int s_min,
    int s_max,double wr,double ws,int tail=56,double tail_mu=0.1,double r_init=-1e6,double r_end=-1e6,double s_init=-1e6,double s_end=-1e6);
/// APPLY VECTOR SHIFT TO A VECTOR DATA
vector<double> apply_shift(vector <double> &v,vector<double> &shift);
/// JOIN INDICATOR VALUES BY DATE
int joint_indicators_by_date(vector<string> &date0,vector<double> &i0,vector<string> &date1,vector<double> &i1,
                             vector<string> &date,vector<double> &f,vector<double> &g);
/// BASIC LINEAR INTERPOLATION. IT RETURNS THE INTERPOLATION VALUE (USING NEUMAN BOUNDARY CONDITION)
double linear_interpolation2(vector<double> &N,double t);
/// FORECASTING THE INDICATOR f(t) FROM THE INDICATOR g(t)
void indicators_forecast(vector<string> &dates,vector<double> &f,vector<double> &g,vector<double> &r0,
                         vector<double> &s0,double r_forecast_end,double s_forecast_end,int NforecastDays,int NdaysHermite,int NdaysHermiteDerivada);
/// MEAN OF THE LAST ELEMENTS OF A VECTOR
double mean_tail(vector<double> &i,int Ndays);
/// GEOMETRIC MEAN OF THE LAST ELEMENTS OF A VECTOR
double geometric_mean_tail(vector<double> &i,int Ndays);
