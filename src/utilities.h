/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file utilities.h
 * @brief auxiliary functions for the EpiInvert method to compute the reproduction number and a restored incidence curve
 *
 * @author Luis Alvarez <lalvarez@ulpgc.es>
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

#define R_COMPILE

using namespace std;

///  MACROS TO FIT THE SIMILARITY MEASURE BETWEEN INCIDENCE CURVE FOR FORECASTING
#define CORRELATION 1
#define L_ONE 0
/// FUNCTION TO READ FROM DISK A DATAFRAME IN A STRING FORM
vector< vector<string> > read_dataframe(char filename[],char sep);
/// HERMITE EXTRAPOLATION
void HermiteInterpolation(vector<double> &r,double r1,int NNewDays,int NNewDaysDerivative);
/// GAUSS CONVOLUTION
void gauss_conv(vector<double> &v,double sigma,int BoundaryCondition=1);
/// WEIGHTED MEDIAN OF A VECTOR x WITH WEIGHTS W
double weightedMedian(vector<double> x,vector<double> W);
/// LOG-NORMAL EVALUATION
double log_normal(double x,double mu,double sigma,double shift);
/// FUNCTION TO COMPUTE THE LINEAR REGRESSION (y=mx+n) BETWEEN VECTORS x AND y
/// IT RETURNS THE CORRELATION FACTOR
double linear_regression(vector<double> &x,vector<double> &y,double &m,double &n);
/// ERROR MANAGEMENT
void fprintf_demo_failure(char mes[300]);
/// BASIC STATISTICS: MEAN AND STANDARD DEVIATION (THE INPUT IS A SINGLE ARRAY)
void basic_statistics(vector<double> &i,double &mean,double &sd);
/// BASIC STATISTICS: MEAN AND STANDARD DEVIATION (THE INPUT A MATRIX)
void basic_statistics(vector < vector<double> > &i,double &mean,double &sd);
/// PERCENTIL COMPUTATION
double percentil(const int k,const vector<double> &x);
/// FUNCTION TO READ A DATABASE OF COVID-19 SEQUENCES
vector< vector<double> > read_matrix(const char name[]);
/// ESTIMATION OF A LINEAR REGRESSION INTERPOLATING IN THE LAST 7 DAYS
vector<double> last_week_regression_interpolation(const vector<double> &d);
/// EVALUATION OF THE INTERPOLATION POLYNOMIAL CENTERED IN THE LAST WEEK
double last_week_polynomial_evaluation(const int pos,const vector<double> &d, const vector<double> &P);
/// STANDARD GAUSS METHOD FOR SOLVING A LINEAR SYSTEM Au=b
vector< double > linear_system_solution(vector< vector<long double> > &A,vector< long double > &b);
/// FUNCTION TO READ COUNTRY DATA FROM THE FILE "owid-covid-data.csv". RETURNS  A VECTOR WITH THE NEW CASES.
vector<double>  read_country(const char C[],char date[],vector<double> &death);
/// FUNCTION TO READ MULTIPLE DATA FROM THE FILE name[]
vector< vector<double> > read_data_multiple(const char name[],time_t &current_time,bool death=false);
/// READING THE SERIAL INTERVAL FROM A FILE. IT RETURNS f0 (THE POSITION O ZERO IN THE SERIAL INTERVAL)
int read_si_distr(const char name[],vector<double> &si_distr);
/// BUILD A PARAMETRIC SERIAL INTERVAL FROM A SHIFTED LOG-NORMAL
int parametric_si_distr(double mean,double sd,double shift,vector<double> &si_distr);
/// BASIC LINEAR INTERPOLATION
double linear_interpolation(vector<double> &N,double t);
/// L2 ERROR (USED TO COMPUTE THE BEST MATCHING BETWEEN EPIESTIM AND EPIINVERT)
double L2(vector<double> &c1,double t1,vector<double> &c2,double t2,int kmin,int kmax);
/// FUNCTION TO READ THE DATA
void read_data_single(char name[],vector<double> &c);
/// FUNCTION TO COMPUTE THE PERCENTIL IN THE LAST DAYS
vector<double> back_percentil(vector<double> &i,int radius);
/// FUNCTION TO CONVERT A CHAR ARRAY WITH THE FORMAT "YYYY-MM-DD" TO time_t
time_t string2date(const char *date);
/// FUNCTION TO CONVERT time_t TO STRING FORMAT "YYYY-MM-DD"
string date2string(time_t time);
/// FUNCTION TO READ DATA CONTAINING THE DATE AND THE INCIDENCE VALUE
vector<double>  read_data_single_and_date(const char filename[],time_t &current_time);
/// FUNCTION TO READ THE FILE WITH THE FESTIVE DAYS
void read_festive_days(char name[],vector<string> &festive_dates);
/// FUNCTION TO INITIALIZE THE VECTOR TO DETERMINE THE FESTIVE DAYS
vector<int> daily_festive_day_initialization(time_t current_time, int i_size,vector<string> &festive_dates);
/// FUNCTION TO GET ONE OF THE STORED FESTIVE DAY SEQUENCES
vector<string> get_stored_festive_days(vector<double> &i);

double evaluation_init_extrapolation_14(int t,vector <double> &C);

double linear_regression_14(vector<double> &i, vector <double> &C);

double exponential_approximation_14(vector<double> &i,vector <double> &C);

