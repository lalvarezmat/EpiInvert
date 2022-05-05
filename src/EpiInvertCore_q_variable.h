/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file EpiInvertCore.h
 * @brief EpiInvert method to compute the reproduction number and a restored incidence curve
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
#include <omp.h>
//#include <random>
#include "utilities.h"

#define CASE false
#define INST true


using namespace std;

/// ESTIMATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE.
vector <double> initial_incidence_growth_estimation(vector<double> &i);
/// EVALUATION OF THE RENEWAL EQUATION FORMULA
double RenewalEquation(const int t,const vector<double> &si_distr,const int k0,const vector<double> &R,const vector<double> &i,
                       vector<double> &P,const bool EpiEst);
/// EpiInvert ESTIMATION. IT RETURNS THE TIME WHERE Rt STARTS TO BE COMPUTED (A NEGATIVE VALUE IN CASE OF FAILURE)
int EpiInvert(
    vector<double> &i0,vector<double> &i1,const char incidence_filename[], const char si_distr_filename[],const double Rt_regularization_weight,const bool RenewalEquationModel,
    const bool WeeklyBiasCorrection,const int max_time_interval,
    vector<double> &i_renewal_equation,vector<double> &Rt,vector<double> &seasonality,vector<double> &V,double &RMSE_factor,int &iter_alternate_optimization,
    vector <string> &festive_days,time_t &current_day,int NweeksToKeepIncidenceSum,double seasonality_regularization_weight);
/// EpiInvert ESTIMATION. VERSION FOR THE R PACKAGE
void EpiInvertEstimate(
    /// INPUT DATA
    vector<double> &i_original /** ORIGINAL INCIDENCE CURVE*/,
    const string last_incidence_date /** DATE OF THE LAST DATA IN THE FORMAT YYYY-MM-DD */,
    vector <string> &festive_days /** VECTOR OF FESTIVE OR ANOMALOUS DAYS IN THE FORMAT YYYY-MM-DD*/,
    vector<double> &si_distr /** NON-PARAMETRIC SERIAL INTERVAL (IF THIS VECTOR IS NOT EMPTY WE USE IT AS SERIAL INTERVAL)*/,
    int &shift_si_np /** SHIFT OF THE NON-PARAMETRIC SERIAL INTERVAL */,
    /// OUTPUTS
    vector<double> &i_festive /** FESTIVE BIAS CORRECTED INCIDENCE CURVE*/,
    vector<double> &i_bias_free /** BIAS FREE INCIDENCE CURVE*/,
    vector<double> &i_restored /** FILTERED INCIDENCE CURVE WITH WEEKLY BIAS CORRECTION AND APPLICATION OF THE RENEWAL EQUATION*/,
    vector<double> &Rt /** EpiInvert ESTIMATION OF THE REPRODUCTION NUMBER*/,
    vector<double> &seasonality /** VECTOR WITH THE 7-DAY WEEKLY CORRECTION FACTORS (THE LAST VALUE OF seasonality CORRESPONDS TO THE LAST VALUE OF i0)*/,
    vector<double> &Rt_CI95 /** 95% CONFIDENCE INTERVAL RADIUS FOR THE Rt ESTIMATE */,
    vector<string> &dates /** DATE ASSOCIATED TO EACH INCIDENCE DATUM */,
    vector<bool> &festive /** BOOL TO MANAGE WHETHER EACH INCIDENCE DATUM IS FESTIVE */,
    double &RMSE_factor /** REDUCTION FACTOR OF THE RMSE BETWEEN THE INCIDENCE CURVE AND THE RENEWAL EQUATION BEFORE AND AFTER THE WEEKLY BIAS CORRECTION */,
    int &iter_alternate_optimization /** NUMBER OF ITERATIONS OF THE ALTERNATE ALGORITHM TO COMPUTE Rt AND seasonality */,
    double &a /** ESTIMATED POWER  IN THE RELATION i_bias_free[k] = i_restored[k] + epsilon[k]*i_restored[k]^a */,
    vector<double> &epsilon /** ERROR DISTRIBUTION GIVEN BY  (i_bias_free[k] - i_restored[k])/i_restored[k]^a */,
    /// INPUT PARAMETERS
    const double mean_si=12.267893 /** MEAN OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double sd_si=5.667547 /** STANDARD DEVIATION OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double shift_si=-5. /** SHIFT OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double Rt_regularization_weight=5. /** REGULARIZATION WEIGHT PARAMETER OF EpiInvert METHOD (DEFAULT VALUE: 5)*/,
    const double seasonality_regularization_weight=5. /** WEIGHT PARAMETER OF THE REGULARIZATION  TERM FOR THE SEASONALITY q (DEFAULT VALUE 5) */,
    const int max_time_interval=9999 /** MAX SIZE OF THE INCIDENCE DATA USED TO COMPUTE Rt (DEFAULT VALUE: 9999). THIS PARAMETER IS USED TO REDUCE HE COMPUTATIONAL COST OF THE ALGORITHM WHEN WE ARE JUST INTERESTED IN THE LAST PART OF THE SEQUENCE */,
    const int NweeksToKeepIncidenceSum=2 /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/
);
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES
void IncidenceExtrapolationByLearning(vector<double> &i,const vector< vector <double > > &i42,const vector< vector <double > > &i56,const int NweeksToKeepIncidenceSumeeksBackToForeCast,double sigma,int COMPARISON_TYPE,int index);
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES (USING THE MEDIAN IN THE LAST 5 WEEKS ESTIMATION)
void IncidenceExtrapolationByLearningMedian5Weeks(vector<double> &i,const vector< vector <double > > &i42,const vector< vector <double > > &i56,int COMPARISON_TYPE);
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES (USING THE MEDIAN IN THE LAST 3 WEEKS ESTIMATION)
void IncidenceExtrapolationByLearningMedian3Weeks(vector<double> &i,const vector< vector <double > > &i42,const vector< vector <double > > &i56,int COMPARISON_TYPE);
/// ESTIMATION OF THE EFFECTIVE REPRODUCTIVE NUMBER
void Rt_estimation(vector<double> &c,vector<double> &P,const vector<double> &si_distr,
                   const int k0,const double Rt_regularization_weight,const vector<double> &nf,const bool EpiEst,vector<double> &R,const vector<double> &i0,int NweeksToKeepIncidenceSum);
/// ESTIMATION OF THE WEEKLY BIAS CORRECTION FACTORS (IT RETURNS A VECTOR WITH THE FACTORS
vector<double> periodic_7days(const vector<double> &i,const vector<double> &R,const vector<double> &si_distr,const int k0,
                              vector<double> &P,const bool EpiEst,const int tmin,const int tmax);
/// Rt AND WEEKLY CORRECTION FACTORS q ESTIMATION RETURNS I (THE RMSE RATIO AFTER WEEKLY BIAS CORRECTION)
double Rt_q_estimation(vector<double> &i,vector<double> &P,const vector<double> &si_distr,
                       const int f0,const double Rt_regularization_weight,const vector<double> &nf,const bool EpiEst,vector<double> &R,vector<double> &q,
                       vector<double> &iRenEq,int &iter_alternate_optimization,const int MaxIter,bool WeeklyBiasCorrection,vector<int> &daily_festive_day,const vector<double> &i0,int NweeksToKeepIncidenceSum,double seasonality_regularization_weight);
/// DATA PRE-PROCESSING. DEALING WITH LACK OF DATA OR NEGATIVE INCIDENCE VALUES
void lack_of_data_processing (vector<double> &c);
/// COMPUTATION OF THE EPIESTIM ESTIMATION OF Rt
vector<double> EpiEstim(vector<double> &c0,double a,double b,int tau,char si_distr_filename[]);
/// DATA PREPROCESSING
vector<double> data_pre_processing(const vector< vector<double> >  &cV0,int NLastDaysToRemove,int max_time_interval);
/// WE WRITE THE RESULTS IN THE FILES SummaryExecutionOutcomes.txt and Rt.csv
void write_files(vector<double> &R,vector<double> &Epi,vector<double> &c,vector<double> &i0,vector<double> &i1,vector<double> &V,bool WeeklyBiasCorrection,
                 vector<double> &q,double error_min,double MinE,double tmin,int iter_alternate_optimization,bool RenewalEquationModel,
                 time_t current_day,vector<string> festive_days);

