#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


#include "EpiInvertCore_q_variable.h"
///

// [[Rcpp::export]]
List EpiInvertC(
    NumericVector i_original0,
    String last_incidence_date,
    CharacterVector festive_days,
    NumericVector si_distr0,
    int shift_si_distr=0,
    int max_time_interval=50,
    double mean_si=12.267893,
    double sd_si=5.667547,
    double shift_si=-5.,
    double Rt_regularization_weight=5.,
    double seasonality_regularization_weight=5.
){
  clock_t t=clock();
  
  vector<double> i_original(i_original0.size());
  for(int k=0;k<i_original.size();k++) i_original[k]=i_original0[k];
  
  vector<double> si_distr;
  if(si_distr0.size()>0){
    si_distr=vector<double>(si_distr0.size());
    for(int k=0;k<si_distr.size();k++){
      si_distr[k]=si_distr0[k];
      //printf("si_distr[%d]=%lf\n",k,si_distr[k]);
    }
  }
  
  /// INPUT VARIABLES
  string last_incidence_dateC=string(last_incidence_date.get_cstring());/** DATE OF THE LAST DATA IN THE FORMAT YYYY-MM-DD */;
  vector <string> festive_daysC /** VECTOR OF FESTIVE OR ANOMALOUS DAYS IN THE FORMAT YYYY-MM-DD*/;
  for(int k=0;k<festive_days.size();k++) {
    if(strlen(festive_days[k])==10 && (festive_days[k][0]=='2' || festive_days[k][0]=='1'))
      festive_daysC.push_back(string(festive_days[k]));
  }
  
  //vector<double> si_distr;
  
  //int shift_si_distr /** SHIFT OF THE NON-PARAMETRIC SERIAL INTERVAL */;
  /// OUTPUTS
  vector<double> i_festive /** FESTIVE BIAS CORRECTED INCIDENCE CURVE*/;
  vector<double> i_bias_free /** BIAS FREE INCIDENCE CURVE*/;
  vector<double> i_restored /** FILTERED INCIDENCE CURVE WITH WEEKLY BIAS CORRECTION AND APPLICATION OF THE RENEWAL EQUATION*/;
  vector<double> Rt /** EpiInvert ESTIMATION OF THE REPRODUCTION NUMBER*/;
  vector<double> seasonality /** VECTOR WITH THE 7-DAY WEEKLY CORRECTION FACTORS (THE LAST VALUE OF seasonality CORRESPONDS TO THE LAST VALUE OF i0)*/;
  vector<double> Rt_CI95 /** 95% CONFIDENCE INTERVAL RADIUS FOR THE Rt ESTIMATE */;
  vector<string> dates /** DATE ASSOCIATED TO EACH INCIDENCE DATUM */;
  vector<bool> festive /** BOOL TO MANAGE WHETHER EACH INCIDENCE DATUM IS FESTIVE */;
  double RMSE_factor /** REDUCTION FACTOR OF THE RMSE BETWEEN THE INCIDENCE CURVE AND THE RENEWAL EQUATION BEFORE AND AFTER THE WEEKLY BIAS CORRECTION */;
  int iter_alternate_optimization /** NUMBER OF ITERATIONS OF THE ALTERNATE ALGORITHM TO COMPUTE Rt AND seasonality */;
  double power_a /** ESTIMATED POWER  IN THE RELATION i_bias_free[k] = i_restored[k] + epsilon[k]*i_restored[k]^a */;
  vector<double> epsilon /** ERROR DISTRIBUTION GIVEN BY  (i_bias_free[k] - i_restored[k])/i_restored[k]^a */;
  /// INPUT PARAMETERS
  int NweeksToKeepIncidenceSum=2 /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/;
  
  /// CHECK INPUT VALUES
  printf("Last Incidence Values\n");
  for(int k=i_original.size()-6;k<i_original.size();k++) printf("i[%d]=%1.0lf, ",k,i_original[k]);
  printf("\n");
  if(strlen(last_incidence_dateC.c_str())!=10 || last_incidence_dateC.c_str()[4]!='-' || last_incidence_dateC.c_str()[7]!='-'){
    stop("EpiInvert second argument must be a date in the format YYYY-MM-DD");
    return List::create();
  }
  printf("last_incidence_date %s\n",last_incidence_dateC.c_str());
  printf("Last festive days\n");
  for(int k=festive_daysC.size()-6;k<(int) festive_daysC.size();k++){
    if(k>=0) { printf("f[%d]=%s, ",k,festive_daysC[k].c_str()); }
  }
  printf("\n");
  printf("Last si_distr\n");
  for(int k=si_distr.size()-6;k<(int) si_distr.size();k++){
    if(k>=0) { printf("si_distr[%d]=%lf, ",k,si_distr[k]); }
  }
  
  
  printf("shift_si_distr=%d\n",shift_si_distr);
  printf("max_time_interval=%d\n",max_time_interval);
  printf("mean_si=%lf\n",mean_si);
  printf("sd_si=%lf\n",sd_si);
  printf("shift_si=%lf\n",shift_si);
  printf("Rt_regularization_weight=%lf\n",Rt_regularization_weight);
  printf("seasonality_regularization_weight=%lf\n",seasonality_regularization_weight);
  
  
  
  
  
  EpiInvertEstimate(
    i_original /** ORIGINAL INCIDENCE CURVE*/,
    last_incidence_dateC /** DATE OF THE LAST DATA IN THE FORMAT YYYY-MM-DD */,
    festive_daysC /** VECTOR OF FESTIVE OR ANOMALOUS DAYS IN THE FORMAT YYYY-MM-DD*/,
    si_distr /** NON-PARAMETRIC SERIAL INTERVAL (IF THIS VECTOR IS NOT EMPTY WE USE IT AS SERIAL INTERVAL)*/,
    shift_si_distr /** SHIFT OF THE NON-PARAMETRIC SERIAL INTERVAL */,
    /// OUTPUTS
    i_festive /** FESTIVE BIAS CORRECTED INCIDENCE CURVE*/,
    i_bias_free /** WEEKLY AND FESTIVE BIAS CORRECTED INCIDENCE CURVE*/,
    i_restored /** FILTERED INCIDENCE CURVE WITH WEEKLY BIAS CORRECTION AND APPLICATION OF THE RENEWAL EQUATION*/,
    Rt /** EpiInvert ESTIMATION OF THE REPRODUCTION NUMBER*/,
    seasonality /** VECTOR WITH THE 7-DAY WEEKLY CORRECTION FACTORS */,
    Rt_CI95 /** 95% CONFIDENCE INTERVAL RADIUS FOR THE Rt ESTIMATE */,
    dates /** DATE ASSOCIATED TO EACH INCIDENCE DATUM */,
    festive /** BOOL TO MANAGE WHETHER EACH INCIDENCE DATUM IS FESTIVE */,
    RMSE_factor /** REDUCTION FACTOR OF THE RMSE BETWEEN THE INCIDENCE CURVE AND THE RENEWAL EQUATION BEFORE AND AFTER THE WEEKLY BIAS CORRECTION */,
    iter_alternate_optimization /** NUMBER OF ITERATIONS OF THE ALTERNATE ALGORITHM TO COMPUTE Rt AND seasonality */,
    power_a /** ESTIMATED POWER  IN THE RELATION i_bias_free[k] = i_restored[k] + epsilon[k]*i_restored[k]^a */,
    epsilon /** ERROR DISTRIBUTION GIVEN BY  (i_bias_free[k] - i_restored[k])/i_restored[k]^a */,
    /// INPUT PARAMETERS
    mean_si /** MEAN OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    sd_si /** STANDARD DEVIATION OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    shift_si /** SHIFT OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    Rt_regularization_weight /** REGULARIZATION WEIGHT PARAMETER OF EpiInvert METHOD (DEFAULT VALUE: 5)*/,
    seasonality_regularization_weight /** WEIGHT PARAMETER OF THE REGULARIZATION  TERM FOR THE SEASONALITY q (DEFAULT VALUE 5) */,
    max_time_interval /** MAX SIZE OF THE INCIDENCE DATA USED TO COMPUTE Rt (DEFAULT VALUE: 9999). THIS PARAMETER IS USED TO REDUCE HE COMPUTATIONAL COST OF THE ALGORITHM WHEN WE ARE JUST INTERESTED IN THE LAST PART OF THE SEQUENCE */,
    NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/
  );
  
  
  printf("-> POWER OF THE RESTORED INCIDENCE IN THE INCIDENCE MODEL : %lf\n",power_a);
  
  /// WRITE OUTPUT
  if(true){
    FILE *g;
    g=fopen ("EpiInvertEstimateResults.csv", "w");
    if(g==NULL){
      printf("\nProblems opening Rt.csv. Maybe it is already opened ?\n\n");
      system("pause");
      exit(0);
    }
    
    fprintf(g,"date;festive;i_originalC;i_festive;i_bias_free;i_restored;Rt;Rt_CI95;seasonality;(i_bias_free-i_restored)/i_restored ^ %lf\n",power_a);
    
    for(int k=0;k<(int) i_original.size();k++){
      fprintf(g,"%s;",dates[k].c_str());
      if(festive[k]==false) fprintf(g,"FALSE;");
      else fprintf(g,"TRUE;");
      fprintf(g,"%lf; %lf; %lf; %lf; %lf; %lf; %lf; %lf\n",i_original[k],i_festive[k],i_bias_free[k],i_restored[k],Rt[k],Rt_CI95[k],seasonality[k],epsilon[k]);
    }
    fclose(g);
  }
  
  
  printf("-> END OF PROGRAM EXECUTION\n");
  t = clock() - t;
  printf("-> EXECUTION TIME : %f SECONDS\n",((float)t)/CLOCKS_PER_SEC);
  
  List results = List::create(
    Named("i_original") = i_original ,
    Named("i_festive") = i_festive ,
    Named("i_bias_free") = i_bias_free ,
    Named("i_restored")  = i_restored,
    Named("Rt")  = Rt,
    Named("seasonality")  = seasonality,
    Named("Rt_CI95")  = Rt_CI95,
    Named("dates")  = dates,
    Named("festive")  = festive,
    Named("epsilon")  = epsilon,
    Named("power_a")  = power_a,
    Named("si_distr") = si_distr,
    Named("shift_si_distr") = shift_si_distr
  );
  return results;
}

