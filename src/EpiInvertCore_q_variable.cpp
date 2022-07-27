/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file EpiInvertCore.cpp
 * @brief EpiInvert method to compute the reproduction number and a restored incidence curve
 *
 * @author Luis Alvarez <lalvarez@ulpgc.es>
 */

#define R_COMPILE

#include "EpiInvertCore_q_variable.h"


///----------------------------------------------------------------------------------------------------
/// DATA PREPROCESSING : THE INPUT IS A MATRIX OF INCIDENCE VALUES (1 ROW PER COUNTRY)
/// THE OUTPUT IS A SINGLE CUMMULATIVE INCIDENCE. THE PRE-PROCESSING INCLUDES THE FOLLOWING STEPS:
/// (1) WE OPTIONALLY REMOVE THE LAST VALUES OF THE INCIDENCE SEQUENCE "USED TO COMPUTE THE DAILY VARIABILITY OF Rt IN THE LAST DAYS"
/// (2) FOR EACH COUNTRY WE MANAGE NON-POSITIVE VALUES.
/// (3) WE ACCUMULATE THE DIFFERENT COUNTRY INCENDENCE IN A SINGLE VECTOR
/// (4) WE REMOVE THE FIRST SMALL VALUES UP TO OBTAIN THAT THE CUMMULATIVE INCIDENCE IS LARGER THAN 100
/// (5) WE ADJUST THE FINAL SIZE OF THE INCIDENCE CURVE.
vector<double> data_pre_processing(
    const vector< vector<double> >  &cV0, /// ORIGINAL DATA
    int NLastDaysToRemove, /// NUMBER OF DAYS TO REMOVE FROM THE END OF THE INCIDENCE SEQUENCE "USED TO COMPUTE THE DAILY VARIABILITY OF Rt IN THE LAST DAYS"
    int max_time_interval) /// MAX SIZE OF THE OUTPUT INCIDENCE CURVE ALLOWED
{
  /// WE COMPUTE THE INITIAL INCIDENCE SIZE REMOVING NLastDaysToRemove FROM THE END OF THE SEQUENCE
  int i_size = cV0[0].size()- NLastDaysToRemove;
  
  /// OUTPUT INCIDENCE CURVE
  vector<double> i(i_size,0.);
  
  /// PRE-PROCESSING OF COUNTRIES INCIDENCE AND COMPUTATION OF THE CUMMULATIVE INCIDENCE
  for(int k=0;k<(int) cV0.size();k++){
    vector<double> c2=cV0[k];
    if((int) cV0[k].size()>i_size) c2.resize(i_size);
    /// WE MANAGE NON-POSITIVE VALUES.
    lack_of_data_processing (c2);
    for(int n=0;n<(int) c2.size();n++) i[n]+=c2[n];
  }
  /// WE REMOVE NON-POSITIVE VALUES FROM THE END OF THE INCIDENCE SEQUENCE
  while(i.size()>1 && i[i.size()-1]<=0) i.resize(i.size()-1);
  
  
  /// REMOVING SMALL VALUES FROM THE BEGINNING OF THE SEQUENCE
  int t0;
  double sum=0;
  for(t0=0;t0<(int) i.size()-10;t0++){
    sum+=i[t0];
    if(sum>=100) break;
  }
  if(t0>0) i.erase(i.begin(),i.begin()+t0);
  
  if(i.size()<20){
#ifndef R_COMPILE
    //printf("The number of samples is too small : %d samples\n",(int) i.size());
    char mes[300];
    sprintf(mes,"The number of samples is too small : %d samples\n",(int) i.size());
    //printf("%s\n",mes);
    fprintf_demo_failure(mes);
#endif
    return vector<double>();
  }
  
  /// WE REDUCE (IF NECESSARY) THE SIZE OF THE INCIDENCE CURVE
  if((int) i.size()>max_time_interval){
    int NdaysToErase=i.size()-max_time_interval;
    i.erase(i.begin(),i.begin()+NdaysToErase);
  }
  
  return i;
}

///----------------------------------------------------------------------------------------------------
/// DATA PRE-PROCESSING. DEALING WITH LACK OF DATA OR NEGATIVE INCIDENCE VALUES
void lack_of_data_processing (vector<double> &c){
  
  /// WE REMOVE FROM THE SEQUENCE THE LAST NON-POSITIVE VALUES.
  while( c[c.size()-1]<=0. && c.size()>0) c.resize(c.size()-1);
  
  /// IF AN INCIDENCE VALUE IS POSITIVE AND THE PREVIOUS ONES ARE ZERO
  /// WE ASSIGN TO ALL AFECTED VALUES A CONSTANT VALUE PRESERVING THE TOTAL NUMBER OF CASES.
  /// THE MAXIMUM NUMBER OF ALLOWED CONSECUTIVE ZERO VALUES ARE 7.
  for(int k=(int) c.size()-1; k>0;k--){
    if(c[k]<0 || !(c[k-1]==0) ) continue;
    int m=1;
    while(m<7){
      if(k-m-1<0 || !(c[k-m-1]==0) ) break;
      m++;
    }
    if(m==7) break;
    double average=c[k]/(m+1);
    for(int n=k;n>=k-m;n--) c[n]=average;
    k=k-m+1;
  }
  
  /// NEGATIVE INCIDENCE VALUES ARE NOT ALLOWED. WE REPLACE THEN BY 0.
  for(int k=0; k<(int) c.size();k++){
    if(c[k]<0) c[k]=0;
  }
  
  return;
  
}

/// ESTIMATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE. USED TO EXTRAPOLATE THE INCIDENCE BEFORE THE INITIAL TIME
/// IT TAKES THE BEST APPROXIMATION BETWEEN A LINEAR REGRESION APPROXIMATION (i[t]=C[0]*t+C[1]) AND AN EXPONENTIAL APPROXIMATION (i[t]=C[0]*e^(C[1]*t)+C[2])
/// IN THE FIRST 14 DAYS. IT RETURNS THE VECTOR C[]. WE RECOGNIZE THE TYPE OF INTERPOLATION ACCORDINGLY TO THE SIZE OF VECTOR C[]
vector <double> initial_incidence_growth_estimation(
    vector<double> &i /** daily infected */)
{
  
  /// LINEAR REGRESION IN THE FIRST 14 DAYS
  vector<double> Cl; /// vector with the parameters of the linear regression (i[t]=Cl[0]*t+Cl[1])
  double error1=linear_regression_14(i,Cl);
  
  /// EXPONENTIAL APPROXIMATION IN THE FIRST 14 DAYS
  vector<double> Ce; /// vector with the parameters of the exponential approximation (i[t]=Ce[0]*e^(Ce[1]*t)+Ce[2])
  double error2=exponential_approximation_14(i,Ce);
  
  /// WE LOOK FOR THE BEST APPROXIMATION TTT
  if(error1<error2 && evaluation_init_extrapolation_14(0,Cl)>0 && evaluation_init_extrapolation_14(-5,Cl)>0){
    //printf("linear regression has a lower error\n");
    return Cl;
  }
  else{
    //printf("exponential approximation has a lower error\n");
    return Ce;
  }
  
}


/// -----------------------------------------------------------------------------------------
/// EVALUATION OF THE RENEWAL EQUATION FORMULA
double RenewalEquation(
    const int t /** TIME WHERE THE FORMULA IS EVALUATED */,
    const vector<double> &si_distr /** SERIAL INTERVAL */,
    const int k0 /** POSITION OF 0 VALUE IN THE SERIAL INTERVAL  */,
    const vector<double> &R /** Rt */,
    const vector<double> &i /** DAILY TESTED POSITIVE */,
    vector<double> &Pi /** PARAMETERS OF THE EXTRAPOLATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE*/,
    const bool RenewalEquationModel /**THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */)
{
  double RenewalEquationEvaluation=0.;
  
  /// EXTRAPOLATION TOWARDS THE FUTURE OF THE INCIDENCE CURVE USING THE LAST SEVEN DAYS
  vector<double> P;
  if(i.size()>7) P=last_week_regression_interpolation(i);
  
  /// EVALUATION OF THE RENEWAL EQUATION FORMULA
  for(int k=0;k<(int) si_distr.size();k++){
    int ts=t-k+k0; /// value t-s in the formula
    double Rts,its; /// values R(t-s) and i(t-s) in the formula
    
    /// IF ts<0 WE USE THE EXTRAPOLATION OF THE INITIAL GROWTH TO EVALUATE i(t-s)
    if(ts<0){
      its=evaluation_init_extrapolation_14(ts,Pi);
    }
    /// IF ts IS AFTER THE CURRENT TIME, WE USE THE EXTRAPOLATION TOWARDS THE FUTURE TO EVALUATE i(t-s)
    else if(ts>=(int) i.size()){
      if(i.size()>7) its=last_week_polynomial_evaluation(ts,i,P);
      else{
        its=evaluation_init_extrapolation_14(ts,Pi);
      }
    }
    /// OTHERWISE i(t-s) IS COMPUTED DIRECTLY
    else its=i[ts];
    
    /// COMPUTATION OF Rt(t-s)
    /// IF (t-s)<0 Rt(t-s)=Rt(0)
    if(ts<0) Rts=R[0];
    /// IF (t-s)> current time Rt(t-s)=Rt(current time)
    else if(ts>=(int) R.size()) Rts=R[R.size()-1];
    /// OTHERWISE Rt(t-s) IS COMPUTED DIRECTLY
    else{
      Rts=R[ts];
    }
    
    /// WE MANAGE THE CASE OF THE INSTANTANEOUS REPRODUCTIVE NUMBER
    if(RenewalEquationModel==INST){
      if(t<0)  Rts=R[0];
      else if(t>=(int) R.size()) Rts=R[R.size()-1];
      else Rts=R[t];
    }
    
    /// WE ACCUMULATE THE VALUE TO THE RENEWAL EQUATION FORMULA
    RenewalEquationEvaluation+=si_distr[k]*its*Rts;
  }
  
  return RenewalEquationEvaluation;
}


/// -----------------------------------------------------------------------------------------
/// COMPUTATION OF Rt (OUTPUT) USING A LINEAR SYSTEM
void  LinearSystemRt( //III
    const vector<double> &i /** DAILY TESTED POSITIVE */,
    vector<double> &Pi /** PARAMETERS OF THE EXTRAPOLATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE*/,
    const vector<double> &si_distr /** SERIAL INTERVAL */,
    const int k0 /** LOCATION ON 0 IN THE SERIAL INTERVAL */,
    const vector<double> w /** WEIGHTS IN THE REGULARIZATION TERM OF THE SERIAL INTERVAL */,
    const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
    vector<double> &R /** OUTPUT Rt ESTIMATION*/,
    const vector<double> &nf /** MEDIAN OF THE INCIDENCE IN THE LAST 21 DAYS TO NORMALIZE THE ENERGY */,
    const vector<double> &i0 /** DAILY TESTED POSITIVE */,
    int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/)
{
  /// INDEPENDENT TERM OF THE LINEAR SYSTEM
  vector<long double> b;
  for(int k=0; k<(int) i.size();k++){
    b.push_back(i[k]);
  }
  
  /// LINEAR REGRESSION ESTIMATION OF THE INCIDENCE CURVE FOR THE LAST 7 DAYS
  vector<double> P=last_week_regression_interpolation(i);
  
  /// COMPUTATION OF THE MATRIX OF THE LINEAR SYSTEM
  vector < vector<long double> > A(b.size(),vector<long double>(b.size(),0.));
  for(int n=0;n<(int) A.size();n++){
    for(int k=0;k<(int) si_distr.size();k++){
      int t=n-k+k0;
      if(t<0){
        A[n][0]+=si_distr[k]*evaluation_init_extrapolation_14(t,Pi);
        continue;
      }
      if(t>=(int) i.size()){
        double i2=last_week_polynomial_evaluation(t,i,P);
        A[n][A[n].size()-1]+=si_distr[k]*i2;
        continue;
      }
      A[n][t]+=si_distr[k]*i[t];
    }
  }
  
  /// INSTANTANEOUS REPRODUCTION NUMBER (CORI ET AL.)
  if(RenewalEquationModel==INST){
    for(int n=0;n<(int) A.size();n++){
      for(int k=0;k<(int) A.size();k++){
        if(k==n) continue;
        A[n][n]+=A[n][k];
        A[n][k]=0;
      }
    }
  }
  
  /// NONLINEAR NORMALIZATION  EEE
  for(int k=0;k<(int) b.size();k++){
    b[k]/=nf[k];
    for(int m=0;m<(int) b.size();m++) A[k][m]/=nf[k];
  }
  
  /// SOLUTION COMPUTATION
  vector < vector<long double> > B(A[0].size(),vector<long double>(A[0].size(),0.));
  vector<long double> b2(A[0].size(),0.);
  
  for(int i=0;i<(int) B.size();i++){
    for(int k=0;k<(int) A.size();k++){
      b2[i]+=A[k][i]*b[k];
      for(int j=0;j<(int) B.size();j++){
        B[i][j]+=A[k][i]*A[k][j];
      }
    }
  }
  
  /// WE ADD THE REGULARITY TERM TO THE SYSTEM
  B[0][0]+=w[0];  B[B.size()-1][B.size()-1]+=w[B.size()-1];
  B[0][1]-=w[0];  B[B.size()-1][B.size()-2]-=w[B.size()-1];
  
  for(int k=1;k<(int) B.size()-1;k++){
    B[k][k-1]-=w[k];
    B[k][k]+=w[k]+w[k+1];
    B[k][k+1]-=w[k+1];
  }
  
  /// WE ADD THE LAGRANGE MULTIPLIERS
  int f0=-k0,f=si_distr.size()-k0;
  
  /// NUMBER OF PERIODS TO APPLY THE AVERAGE CONSERVATION OF THE INCIDENCE
  int N=i.size()/(7*NweeksToKeepIncidenceSum)-2;
  if(N<0) N=0;
  
  /// WE TAKE MEMORY FOR THE SYSTEM MATRIX TAKING INTO ACCOUNT THE LAGRANGE MULTIPLIERS
  vector < vector<long double> > C(b.size()+1+N,vector<long double>(b.size()+1+N,0.));
  
  /// WE COPY IN THE NEW MATRIX THE ONE BEFORE USING THE LAGRANGE MULTIPLIERS
  for(int n=0;n<(int) i.size();n++)
    for(int m=0;m<(int) i.size();m++)
      C[n][m]=B[n][m];
  
  /// WE ADD THE LAGRANGE MULTIPLIERS CONSTRAINTS TO THE SYSTEM
  for(int m=0;m<=N;m++){
    int t0=0,t1=i.size();
    if(m>0){
      t1=i.size()-(m-1)*7*NweeksToKeepIncidenceSum;
      t0=t1-7*NweeksToKeepIncidenceSum;
    }
    long double sum=0.;
    for(int k=t0;k<t1;k++) sum+=i[k];
    b2.push_back(sum);
    
    for(int s=t0-f;s<t1-f0;s++){
      sum=0;
      for(int t=t0;t<t1;t++){
        int t2=t-s+k0;
        if(t2<0 || t2>=(int) si_distr.size()) continue;
        sum+=si_distr[t2];
      }
      
      if(s<0){
        C[i.size()+m][0] += sum*evaluation_init_extrapolation_14(s,Pi);
      }
      else if(s>=(int) i.size()){
        C[i.size()+m][i.size()-1] += sum*last_week_polynomial_evaluation(s,i,P);
      }
      else{
        C[i.size()+m][s] += sum*i[s];
      }
    }
    
    for(int k=0;k<(int) i.size();k++) C[k][i.size()+m]=C[i.size()+m][k];
  }
  
  /// WE SOLVE THE LINEAR SYSTEM
  R=linear_system_solution(C,b2);
  
  /// WE REMOVE FROM THE SOLUTION THE LAGRANGE MULTIPLIERS.
  R.resize(i.size());
}

/// ----------------------------------------------------------------------------------------------------------
/// COMPUTATION OF it USING A LINEAR SYSTEM TAKING INTO ACCOUNT THE FESTIVE DAYS. IT MODIFIES THE VALUES OF VECTOR i
void  LinearSystem_it( //III
    vector<double> &i0 /** DAILY TESTED POSITIVE ORIGINAL */,
    vector<double> &i /** DAILY TESTED POSITIVE WEEKLY BIAS CORRECTED*/,
    vector<double> q /**  7-DAY WEEKLY BIAS CORRECTION FACTORS */,
    const vector<double> &si_distr /** SERIAL INTERVAL */,
    const int k0 /** LOCATION ON 0 IN THE SERIAL INTERVAL */,
    const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
    const vector<double> &R /** REPRODUCTION NUMBER */,
    const vector<double> &nf /** 21-DAY BACKWARD MEDIAN OF THE INCIDENCE FOR ENERGY NORMALIZATION */,
    const vector<int> daily_festive_day /** CONTROL OF PUBLIC DAILY FESTIVE DAYS (==1 -> FESTIVE DAY, ==0  -> OTHERWISE*/)
{
  
  /// INDEPENDENT TERM OF THE LINEAR SYSTEM
  vector<long double> b;
  for(int k=0; k<(int) i.size();k++){
    if(daily_festive_day[k]==0 || k<=(int) si_distr.size()-k0) b.push_back(i[k]);
    else b.push_back(0.);
  }
  
  /// LINEAR REGRESSION ESTIMATION OF THE INCIDENCE CURVE FOR THE LAST 7 DAYS
  vector<double> P=last_week_regression_interpolation(i);
  
  /// COMPUTATION OF THE MATRIX OF THE LINEAR SYSTEM
  vector < vector<long double> > A(b.size(),vector<long double>(b.size(),0.));
  if(RenewalEquationModel==CASE){
    for(int n=0;n<(int) A.size();n++){
      A[n][n]=1;
      if( daily_festive_day[n]==0  || n<=(int) si_distr.size()-k0){ continue; }
      for(int k=0;k<(int) si_distr.size();k++){
        int t=n-k+k0;
        int m=n-k+k0;
        if(t>=(int) i.size()){
          double i2=last_week_polynomial_evaluation(t,i,P);
          b[n]+=si_distr[k]*i2*R[R.size()-1];
          continue;
        }
        A[n][m]-=si_distr[k]*R[t];
      }
    }
  }
  else{
    for(int n=0;n<(int) A.size();n++){
      A[n][n]=1;
      if( daily_festive_day[n]==0  || n<=(int) si_distr.size()-k0){ continue; }
      for(int k=0;k<(int) si_distr.size();k++){
        int t=n-k+k0;
        int m=n-k+k0;
        if(t>=(int) i.size()){
          double i2=last_week_polynomial_evaluation(t,i,P);
          b[n]+=si_distr[k]*i2*R[R.size()-1];
          continue;
        }
        A[n][m]-=si_distr[k]*R[n+k0];
      }
    }
  }
  
  
  /// NONLINEAR NORMALIZATION
  for(int k=0;k<(int) b.size();k++){
    b[k]/=nf[k];
    for(int m=0;m<(int) b.size();m++) A[k][m]/=nf[k];
  }
  
  /// WE ADD TO THE SYSTEM THE FESTIVE DAY ADDITION CONSTRAINT IN THE ENERGY
  double sum=0;
  for(int k=0; k<(int) i.size();k++){
    /// WE CHECK IF k IS A FESTIVE DAY
    if(daily_festive_day[k]==1){
      /// WE COMPUTE THE NUMBER OF AFFECTED DAYS
      int m=1,Ndays=0;
      while(k+m< (int) i.size() && daily_festive_day[k+m++]==1) Ndays++;
      /// WE COMPUTE THE WEIGHT IN THE ENERGY FOR THE FESTIVE DAY CONSTRAINT
      double beta=pow(2.,(double) i.size()-k-4);
      if(beta>1e6) beta=1e6;
      int N=i.size()-k-1;
      if(N>Ndays) N=Ndays;
      /// WE ADD TO THE SYSTEM THE FESTIVE DAY CONSTRAINT
      if(k+1<(int) i.size()){
        A.push_back(vector<long double>(b.size(),0.));
        sum=0.;
        for(int m=0;m<=N;m++){
          A[A.size()-1][k+m]=beta/q[k+m]/nf[k];
          sum+=i0[k+m];
        }
        b.push_back(beta*sum/nf[k]);
        k+=N+1;
      }
    }
  }
  
  /// SOLUTION COMPUTATION
  vector < vector<long double> > B(A[0].size(),vector<long double>(A[0].size(),0.));
  vector<long double> b2(A[0].size(),0.);
  
  for(int m=0;m<(int) B.size();m++){
    for(int k=0;k<(int) A.size();k++){
      b2[m]+=A[k][m]*b[k];
      for(int j=0;j<(int) B.size();j++){
        B[m][j]+=A[k][m]*A[k][j];
      }
    }
  }
  
  vector<double> u=linear_system_solution(B,b2);
  
  for(int k=0;k<(int) u.size();k++) i[k]=u[k];
  
}


///----------------------------------------------------------------------------------------------------
/// ESTIMATION OF THE EFFECTIVE REPRODUCTION NUMBER
void Rt_estimation( //III
    vector<double> &c /** daily infected */,
    vector<double> &Pi /** PARAMETERS OF THE EXTRAPOLATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE*/,
    const vector<double> &si_distr /** serial interval distribution */,
    const int k0 /** position of the value 0 in the serial interval */,
    const double Rt_regularization_weight /** weigth for the energy minimization*/,
    const vector<double> &nf /** median used to normalize the incidence curve */,
    const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
    vector<double> &R /** output computation of Rt*/,
    const vector<double> &i0 /** ORIGINAL DAILY INFECTED */,
    int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/
)
{
  vector<double> d=c;
  
  /// REGULARIZATION WEIGHT OF ENERGY ERROR
  vector<double> w(c.size()+1,Rt_regularization_weight);
  
  /// WE INCREASE THE REGULARIZATION WEIGHT IF THE VALUE OF Rt IS TOO SMALL OR NEGATIVE.
  for(int k=0;k< (int) R.size();k++){
    if(R[k]<0.1) w[k]*=10;
    if(R[k]<0.)  w[k]*=10;
  }
  
  /// COMPUTATION OF THE APPROXIMATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE FOR EXTRAPOLATION PURPOSES
  if(Pi.size()==0) Pi=initial_incidence_growth_estimation(c);
  
  /// COMPUTATION OF Rt
  LinearSystemRt(d,Pi,si_distr,k0,w,RenewalEquationModel,R,nf,i0,NweeksToKeepIncidenceSum);  //III
  
  /// RECOMPUTE Rt IN THE CASE Rt IS NEGATIVE
  for(int k=0;k<5;k++){
    bool recompute=false;
    for(int k=0;k< (int) R.size();k++){
      if(R[k]<0.1) w[k]*=10;
      if(R[k]<0.){
        w[k]*=10;
        recompute=true;
      }
    }
    if(recompute==false) break;
    LinearSystemRt(d,Pi,si_distr,k0,w,RenewalEquationModel,R,nf,i0,NweeksToKeepIncidenceSum);
  }
  
  
  if(R.size()==0){
#ifndef R_COMPILE
    char mes[300];
    sprintf(mes,"Problems computing Rt\n");
    //printf("%s\n",mes);
    fprintf_demo_failure(mes);
#endif
    return;
  }
}

///----------------------------------------------------------------------------------------------------
/// ESTIMATION OF THE WEEKLY BIAS CORRECTION FACTORS (IT RETURNS A VECTOR WITH THE FACTORS
/// WE USE LAGRANGE MULTIPLIERS TO INCLUDE THE INCIDENCE MEAN PRESERVING CONSTRAINT
vector<double> periodic_7days(
    const vector<double> &i /** incidence curve */,
    const vector<double> &R /** Rt estimate with EpiInvert*/,
    const vector<double> &si_distr /** serial interval distribution */,
    const int k0 /** position of the value 0 in the serial interval */,
    vector<double> &Pi /** PARAMETERS OF THE EXTRAPOLATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE*/,
    const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
    const int tmin /** min date to compare the incidence curve (typically i.size()-1-NweeksToKeepIncidenceSumeeks*7)*/,
    const int tmax /** max date to compare the incidence curve (typically i.size()-1)*/,
    const vector<double> &nf /** median used to normalize the incidence curve in the energy*/,
    const vector<double> &i0 /** original incidence curve*/,
    int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/,
    double seasonality_regularization_weight /** weight of the regularization term for the seasonality q */)
{
  
  vector < vector<double> > A(i.size(),vector<double>(i.size(),0.)); /// matrix to store equations
  vector<long double> b(i.size()+1,0.); /// independent vector to the final linear system including Lagrange multipliers
  
  vector<double> P=last_week_regression_interpolation(i);
  /// WE FILL THE MATRIX TO STORE EQUATIONS AND b
  for(int t=0;t<(int) i.size();t++){
    A[t][t]=i[t];
    /// CONTRIBUTION OF RENEWAL EQUATION
    for(int s=0;s<(int) si_distr.size();s++){
      //int m=(n-s+k0+700)%7;
      int t1=t-s+k0;
      double R0,i0;
      if(t1<0) i0=evaluation_init_extrapolation_14(t1,Pi);
      else if(t1>=(int) i.size()){
        i0=last_week_polynomial_evaluation(t1,i,P);
      }
      else i0=i[t1];
      if(t1<0){
        t1=(t1+700)%7;
        R0=R[0];
      }
      else if(t1>=(int) R.size()){
        R0=R[R.size()-1];
        t1=i.size()-7-(t1-R.size())%7;
      }
      else{
        R0=R[t1];
      }
      if(RenewalEquationModel==INST){
        if(t>=(int) R.size()) R0=R[R.size()-1];
        else R0=R[t];
      }
      A[t][t1]-=si_distr[s]*i0*R0;///III
    }
    b[i.size()]+=i0[t];
  }
  
  for(int t=0;t<(int) i.size();t++){
    for(int j=0;j<(int) i.size();j++){
      A[t][j]/=nf[t];
    }
  }
  
  /// CONSTRUCTION OF THE LINEAR SYSTEM (INCLUDING LAGRANGE MULTIPLIER)
  int N=i.size()/(7*NweeksToKeepIncidenceSum)-2;
  if(N<0) N=0;
  //N=0;
  vector < vector<long double> > B(i.size()+1+N,vector<long double>(i.size()+1+N,0.));
  for(int z=0;z<(int) i.size();z++){
    for(int k=0;k<(int) i.size();k++){
      for(int j=0;j<(int) i.size();j++){
        B[z][j]+=2.*A[k][z]*A[k][j];
      }
    }
    B[z][i.size()]=i[z];
    B[i.size()][z]=i[z];
  }
  int cont=i.size()-1;
  for(int k=0;k<N;k++){
    double sum=0;
    for(int m=0;m<7*NweeksToKeepIncidenceSum;m++){
      if(cont<0) break;
      sum+=i0[cont];
      B[cont][i.size()+1+k]=i[cont];
      B[i.size()+1+k][cont]=i[cont];
      cont--;
    }
    b.push_back(sum);
  }
  
  /// q Regularity ///NEW
  double aux=1e2;
  for(int t=0;t<(int) i.size();t++){
    if(t<7){
      B[t][t]+=aux*seasonality_regularization_weight;
      B[t][t+7]-=aux*seasonality_regularization_weight;
    }
    else if(t>=(int) i.size()-7){
      B[t][t]+=seasonality_regularization_weight;
      B[t][t-7]-=seasonality_regularization_weight;
    }
    else{
      if(t<56){
        B[t][t]+=aux*2.*seasonality_regularization_weight;
        B[t][t-7]-=aux*seasonality_regularization_weight;
        B[t][t+7]-=aux*seasonality_regularization_weight;
      }
      else{
        B[t][t]+=2.*seasonality_regularization_weight;
        B[t][t-7]-=seasonality_regularization_weight;
        B[t][t+7]-=seasonality_regularization_weight;
      }
    }
  }
  
  /// SOLUTION COMPUTATION
  vector<double> u=linear_system_solution(B,b);
  
  if(u.size()==0){
    return vector<double>(i.size(),1.);
  }
  u.resize(u.size()-1); /// we remove from the solution the Lagrange multiplier.
  
  return u;
}

///----------------------------------------------------------------------------------------------------
/// COMPUTATION OF THE EPIESTIM ESTIMATION OF Rt
vector<double> EpiEstim(
    vector<double> &c0 /** incidence curve */,
    double a,double b,int tau /** EPIESTIM parameters*/,
    char si_distr_filename[])
{
  //printf("READING THE SERIAL INTERVAL\n\n");
  vector<double> si_distr;
  int k0=read_si_distr(si_distr_filename,si_distr);
  
  vector<double> c(c0.size(),0.);
  
  for(int k=0;k<(int) c.size();k++){
    for(int m=0;m<tau;m++){
      c[k]+=c0[k>m?k-m:0];
    }
    c[k]/=tau;
  }
  
  vector<double> Epi(c.size(),0.);
  for(int n=0;n<(int) c.size();n++){
    double den=1./(b*tau);
    for(int k=0;k<(int) si_distr.size();k++){
      if(k<=k0){
        den+=si_distr[k]*c[n];
      }
      else{
        int t=n-(k-k0);
        den+=si_distr[k]*c[t>0?t:0];
      }
    }
    Epi[n]=(a/tau+c[n])/den;
  }
  return Epi;
}


///----------------------------------------------------------------------------------------------------
/// Rt AND WEEKLY CORRECTION FACTORS q ESTIMATION
/// RETURNS I (THE RMSE RATIO AFTER WEEKLY BIAS CORRECTION)
double Rt_q_estimation(
    vector<double> &i /** daily infected */,
    vector<double> &Pi /** PARAMETERS OF THE EXTRAPOLATION OF THE INITIAL GROWTH OF THE INCIDENCE CURVE*/,
    const vector<double> &si_distr /** serial interval distribution */,
    const int f0 /** position of the value 0 in the serial interval */,
    const double Rt_regularization_weight /** weigth for the energy minimization*/,
    const vector<double> &nf /** median used to normalize the incidence curve in the energy*/,
    const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
    vector<double> &R /** Rt estimate. We assume it is initialized using the original incidence curve*/,
    vector<double> &q /** weekly bias correction factors */,
    vector<double> &iRenEq /** final incidence curve using the renewal equation and correcting the weekly bias */,
    int &iter_alternate_optimization /** number of iterations of the alternate optimization algorithm to compute R and q*/,
    const int MaxIter /** max number of iterations in the alternate minimization algorithm */,
    bool WeeklyBiasCorrection /** if TRUE the weekly is corrected*/,
    vector<int> &daily_festive_day /** if this vector is not empty we take into account the daily_festive_days when minimizing the energy*/,
    const vector<double> &i0 /** original daily infected */,
    int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/,
    double seasonality_regularization_weight /** weight of the regularization term for the seasonality q */
){
  q=vector<double>(i.size(),1.); /// q initialization
  int tmin=0;
  
  /// FIRST INCIDENCE CURVE OBTAINED USING THE RENEWAL EQUATION WITHOUT WEEKLY CORRECTION FACTOR
  iRenEq=i;
  for(int k=0;k<(int) i.size();k++){
    iRenEq[k]=RenewalEquation(k,si_distr,f0,R,i,Pi,RenewalEquationModel);
  }
  
  /// ESTIMATION OF THE CUMMULATIVE SQUARED ERROR BETWEEN THE ORIGINAL INCIDENCE CURVE AND THE
  /// OBTAINED USING THE RENEWAL EQUATION
  double error0=0.;
  for(int k=tmin;k<(int) i.size();k++)  error0+=(i[k]-iRenEq[k])*(i[k]-iRenEq[k]);
  
  /// ALTERNATE ALGORITHM TO OPTIMIZE R AND q
  double error_min=1e20;
#ifndef R_COMPILE
  //printf("iter : ");
#endif
  for(iter_alternate_optimization=1;iter_alternate_optimization<=MaxIter;iter_alternate_optimization++){
#ifndef R_COMPILE
    //printf("%d,",iter_alternate_optimization);
#endif
    
    /// COMPUTATION WEEKLY BIAS CORRECTION FACTORS
    vector<double> u=vector<double>(i.size(),1.);
    if(WeeklyBiasCorrection==true) u=periodic_7days(i,R,si_distr,f0,Pi,RenewalEquationModel,tmin,i.size()-1,nf,i,NweeksToKeepIncidenceSum,seasonality_regularization_weight);
    vector<double> iBiasCor=i;
    for(int k=i.size()-1;k>=0;k--){
      iBiasCor[k]=i[k]*u[k];
    }
    
    /// WE TAKE INTO ACCOUNT THE DAILY FESTIVE DAYS
    if(daily_festive_day.size()>=i.size()){
      vector<int> daily_festive_day2=daily_festive_day;
      for(int k=0;k<(int) daily_festive_day.size();k++){
        if(daily_festive_day[k]==1 && k+1<(int) daily_festive_day.size()) daily_festive_day2[k+1]=1;
        if(daily_festive_day[k]==1 && k+2<(int) daily_festive_day.size()) daily_festive_day2[k+2]=1;
        //if(daily_festive_day[k]==1 && k+3<(int) daily_festive_day.size()) daily_festive_day2[k+3]=1; ///NEW
        //if(daily_festive_day[k]==1 && k+4<(int) daily_festive_day.size()) daily_festive_day2[k+4]=1; ///NEW
      }
      
      LinearSystem_it(i,iBiasCor,u,si_distr,f0,RenewalEquationModel,R,nf,daily_festive_day2);
      for(int k=0;k<(int) i.size();k++){
        if(daily_festive_day2[k]==1){
          i[k]=iBiasCor[k]/u[k];
        }
      }
    }
    
    /// UPDATE OF THE INCIDENCE CURVE USING THE RENEWAL EQUATION
    vector<double> iRenEqNew=i;
    for(int k=0;k<(int) i.size();k++){
      iRenEqNew[k]=RenewalEquation(k,si_distr,f0,R,iBiasCor,Pi,RenewalEquationModel);
    }
    
    /// COMPUTING RMSE RATIO AFTER WEEKLY BIAS CORRECTION ERROR
    double error1=0.;
    for(int k=tmin;k<(int) i.size();k++) error1+=(iBiasCor[k]-iRenEqNew[k])*(iBiasCor[k]-iRenEqNew[k]);
    double error2=sqrt(error1/error0);
    
    /// CHECKING IF RMSE RATIO AFTER WEEKLY BIAS CORRECTION ERROR IS REDUCED. OTHERWISE WE STOP ITERATIONS
    if(error2<error_min){
      error_min=error2;
      q=u;
      q.resize(i.size());
      iRenEq=iRenEqNew;
    }
    else break;
    
    /// R UPDATE USING THE NEW WEEKLY BIAS CORRECTION FACTORS
    Rt_estimation(iBiasCor,Pi,si_distr,f0,Rt_regularization_weight,nf,RenewalEquationModel,R,iBiasCor,NweeksToKeepIncidenceSum);
  }
#ifndef R_COMPILE
  //printf("\n");
#endif
  return error_min;
}


///----------------------------------------------------------------------------------------
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
    const double mean_si /** MEAN OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double sd_si /** STANDARD DEVIATION OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double shift_si /** SHIFT OF THE LOG-NORMAL PARAMETRIC SERIAL INTERVAL */,
    const double Rt_regularization_weight /** REGULARIZATION WEIGHT PARAMETER OF EpiInvert METHOD (DEFAULT VALUE: 5)*/,
    const double seasonality_regularization_weight /** WEIGHT PARAMETER OF THE REGULARIZATION  TERM FOR THE SEASONALITY q (DEFAULT VALUE 5) */,
    const int max_time_interval /** MAX SIZE OF THE INCIDENCE DATA USED TO COMPUTE Rt (DEFAULT VALUE: 9999). THIS PARAMETER IS USED TO REDUCE HE COMPUTATIONAL COST OF THE ALGORITHM WHEN WE ARE JUST INTERESTED IN THE LAST PART OF THE SEQUENCE */,
    const int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/,
    bool weekly_aggregated_incidence /** IF TRUE, EACH INCIDENCE VALUE CORRESPONDS TO THE LAST 7-DAY AGGREGATED INCIDENCE */
){
  /// WE MANAGE THE CASE OF WEEKLY AGGREGATED INCIDENCE
  if(weekly_aggregated_incidence==true){
    vector<double> i_original2(7*i_original.size());
    for(int k=0,m=0;k<(int) i_original.size();k++){
      double value=i_original[k]/7.;
      for(int n=0;n<7;n++) i_original2[m++]=value;
    }
    i_original=i_original2;
  }
  
  /// RENEWAL EQUATION MODEL
  bool RenewalEquationModel = CASE;
  
  /// WEEKLY BIAS CORRECTION
  bool WeeklyBiasCorrection = true;
  
  /// SERIAL INTERVAL
  if(si_distr.size()==0){
    shift_si_np = parametric_si_distr(mean_si,sd_si,shift_si,si_distr);
  }
  int k0=-shift_si_np;
  
  /// VECTOR TO STORE ALL INITIAL DATA SEQUENCES CORRESPONDING TO SEVERAL COUNTRIES (IN GENERAL WE USE A SINGLE COUNTRY BUT WE CAN ALSO ACCUMULATE THE DATA OF SEVERAL CONTRIES)
  vector< vector<double> >  iV0;
  
  /// WE GET THE INCIDENCE DATA
  iV0.push_back(i_original);
  
#ifndef R_COMPILE
  //printf("-> NUMBER OF SAMPLES OF THE ORIGINAL INCIDENCE CURVE : %d \n",(int) i_original.size());
  //printf("-> PRE-PROCESSING OF THE ORIGINAL INCIDENCE CURVE\n");
#endif
  
  i_original=data_pre_processing(iV0,0,max_time_interval);
  
#ifndef R_COMPILE
  //printf("-> NUMBER OF SAMPLES OF THE ORIGINAL INCIDENCE CURVE AFTER PRE-PROCESSING: %d \n",(int) i_original.size());
  //printf("-> NUMBER OF FESTIVE OR ANOMALOUS DAYS USED : %d \n",(int) festive_days.size());
#endif
  
  
  i_festive=i_original;
  
  
  /// CURRENT DAY
  time_t current_day = string2date(last_incidence_date.c_str());
  
  ///FROM THE FESTIVE DAY SEQUENCE WE BUILT AN INDEX VECTOR WHICH ASSOCIATES THE VALUE 1 TO A FESTIVE DAY AND 0 OTHERWISE
  vector<int> daily_festive_day=daily_festive_day_initialization(current_day, i_festive.size(),festive_days);
  
  /// NUMBER OF DAYS IN THE PAST USED TO COMPUTE THE Rt VARIABILITY
  int NdaysEmpiricalVariability=4;
  
  /// WE ESTIMATE THE FULL Rt SEQUENCE NdaysEmpiricalVariability TIMES. EACH TIME WE REMOVE AN ADDITIONAL DAY TO THE SEQUENCE
  vector< vector<double> > RtV(NdaysEmpiricalVariability); /// RtV[m] IS THE Rt ESTIMATE REMOVING m DAYS OF THE DATA
  vector<double> RtLastV(NdaysEmpiricalVariability); /// Last value of RtV[m]
  vector< vector<double> > iV(NdaysEmpiricalVariability);  /// iV[m] IS THE FINAL VERSION OF INCIDENCE CURVE USING THE RENEWAL EQUATION AND WEEKLY BIAS CORRECTION
  vector< vector<double> > qV(NdaysEmpiricalVariability,vector<double>(7,1.)); ///qV[m] IS THE VECTOR WITH THE PERIODIC BIAS CORRECTION FACTOR
  vector<double> RMSE_factorV(NdaysEmpiricalVariability,1e20);  ///RMSE_factorV[m] is the RMSE reduction in the renewal equation error after weekly bias correction
  vector<int> iter_alternate_optimizationV(NdaysEmpiricalVariability); /** number of iterations of the alternate optimization algorithm to compute RtV[m] and qV[m]*/
  
  /// COMPUTATION OF THE PERCENTIL IN THE LAST 21 DAYS TO NORMALIZE THE INCIDENCE CURVE
  vector<double> nf=back_percentil(i_festive,21);
  
#ifndef R_COMPILE
  //printf("-> EpiInvert COMPUTATION\n");
#endif
  
  
  /// PARALLEL COMPUTATION OF EpiInvert FOR THE LAST "NdaysEmpiricalVariability" DAYS  ///TTT
#ifdef _OPENMP
#pragma omp parallel for shared(i_festive,iV0,RtV,RtLastV,iV,qV,RMSE_factorV,iter_alternate_optimizationV,si_distr,shift_si_np,nf)
#endif
  for(int m=0;m<NdaysEmpiricalVariability;m++){
#ifndef R_COMPILE
    //printf("  -> current day - %d : ",m);
#endif
    
    /// LOCAL VARIABLES
    vector<double> i,R,q(7,1.),i0n;
    int iter_alternate_optimization;
    double RMSE_factor=1e20;
    vector<double> Pn;
    
    /// WE USE THE DATA UP TO m DAYS IN ADVANCE WITH RESPECT TO THE LAST AVAILABLE DAY (USED TO COMPUTE THE Rt VARIABILITY)
    i0n=m==0?i_festive:data_pre_processing(iV0,m,max_time_interval);
    
    ///FROM THE FESTIVE DAY SEQUENCE WE BUILT AN INDEX VECTOR WHICH ASSOCIATES THE VALUE 1 TO A FESTIVE DAY AND 0 OTHERWISE
    vector<int> daily_festive_day2=daily_festive_day_initialization(current_day-m*86400, i0n.size(),festive_days);
    
    /// EpiInvert COMPUTATION AFTER REMOVING m DAYS FROM THE ORIGINAL SEQUENCE
    Rt_estimation(i0n,Pn,si_distr,k0,Rt_regularization_weight,nf,RenewalEquationModel,R,i0n,NweeksToKeepIncidenceSum); ///TTT
    
    /// MANAGEMENT OF WEEKLY BIAS CORRECTION AND FESTIVE DAYS
    if(WeeklyBiasCorrection==true || daily_festive_day2.size()>0){
      ///printf("ITERATIVE OPTIMIZATION OF Rt, THE WEEKLY CORRECTION FACTORS AND THE INCIDENCE IN THE FESTIVE DAYS\n");
      RMSE_factor=Rt_q_estimation(i0n,Pn,si_distr,k0,Rt_regularization_weight,nf,RenewalEquationModel,
                                  R,q,i,iter_alternate_optimization,(m==0?20:2),WeeklyBiasCorrection,daily_festive_day2,i0n,NweeksToKeepIncidenceSum,seasonality_regularization_weight);
      
      
      /// WE UPDATE THE INCIDENCE CURVE IN THE CASE FESTIVE DAYS ARE USED
      if(m==0 && daily_festive_day2.size()>0 ) i_festive=i0n;
    }
    
    /// OTHERWISE WE COMPUTE MANUALLY THE RESTORED INCIDENCE CURVE
    else{
      i=vector<double>(i0n.size());
      for(int k=0;k<(int) i.size();k++){
        i[k]=RenewalEquation(k,si_distr,k0,R,i0n,Pn,RenewalEquationModel);
      }
    }
    
    ///printf("WE STORE EpiInvert RESULTS FOR ITERATION m\n");
    RtLastV[m]=R[R.size()-1];
    RMSE_factorV[m]=RMSE_factor;
    RtV[m]=R;
    qV[m]=q;
    iV[m]=i;
    iter_alternate_optimizationV[m]=iter_alternate_optimization;
  }
  
  //printf("EXTRAPOLATION UP TO THE CURRENT DATE OF RtV[m] (TO EVALUATE THE EMPIRIC VARIABILITY UP TO THE CURRENT DATE)\n");
  for(int m=1;m<NdaysEmpiricalVariability;m++){
    double deriv=RtLastV[m]-RtV[m][RtV[m].size()-2];
    for(int k=1;RtV[m].size()<RtV[0].size();k++)  RtV[m].push_back(RtLastV[m]+k*deriv);
  }
  
  /// COMPUTATION OF THE Rt VARIABILITY UP TO THE CURRENT DATE
  vector<double> V = vector<double>(RtV[0].size(),0.);
  for(int k=0;k<(int) V.size();k++){
    for(int m=1;m<NdaysEmpiricalVariability;m++) V[k]+=(RtV[0][k]-RtV[m][k])*(RtV[0][k]-RtV[m][k]);
    V[k]=sqrt(V[k]/(NdaysEmpiricalVariability-1.));
  }
  
  /// POST-PROCESSING V
  {
    int N=V.size()-1;
    while(N>0 && V[N]>V[N-1]) N--;
    for(int k=0;k<N;k++) V[k]=0;
  }
  
  /// EMPIRICAL PRECOMPUTED VALUES TO DEFINE EMPIRIC CONFIDENCE INTERVALS FOR Rt
  Rt_CI95=V;
  double B95=0.24,C95=0.03;
  
  for(int k=0;k<(int) V.size();k++){
    double aux95=B95+C95*(k-(V.size()-1.));
    if(aux95>0) Rt_CI95[k]+=aux95;
  }
  
  
  /// WE SET THE OUTPUT VARIABLES
  i_restored=iV[0];
  Rt=RtV[0];
  seasonality=qV[0];
  RMSE_factor=RMSE_factorV[0];
  iter_alternate_optimization=iter_alternate_optimizationV[0];
  i_bias_free=vector<double>(i_festive.size());
  for(int k=0;k< (int) i_bias_free.size();k++) i_bias_free[k]=i_festive[k]*seasonality[k];
  
  //for(int k=0;k<i_restored.size();k++) printf("%lf\n",i_restored[k]);
  
  /// ERROR ANALYSIS
  vector<double> log_i_rest, log_abs_dif;
  for(int k=21;k< (int) i_bias_free.size();k++){
    if(i_restored[k]<100) continue;
    log_i_rest.push_back(log(i_restored[k]));
    log_abs_dif.push_back(log(fabs(i_bias_free[k]-i_restored[k])));
  }
  double b;
  linear_regression(log_i_rest,log_abs_dif,a,b);
  //printf("a=%lf  b=%lf r=%lf size=%d\n",a,b,r, (int) log_i_rest.size());
  
  vector<double> x,y,error(log_i_rest.size());
  for(int k=0;k< (int) error.size();k++){ error[k]=fabs(a*log_i_rest[k]+b-log_abs_dif[k]);}
  int perc=error.size()*0.95;
  double Q95=percentil(perc,error);
  for(int k=0;k< (int) error.size();k++){
    if(error[k]<Q95){
      x.push_back(log_i_rest[k]);
      y.push_back(log_abs_dif[k]);
    }
  }
  linear_regression(x,y,a,b);
  
  epsilon=vector<double>(i_bias_free.size(),0.);
  
  for(int k=0;k< (int) i_bias_free.size();k++){
    if(i_restored[k]<=0) continue;
    epsilon[k]=(i_bias_free[k]-i_restored[k])/pow(i_restored[k],a);
  }
  
  dates=vector<string>(i_original.size());
  festive=vector<bool>(i_original.size());
  //printf("daily_festive_day.size()=%d, festive.size()=%d\n",(int) daily_festive_day.size(),(int) festive.size());
  for(int k=0;k<(int) i_original.size();k++){
    time_t t2=current_day-(i_original.size()-1-k)*86400+86400/2;
    struct tm * timeinfo;
    timeinfo = localtime (&t2);
    char buffer [80];
    strftime (buffer,80,"%Y-%m-%d",timeinfo);
    dates[k]=string(buffer);
    festive[k]=daily_festive_day[k]==0?false:true;
  }
}

/// ----------------------------------------------------------------------------------------
/// FORECAST OF THE RESTORED INCIDENCE USING A LEARNING PROCEDURE
/// ----------------------------------------------------------------------------------------
vector<double> IncidenceForecastByLearning(
    vector<double> &ir /** RESTORED INCIDENCE TO BE FORECASTED */,
    const string last_incidence_date /** DATE OF THE LAST DATA IN THE FORMAT YYYY-MM-DD */,
    vector<double> &q /** 7-DAY QUASI-PERIODIC WEKLY BIAS CORRECTION FACTORS */,
    vector< vector <double > > &ir_database,
    double lambda /** PARAMETER IN THE WEIGHTED AVERAGE OF INCIDENCE SEQUENCES (IF NEGATIVE WE USE THE PRE-ESTIMATED ONES */,
                                                                                double mu /** PARAMETER TO ESTIMATE THE DISTANCE BETWEEN RESTORED INCIDENCE CURVES*/,
                                                                                vector <double> &CI50 /** 50% CONFIDENCE INTERVAL RADIUS FOR THE FORECAST OF THE RESTORED INCIDENCE */,
                                                                                vector <double> &CI75 /** 75% CONFIDENCE INTERVAL RADIUS FOR THE FORECAST OF THE RESTORED INCIDENCE */,
                                                                                vector <double> &CI90 /** 90% CONFIDENCE INTERVAL RADIUS FOR THE FORECAST OF THE RESTORED INCIDENCE */,
                                                                                vector <double> &CI95 /** 95% CONFIDENCE INTERVAL RADIUS FOR THE FORECAST OF THE RESTORED INCIDENCE */,
                                                                                vector <double> &i0_forecast /** FORECAST OF THE ORIGINAL INCIDENCE */,
                                                                                vector<string> &dates /** DATE ASSOCIATED TO EACH INCIDENCE DATUM */
)
{
  /// WE CHECK THE INPUT
  if( ir.size()<28 || q.size()!=ir.size() ) return vector<double>();
  
  /// REFERENCE 50% CONFIDENCE INTERVAL RADIUS FOLLOWING THE FORECAST DAY
  double c50[28]={0.042508,0.054770,0.068168,0.081978,0.096363,0.111018,0.125447,0.139321,0.153157,0.167250,0.181470,0.195171,0.209277,0.223595,0.237681,0.251021,
                  0.265775,0.279932,0.296414,0.312585,0.328543,0.345365,0.362289,0.380731,0.397672,0.415894,0.434526,0.453136};
  /// REFERENCE 75% CONFIDENCE INTERVAL RADIUS FOLLOWING THE FORECAST DAY
  double c75[28]={0.081978,0.105534,0.131015,0.157528,0.184099,0.211136,0.236895,0.261647,0.286255,0.310039,0.332936,0.355424,0.377950,0.399522,0.420463,0.441227,0.464055,
                  0.486819,0.508160,0.528966,0.550036,0.570330,0.591685,0.611198,0.630934,0.650795,0.671670,0.690735};
  /// REFERENCE 90% CONFIDENCE INTERVAL RADIUS FOLLOWING THE FORECAST DAY
  double c90[28]={0.139130,0.178687,0.222212,0.268017,0.313669,0.358350,0.401996,0.444749,0.485644,0.524075,0.560992,0.592548,0.620591,0.648606,0.674875,0.699695,0.720586,
                  0.744659,0.764860,0.784399,0.805348,0.825345,0.843817,0.865360,0.885556,0.903012,0.917678,0.936678};
  /// REFERENCE 95% CONFIDENCE INTERVAL RADIUS FOLLOWING THE FORECAST DAY
  double c95[28]={0.194120,0.246874,0.305963,0.370376,0.438343,0.503850,0.568111,0.629777,0.691450,0.750809,0.812840,0.878160,0.933841,0.987263,1.044576,1.097479,1.164950,1.233489,
                  1.287934,1.344992,1.422996,1.504297,1.587260,1.643302,1.721445,1.784263,1.845480,1.906459};
  
  /// NORMALIZATION OF THE LAST 28 DAYS OR THE RESTORED INCIDENCE
  vector<double> u(28);
  double sum=0;
  for(int k=0;k<28;k++){
    u[k]=ir[ir.size()-28+k];
    sum+=u[k];
  }
  sum/=28;
  double u0=u[27];
  for(int k=0;k<(int) u.size();k++) u[k]/=sum;
  
  /// 28-DAY FORECAST OF THE RESTORED INCIDENCE USING A WEIGTHED AVERAGE OF THE RESTORED INCIDENCE DATABASE
  vector<double> v(56,0.);
  for(int k=0;k<(int) ir_database.size();k++){
    /// COMPUTATION OF THE AVERAGE DISTANCE BETWEEN u AND THE DATABASE INSTANCE
    double L1=0;
    for(int m=0;m<28;m++) L1+=fabs(u[m]-ir_database[k][m])*exp(-mu*(27-m));
    L1/=28;
    /// COMPUTATION OF THE WEIGHT
    double w=exp(-lambda*L1);
    /// COMPUTATION OF THE CONTRIBUTION OF THIS DATABASE INSTANCE TO THE FORECAST
    for(int m=0;m<(int) v.size();m++) v[m]+=w*ir_database[k][m];
  }
  
  /// WE SCALE THE FORECAST TO u[27]=v[27] (THAT IS, WE IMPOSE CONTINUITY IN THE FORECAST ESTIMATION)
  double v0=v[27];
  for(int m=0;m<(int) v.size();m++) v[m]*=u0/v0;
  
  /// WE STORE THE RESULTS
  vector<double> ir_forecast(28);
  i0_forecast=vector<double>(28);
  CI50=vector<double>(28);
  CI75=vector<double>(28);
  CI90=vector<double>(28);
  CI95=vector<double>(28);
  dates=vector<string>(28);
  
  for(int k=0;k<28;k++){
    ir_forecast[k]=v[k+28];
    i0_forecast[k]=ir_forecast[k]/q[q.size()-7+k%7];
    CI50[k]=c50[k]*ir_forecast[k];
    CI75[k]=c75[k]*ir_forecast[k];
    CI90[k]=c90[k]*ir_forecast[k];
    CI95[k]=c95[k]*ir_forecast[k];
    time_t current_day = string2date(last_incidence_date.c_str());
    time_t t2=current_day + (k+1)*86400+86400/2;
    struct tm * timeinfo;
    timeinfo = localtime (&t2);
    char buffer [80];
    strftime (buffer,80,"%Y-%m-%d",timeinfo);
    dates[k]=string(buffer);
  }
  
  //  /// TESTING
  //  FILE *f;
  //  f=fopen ("forecast.csv", "w");
  //  fprintf(f,"date;ir_forecast;i0_forecast;CI95_ir_forecast -;CI95_ir_forecast +\n");
  //  for(int k=0;k<28;k++){
  //      fprintf(f,";%lf;%lf\n",v[k],ir[ir.size()-28+k]);
  //  }
  //  for(int k=0;k<28;k++){
  //    fprintf(f,"%s;%lf;%lf;%lf;%lf\n",dates[k].c_str(),ir_forecast[k],i0_forecast[k],
  //            ir_forecast[k]-CI95[k],ir_forecast[k]+CI95[k]);
  //  }
  //  fclose(f);
  
  
  return ir_forecast;
  
}

///----------------------------------------------------------------------------------------
/// EpiInvert ESTIMATION. IT RETURNS THE TIME WHERE Rt STARTS TO BE COMPUTED (A NEGATIVE VALUE IN CASE OF FAILURE)
int EpiInvert(
    vector<double> &i0 /** ORIGINAL INCIDENCE CURVE (IF EMPTY WE READ THE INCIDENCE FROM THE FILE incidence_filename*/,
                                                     vector<double> &i1 /** ORIGINAL INCIDENCE CURVE (IF EMPTY WE READ THE INCIDENCE FROM THE FILE incidence_filename*/,
                                                                                                      const char incidence_filename[] /** NAME OF THE FILE WITH THE INCIDENCE DATA (THE REGISTERED DAILY NEW TESTED POSITIVE)*/,
                                                                                                      const char si_distr_filename[] /** NAME OF THE FILE CONTAINING THE SERIAL INTERVAL (DEFAULT VALUE: "Ma:txt") */,
                                                                                                      const double Rt_regularization_weight /** REGULARIZATION WEIGHT PARAMETER OF EpiInvert METHOD (DEFAULT VALUE: 5)*/,
                                                                                                      const bool RenewalEquationModel /** THE VALUE CAN BE THE MACRO CASE (=false) OR THE MACRO INST (=true) */,
                                                                                                      const bool WeeklyBiasCorrection /** IF TRUE WE CORRECT THE WEEKLY BIAS (DEFAULT VALUE: TRUE) */,
                                                                                                      const int max_time_interval /** MAX SIZE OF THE INCIDENCE DATA USED TO COMPUTE Rt (DEFAULT VALUE: 9999). THIS PARAMETER IS USED TO REDUCE HE COMPUTATIONAL COST OF THE ALGORITHM WHEN WE ARE JUST INTERESTED IN THE LAST PART OF THE SEQUENCE */,
                                                                                                      vector<double> &i_renewal_equation /** FILTERED INCIDENCE CURVE WITH WEEKLY BIAS CORRECTION AND APPLICATION OF THE RENEWAL EQUATION*/,
                                                                                                      vector<double> &Rt /** EpiInvert ESTIMATION OF THE REPRODUCTION NUMBER*/,
                                                                                                      vector<double> &seasonality /** VECTOR WITH THE 7-DAY WEEKLY CORRECTION FACTORS (THE LAST VALUE OF seasonality CORRESPONDS TO THE LAST VALUE OF i0)*/,
                                                                                                      vector<double> &V /** VARIABILITY MEASURE OF THE Rt ESTIMATE USING THE LAST 4 DAYS */,
                                                                                                      double &RMSE_factor /** REDUCTION FACTOR OF THE RMSE BETWEEN THE INCIDENCE CURVE AND THE RENEWAL EQUATION BEFORE AND AFTER THE WEEKLY BIAS CORRECTION */,
                                                                                                      int &iter_alternate_optimization /** NUMBER OF ITERATIONS OF THE ALTERNATE ALGORITHM TO COMPUTE Rt AND seasonality */,
                                                                                                      vector <string> &festive_days /** control of incidence day which corresponds to daily_festive_days (==1 daily_festive_dayS, ==0 WORKING DAY )*/,
                                                                                                      time_t &current_day, ///IN INPUT, IT CAN BE : (=-1 TO INDICATE THAT WE USE STORED FESTIVE DAYS )
                                                                                                      /// ( >=0 TO INDICATE THE NUMBER OF DAYS PASSED FROM THE LAST AVAILABLE DATA OF THE LAST FESTIVE DAY)
                                                                                                      /// IN OUTPUT REPRESENTS THE DATE OF THE LAST AVAILABLE INCIDENCE DATA
                                                                                                      int NweeksToKeepIncidenceSum /** WE CONSTRAINT ALL THE ESTIMATED INCIDENCE CURVE TO KEEP THE ADDITION OF THE ORIGINAL INCIDENCE IN INTERVALS OF SIZE NweeksToKeepIncidenceSum*7 DAYS*/,
                                                                                                      double seasonality_regularization_weight /** weight of the regularization term for the seasonality q */
){
#ifndef R_COMPILE
  //printf("-> READING THE SERIAL INTERVAL FROM %s\n",si_distr_filename);
#endif
  
  vector<double> si_distr;
  int f0=read_si_distr(si_distr_filename,si_distr);
  
  /// WE STORE THE VALUE OF THE CURRENT DAY AS INTEGER
  int current_day_int=current_day;
  
  /// VECTOR TO STORE ALL INITIAL DATA SEQUENCES CORRESPONDING TO SEVERAL COUNTRIES (IN GENERAL WE USE A SINGLE COUNTRY BUT WE CAN ALSO ACCUMULATE THE DATA OF SEVERAL CONTRIES)
  vector< vector<double> >  iV0;
  
  /// WE GET THE INCIDENCE DATA
  if(i0.size()>0) iV0.push_back(i0); /// the incidence is directly provided as parameter
  else{
#ifndef R_COMPILE
    //printf("-> READING THE DATA FROM %s\n",incidence_filename);
#endif
    
    iV0=read_data_multiple(incidence_filename,current_day);
    if(iV0.size()==0 || iV0[0].size()<25) return -1;
  }
  
  /// PRE-PROCESSING OF THE INITIAL INCIDENCE CURVE
  i0=data_pre_processing(iV0,0,max_time_interval);
  
  i1=i0;
  
  /// WE MANAGE THE FESTIVE DAYS IN THE CASE IT IS NOT PROVIDED BY A USER FILE
  if(current_day_int==-1){ /// we use an stored sequence of festive days (only implemented for the USA, France, Germany and Spain.
    vector <string> festive_days2=get_stored_festive_days(iV0[0]);
    for(int k=0;k< (int) festive_days2.size();k++) festive_days.push_back(festive_days2[k]);
  }
  //  if(festive_days.size()==0){
  //    if(current_day_int>=0){ /// we use as single festive day "current_day_int" days before the current day.
  //      time_t t2=current_day-current_day_int*86400;
  //      struct tm * timeinfo;
  //      timeinfo = localtime (&t2);
  //      char buffer [80];
  //      strftime (buffer,80,"%Y-%m-%d",timeinfo);
  //      festive_days.clear();
  //      string st(buffer);
  //      festive_days.push_back(st);
  //    }
  //    else if(current_day_int==-1){ /// we use an stored sequence of festive days (only implemented for the USA, France, Germany and Spain.
  //      festive_days.clear();
  //      festive_days=get_stored_festive_days(iV0[0]);
  //    }
  //  }
  
  ///FROM THE FESTIVE DAY SEQUENCE WE BUILT AN INDEX VECTOR WHICH ASSOCIATES THE VALUE 1 TO A FESTIVE DAY AND 0 OTHERWISE
  vector<int> daily_festive_day=daily_festive_day_initialization(current_day, i0.size(),festive_days);
  
  /// NUMBER OF DAYS IN THE PAST USED TO COMPUTE THE Rt VARIABILITY
  int NdaysEmpiricalVariability=4;
  
  /// WE ESTIMATE THE FULL Rt SEQUENCE NdaysEmpiricalVariability TIMES. EACH TIME WE REMOVE AN ADDITIONAL DAY TO THE SEQUENCE
  vector< vector<double> > RtV(NdaysEmpiricalVariability); /// RtV[m] IS THE Rt ESTIMATE REMOVING m DAYS OF THE DATA
  vector<double> RtLastV(NdaysEmpiricalVariability); /// Last value of RtV[m]
  vector< vector<double> > iV(NdaysEmpiricalVariability);  /// iV[m] IS THE FINAL VERSION OF INCIDENCE CURVE USING THE RENEWAL EQUATION AND WEEKLY BIAS CORRECTION
  vector< vector<double> > qV(NdaysEmpiricalVariability,vector<double>(7,1.)); ///qV[m] IS THE VECTOR WITH THE PERIODIC BIAS CORRECTION FACTOR
  vector<double> RMSE_factorV(NdaysEmpiricalVariability,1e20);  ///RMSE_factorV[m] is the RMSE reduction in the renewal equation error after weekly bias correction
  vector<int> iter_alternate_optimizationV(NdaysEmpiricalVariability); /** number of iterations of the alternate optimization algorithm to compute RtV[m] and qV[m]*/
  
  /// COMPUTATION OF THE PERCENTIL IN THE LAST 21 DAYS TO NORMALIZE THE INCIDENCE CURVE
  vector<double> nf=back_percentil(i0,21);
#ifndef R_COMPILE
  //printf("-> EpiInvert COMPUTATION\n");
#endif
  
  
  /// PARALLEL COMPUTATION OF EpiInvert FOR THE LAST "NdaysEmpiricalVariability" DAYS  ///TTT
#ifdef _OPENMP
#pragma omp parallel for shared(i0,iV0,RtV,RtLastV,iV,qV,RMSE_factorV,iter_alternate_optimizationV,si_distr,f0,nf)
#endif
  for(int m=0;m<NdaysEmpiricalVariability;m++){
    /// LOCAL VARIABLES
    vector<double> i,R,q(7,1.),i0n;
    int iter_alternate_optimization;
    double RMSE_factor=1e20;
    vector<double> Pn;
    
    /// WE USE THE DATA UP TO m DAYS IN ADVANCE WITH RESPECT TO THE LAST AVAILABLE DAY (USED TO COMPUTE THE Rt VARIABILITY)
    i0n=m==0?i0:data_pre_processing(iV0,m,max_time_interval);
    
    ///FROM THE FESTIVE DAY SEQUENCE WE BUILT AN INDEX VECTOR WHICH ASSOCIATES THE VALUE 1 TO A FESTIVE DAY AND 0 OTHERWISE
    vector<int> daily_festive_day2=daily_festive_day_initialization(current_day-m*86400, i0n.size(),festive_days);
    
    /// EpiInvert COMPUTATION AFTER REMOVING m DAYS FROM THE ORIGINAL SEQUENCE
    
    /// FIRST Rt ESTIMATE
    Rt_estimation(i0n,Pn,si_distr,f0,Rt_regularization_weight,nf,RenewalEquationModel,R,i0n,NweeksToKeepIncidenceSum);
    
    /// MANAGEMENT OF WEEKLY BIAS CORRECTION AND FESTIVE DAYS
    if(WeeklyBiasCorrection==true || daily_festive_day2.size()>0){
      /// ITERATIVE OPTIMIZATION OF Rt, THE WEEKLY CORRECTION FACTORS AND THE INCIDENCE IN THE FESTIVE DAYS
      RMSE_factor=Rt_q_estimation(i0n,Pn,si_distr,f0,Rt_regularization_weight,nf,RenewalEquationModel,
                                  R,q,i,iter_alternate_optimization,(m==0?20:2),WeeklyBiasCorrection,daily_festive_day2,i0n,NweeksToKeepIncidenceSum,seasonality_regularization_weight);
      
      
      /// WE UPDATE THE INCIDENCE CURVE IN THE CASE FESTIVE DAYS ARE USED
      if(m==0 && daily_festive_day2.size()>0 ) i0=i0n;
    }
    
    /// OTHERWISE WE COMPUTE MANUALLY THE RESTORED INCIDENCE CURVE
    else{
      i=vector<double>(i0n.size());
      for(int k=0;k<(int) i.size();k++){
        i[k]=RenewalEquation(k,si_distr,f0,R,i0n,Pn,RenewalEquationModel);
      }
    }
    
    /// WE STORE EpiInvert RESULTS FOR ITERATION m
    RtLastV[m]=R[R.size()-1];
    RMSE_factorV[m]=RMSE_factor;
    RtV[m]=R;
    qV[m]=q;
    iV[m]=i;
    iter_alternate_optimizationV[m]=iter_alternate_optimization;
  }
  
  /// EXTRAPOLATION UP TO THE CURRENT DATE OF RtV[m] (TO EVALUATE THE EMPIRIC VARIABILITY UP TO THE CURRENT DATE)
  for(int m=1;m<NdaysEmpiricalVariability;m++){
    double deriv=RtLastV[m]-RtV[m][RtV[m].size()-2];
    for(int k=1;RtV[m].size()<RtV[0].size();k++)  RtV[m].push_back(RtLastV[m]+k*deriv);
  }
  
  /// COMPUTATION OF THE Rt VARIABILITY UP TO THE CURRENT DATE
  V = vector<double>(RtV[0].size(),0.);
  for(int k=0;k<(int) V.size();k++){
    for(int m=1;m<NdaysEmpiricalVariability;m++) V[k]+=(RtV[0][k]-RtV[m][k])*(RtV[0][k]-RtV[m][k]);
    V[k]=sqrt(V[k]/(NdaysEmpiricalVariability-1.));
  }
  
  /// POST-PROCESSING V
  {
    int N=V.size()-1;
    while(N>0 && V[N]>V[N-1]) N--;
    for(int k=0;k<N;k++) V[k]=0;
  }
  
  /// WE SET THE OUTPUT VARIABLES
  i_renewal_equation=iV[0];
  Rt=RtV[0];
  seasonality=qV[0];
  RMSE_factor=RMSE_factorV[0];
  iter_alternate_optimization=iter_alternate_optimizationV[0];
  
  return 0;
  
}

/// ------------------------------------------------------------------------------------------
/// COMPUTATION OF THE L1 NORM BETWEEN AN INCIDENCE COVID SEQUENCE AND A DATABASE OF SEQUENCES
/// (USED IN THE LEARNIG FORECAST PROCEDURE)
vector <double> incidence_comparison(
    const vector <double> &i /** COVID-19 SEQUENCE */,
    const vector < vector <double> > &iV /** DATABASE OF COVID-19 SEQUENCES*/,
    int COMPARISON_TYPE /** TYPE OF COMPARISON */){
  
  vector <double> L1(iV.size(),0.);
  if(COMPARISON_TYPE==L_ONE){
    for(int k=0;k<(int) iV.size();k++){
      for(int m=0;m<(int) i.size();m++) {
        L1[k]+=fabs(i[m]-iV[k][m]);
      }
    }
  }
  else{
    for(int k=0;k<(int) iV.size();k++){
      double sum1=0,sum2=0,sum=0;
      for(int m=0;m<(int) i.size();m++) {
        sum+=i[m]*iV[k][m];
        sum1+=i[m]*i[m];
        sum2+=iV[k][m]*iV[k][m];
      }
      L1[k]=1.-sum/sqrt(sum1*sum2);
    }
  }
  return L1;
  
}


/// ----------------------------------------------------------------------------------------
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES
void IncidenceExtrapolationByLearning(
    vector<double> &i /** DAILY TESTED POSITIVE BY APPLYING THE RENEWAL EQUATION TO THE WEEKLY BIAS FREE INCIDENCE*/,
    const vector< vector <double > > &i42,
    const vector< vector <double > > &i56,
    const int NweeksToKeepIncidenceSumeeksBackToForeCast /** NUMBER OF WEEKS IN THE PAST OF THE INCIDENCE VALUES TO FORECAST THE INCIDENCE CURVE 14 DAYS*/,
    double lambda /** PARAMETER IN THE WEIGHTED AVERAGE OF INCIDENCE SEQUENCES (IF NEGATIVE WE USE THE PRE-ESTIMATED ONES */,
                                                                                int COMPARISON_TYPE /** SIMILARITY CRITERIUM BETWEEN INCIDENCE CURVE. IT CAN BE CORRELATION OF L1 NORM */,
                                                                                int index_to_remove_from_database /** INDEX TO REMOVE FROM THE DATASET OF INCIDENCE CURVE IN THE LEARNING STEP */
)
{
  /// CHECKING THE NUMBER OF WEEKS USING IN THE PAST TO FORECAST THE INCIDENCE DATA IN THE FUTURE
  if(NweeksToKeepIncidenceSumeeksBackToForeCast<1 || NweeksToKeepIncidenceSumeeksBackToForeCast>6) return;
  
  /// IN THE CASE OF A 14-DAY FORECASTING : OPTIMAL VALUES OF LAMBDA ACCORDINGLY WITH THE NUMBER
  /// OF WEEKS TAKING IN THE PASS AND THE SIMILARITY METHOD USED (CORRELATION OR L1 NORM).
  double lambda_COR[6]={4290.,8540.,4300.,2190.,1520.,3000.};
  double lambda_L1[6]={285.,200.,220.,190.,160.,160.};
  
  /// WE FIX THE VALUE OF lambda
  if(lambda<0){
    if( COMPARISON_TYPE==CORRELATION ) lambda=lambda_COR[NweeksToKeepIncidenceSumeeksBackToForeCast-1];
    else lambda=lambda_L1[NweeksToKeepIncidenceSumeeksBackToForeCast-1];
  }
  
  /// NUMBER OF DAYS USED IN THE FORECAST PROCEDURE
  int Ndays=7*NweeksToKeepIncidenceSumeeksBackToForeCast; /// number of days taking into account in the past
  int Ndays14=7*NweeksToKeepIncidenceSumeeksBackToForeCast+14; /// number of days taking into account in the past + 14 days in the future.
  
  /// MATRIX TO STORE THE DATABASE OF INCIDENCE CURVES
  vector< vector<double> > iV2(i42.size(),vector<double>(Ndays)); /// database of original incidence curve
  vector< vector<double> > iV1(i42.size(),vector<double>(Ndays14)); /// database of incidence curve using data 14 days later.
  
  ///REMOVE INIT POINTS AND INCIDENCE CURVE NORMALIZATION
  for(int m=0;m<(int) iV1.size();m++){
    for(int n=0;n<Ndays;n++){ iV2[m][n]=i42[m][i42[m].size()-Ndays+n];}
    for(int n=0;n<Ndays14;n++){ iV1[m][n]=i56[m][i56[m].size()-Ndays14+n];}
    double sum2=0,sum1=0;
    for(int n=0;n<Ndays;n++){
      sum1+=iV1[m][n];
      sum2+=iV2[m][n];
    }
    
    for(int n=0;n<Ndays;n++){ iV2[m][n]/=sum2;}
    for(int n=0;n<Ndays14;n++){ iV1[m][n]/=sum1;}
  }
  
  /// REMOVE A NEIGHBORHOOD OF AN INDEX IN THE DATABASE (IN THE LEARNING STEP)
  if(index_to_remove_from_database>=0){
    int kmin=index_to_remove_from_database<7?0:index_to_remove_from_database-7;
    int kmax=index_to_remove_from_database>(int) iV2.size()-8?iV2.size()-1:index_to_remove_from_database+7;
    iV2.erase(iV2.begin()+kmin,iV2.begin()+kmax);
    iV1.erase(iV1.begin()+kmin,iV1.begin()+kmax);
  }
  
  /// WE EXTRACT THE LAST "7*NweeksToKeepIncidenceSumeeksBackToForeCast" DAYS FROM THE INCIDENCE CURVE
  vector<double> c(iV2[0].size());
  int dis=i.size()-c.size();
  double sum=0;
  for(int k=0;k<(int) c.size();k++){
    c[k]=i[k+dis];
    sum+=c[k];
  }
  
  /// L1 NORMALIZATION OF THE LAST "7*NweeksToKeepIncidenceSumeeksBackToForeCast" DAYS OF THE INCIDENCE CURVE
  for(int k=0;k<(int) c.size();k++) c[k]/=sum;
  
  /// WE COMPUTE THE SIMILARITY MEASURE OF THE LAST "7*NweeksToKeepIncidenceSumeeksBackToForeCast" DAYS OF THE INCIDENCE CURVE WITH THE DATABASE
  vector <double> L1=incidence_comparison(c,iV2,COMPARISON_TYPE);
  
  /// COMPUTING THE WEIGHT OF EACH INCIDENCE SEQUENCE DATABASE IN THE AVERAGE COMBINATION OF ALL SEQUENCES
  for(int k=0;k<(int) L1.size();k++) L1[k]=exp(-lambda*L1[k]);
  
  /// COMPUTING THE WEIGTHED COMBINATIONS OF THE "7*NweeksToKeepIncidenceSumeeksBackToForeCast + 14" DAY SEQUENCES
  vector <double> in(iV1[0].size(),0.);
  for(int k=1;k<(int) iV1.size();k++){
    for(int m=0;m<(int) in.size();m++) in[m]+=L1[k]*iV1[k][m];
  }
  /// NORMALIZATION OF THE OBTAINED SEQUENCE TO SCALE ITS VALUES TO THE ORIGINAL INCIDENCE CURVE
  double sum2=0;
  for(int m=0;m<(int) iV2[0].size();m++) sum2+=in[m];
  for(int m=0;m<(int) in.size();m++) in[m]*=sum/sum2;
  
  /// THE FINAL EXTRAPOLATION SEQUENCE IS MODIFIED TO GET A SMOOTH TRANSITION WITH THE ORIGINAL INCIDENCE CURVE.
  int BackNdays=7;
  if(BackNdays>=(int) iV2[0].size()) BackNdays=iV2[0].size()-1;
  for(int k=c.size()-BackNdays-1;k<(int) c.size();k++){
    if(k<0 || k+dis<0 ) continue;
    double w=(double)(c.size()-1-k)/BackNdays;
    in[k]=(1.-w)*in[k]+w*i[k+dis];
  }
  for(int k=c.size()-BackNdays-1;k<(int) c.size();k++){
    if(k<0 || k+dis<0 ) continue;
    i[k+dis]=in[k];
  }
  
  /// WE ADD THE EXTRAPOLATED VALUES TO THE INCIDENCE CURVE.
  for(int k=c.size();k<(int) in.size();k++){
    i.push_back(in[k]);
  }
}

/// ------------------------------------------------------------------------------------------
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES (USING THE MEDIAN IN THE LAST 5 WEEKS ESTIMATION)
void IncidenceExtrapolationByLearningMedian5Weeks(
    vector<double> &i /** DAILY TESTED POSITIVE BY APPLYING THE RENEWAL EQUATION TO THE WEEKLY BIAS FREE INCIDENCE*/,
    const vector< vector <double > > &i42 /** database of original incidence curve */,
    const vector< vector <double > > &i56 /** database of incidence curve using data 14 days later*/,
    int COMPARISON_TYPE /** SIMILARITY CRITERIUM BETWEEN INCIDENCE CURVE. IT CAN BE CORRELATION OF L1 NORM */
)
{
  vector < vector <double> > iV(5);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k=0;k<5;k++){
    iV[k]=i;
    IncidenceExtrapolationByLearning(iV[k],i42,i56,k+1,-1.,COMPARISON_TYPE,-1);
  }
  
  for(int m=i.size()-36;m<(int) iV[0].size();m++){
    vector<double> A(5);
    for(int k=0;k<5;k++) A[k]=iV[k][m];
    double p=percentil(2,A);
    if(m<(int) i.size()) i[m]=p;
    else i.push_back(p);
  }
}

/// ------------------------------------------------------------------------------------------
/// 14-DAY INCIDENCE EXTRAPOLATION USING A DATABASE OF COVID-19 SEQUENCES (USING THE MEDIAN IN THE LAST 3 WEEKS ESTIMATION)
void IncidenceExtrapolationByLearningMedian3Weeks(
    vector<double> &i /** DAILY TESTED POSITIVE BY APPLYING THE RENEWAL EQUATION TO THE WEEKLY BIAS FREE INCIDENCE*/,
    const vector< vector <double > > &i42 /** database of original incidence curve */,
    const vector< vector <double > > &i56 /** database of incidence curve using data 14 days later*/,
    int COMPARISON_TYPE /** SIMILARITY CRITERIUM BETWEEN INCIDENCE CURVE. IT CAN BE CORRELATION OF L1 NORM */
)
{
  vector < vector <double> > iV(3);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k=0;k<3;k++){
    iV[k]=i;
    IncidenceExtrapolationByLearning(iV[k],i42,i56,k+1,-1.,COMPARISON_TYPE,-1);
  }
  
  for(int m=i.size()-36;m<(int) iV[0].size();m++){
    vector<double> A(3);
    for(int k=0;k<3;k++) A[k]=iV[k][m];
    double p=percentil(1,A);
    if(m<(int) i.size()) i[m]=p;
    else i.push_back(p);
  }
  
}

#ifndef R_COMPILE
///----------------------------------------------------------------------------------------------------
/// WE WRITE THE OUTPUT FILES
void write_files(
    vector<double> &R /** Rt estimate with EpiInvert*/,
    vector<double> &Epi /** Rt estimate with EpiEstim */,
    vector<double> &c /** final incidence curve using the renewal equation and correcting the weekly bias */,
    vector<double> &i0 /** daily infected */,
    vector<double> &i1 /** daily infected */,
    vector<double> &V /** experimental variability of Rt using the estimation in the last days*/,
    bool WeeklyBiasCorrection /** if true we correct the weekly bias */,
    vector<double> &q /** weekly bias correction factors */,
    double error_min /** RMSE reduction factor after weekly bias correction */ ,
    double MinE /** min RMSE comparing EpiInvert and EpiEstim */,
    double tmin /** optimal shift time with EpiEstim */,
    int iter_alternate_optimization /** number of iterations of aternate algorithm to compute Rt */,
    bool RenewalEquationModel /** renewal equation used */,
    time_t current_day /** last date of available data */,
    vector<string> festive_days /** festive days used*/
)
{
  ///COMPUTE DAILY FESTIVE DAY VECTOR. FOR EACH DAY THE VALUE IS 1 FOR FESTIVE DAYS AND 0 OTHERWISE.
  vector<int> daily_festive_day=daily_festive_day_initialization(current_day, (int) i0.size(),festive_days);
  
  /// WE SAVE THE RESULTS
  
  /// ESTIMATION OF THE INCIDENCE CURVE OSCILLATIONS BEFORE AND AFTER WEEKLY BIAS CORRECTION IN THE LAST 8 WEEKS
  if(WeeklyBiasCorrection==true){
    double factor1=0,factor2=0;
    {
      double sum1=0,sum2=0,sum3=0,sum4=0;
      for(int k=i0.size()-1;k>=(int) i0.size()-1-56 && k>=0;k--){
        sum1+=fabs(i0[k]-i0[k-1]);
        sum2+=i0[k];
        sum3+=fabs(i0[k]*q[k]-i0[k-1]*q[k-1]);
        sum4+=i0[k]*q[k];
      }
      factor1=(sum1/sum2);
      factor2=(sum3/sum4);
      
    }
    
    /// WE SAVE MAIN NUMERIC INDICATORS
    FILE *g;
    g=fopen ("SummaryEpiInvertExecution.txt", "w");
    //g=fopen ("SummaryExecutionOutcomes.txt", "w"); ///DEMO
    if(g==NULL){
      printf("\nProblems opening output file. Maybe it is already opened ?\n\n");
      system("pause");
      return;
    }
    fprintf(g,"RMSE Reduction factor I :  %1.3lf\n",error_min);
    fprintf(g,"Variability measure of the original incidence curve  V(i)                                    :  %1.3lf\n",factor1);
    fprintf(g,"Variability measure of the incidence curve after weekly bias correction  V(i_bias_corrected) :  %1.3lf\n",factor2);
    fprintf(g,"Number of iterations of the alternate optimization algorithm to compute Rt and the weekly bias correction : %d \n",iter_alternate_optimization);
    fprintf(g,"Optimal shift with EpiEstim :  %1.2lf\n",tmin);
    fprintf(g,"RMSE comparing with EpiEstim  :  %1.3lf\n",MinE);
    fprintf(g,"EpiEstim last Rt value:  %1.3lf\n",Epi[Epi.size()-1]);
    fprintf(g,"EpiInvert last Rt value: %1.3lf\n",R[R.size()-1]);
    fprintf(g,"\n");
    fclose(g);
    printf("RMSE Reduction factor I :  %1.3lf\n",error_min);
    printf("Variability measure of the original incidence curve  V(i)                                    :  %1.3lf\n",factor1);
    printf("Variability measure of the incidence curve after weekly bias correction  V(i_bias_corrected) :  %1.3lf\n",factor2);
    printf("Number of iterations of the alternate optimization algorithm to compute Rt and the weekly bias correction : %d \n",iter_alternate_optimization);
    printf("Optimal shift with EpiEstim :  %1.2lf\n",tmin);
    printf("RMSE comparing with EpiEstim  :  %1.3lf\n",MinE);
    printf("EpiEstim last Rt value:  %1.3lf\n",Epi[Epi.size()-1]);
    printf("EpiInvert last Rt value: %1.3lf\n",R[R.size()-1]);
  }
  else{
    FILE *g;
    g=fopen ("SummaryExecutionOutcomes.txt", "w");
    fprintf(g,"Optimal shift with EpiEstim :  %1.2lf\n",tmin);
    fprintf(g,"RMSE comparing with EpiEstim :  %1.3lf\n",MinE);
    fprintf(g,"EpiEstim last Rt value :  %1.3lf\n",Epi[Epi.size()-1]);
    fprintf(g,"EpiInvert last Rt value:  %1.3lf\n",R[R.size()-1]);
    
    fprintf(g,"\n");
    fclose(g);
    
    printf("Optimal shift with EpiEstim :  %1.2lf\n",tmin);
    printf("RMSE comparing with EpiEstim :  %1.3lf\n",MinE);
    printf("EpiEstim last Rt value :  %1.3lf\n",Epi[Epi.size()-1]);
    printf("EpiInvert last Rt value:  %1.3lf\n",R[R.size()-1]);
    
  }
  
  
  /// EMPIRICAL PRECOMPUTED VALUES TO DEFINE EMPIRIC CONFIDENCE INTERVALS FOR Rt
  vector<double> V95=V,V90=V;
  double B95=0.24,C95=0.03,B90=0.16,C90=0.022;
  if(RenewalEquationModel==true){
    B95=0.04; C95=0.016; B90=0.02; ;C90=0.009;
  }
  
  for(int k=0;k<(int) V.size();k++){
    double aux95=B95+C95*(k-(V.size()-1.));
    double aux90=B90+C90*(k-(V.size()-1.));
    if(aux90>0) V90[k]+=aux90;
    if(aux95>0) V95[k]+=aux95;
  }
  
  /// EMPIRICAL PRECOMPUTED VALUES (IN PERCENTAGES) TO DEFINE INCIDENCE CONFIDENCE INTERVALS
  double Forecat_Rt_CI95[14]={0.197157,0.238859,0.28364,0.330833,0.372849,0.421777,0.470726,0.517924,0.563459,0.610647,0.656392,0.702293,0.75145,0.792208};
  double Forecat_CI90[14]={0.147196,0.179847,0.21562,0.249973,0.284236,0.319154,0.354655,0.392196,0.427706,0.458855,0.494011,0.532003,0.566549,0.601516};
  
  /// WE SAVE THE SEQUENCES OF VALUES FOR R, Epi, c, c0 AND V.
  FILE *g;
  g=fopen ("Rt.csv", "w");
  if(g==NULL){
    printf("\nProblems opening Rt.csv. Maybe it is already opened ?\n\n");
    system("pause");
    return;
  }
  
  fprintf(g,"q seasonality;date;festive day;");
  /// DEMO
  if(WeeklyBiasCorrection==false){
    fprintf(g,"Rt EpiInvert;Rt EpiEstim;original incidence;festive day bias corrected;renewal equation;Rt variability;CI 95;CI 90");
    if(c.size()==i0.size()) fprintf(g,"\n");
    else fprintf(g,";CI 95 INCIDENCE;CI 90 INCIDENCE\n");
  }
  else {
    fprintf(g,"Rt EpiInvert;Rt EpiEstim;original incidence;festive day bias corrected;weekly bias correction;renewal equation;Rt variability;CI 95;CI 90");
    if(c.size()==i0.size()) fprintf(g,"\n");
    else fprintf(g,";CI 95 INCIDENCE;CI 90 INCIDENCE\n");
  }
  
  for(int k=0,n=0;k<(int) c.size() || k<(int) i0.size();k++){
    if(k<(int) i0.size()) fprintf(g,"%lf;",q[k]);
    else fprintf(g,";");
    if(current_day<1) fprintf(g,";;");
    else{
      time_t t2=current_day-(i0.size()-1-k)*86400+86400/2;
      struct tm * timeinfo;
      timeinfo = localtime (&t2);
      char buffer [80];
      strftime (buffer,80,"%Y-%m-%d",timeinfo);
      fprintf(g,"%s;",buffer);
      if(k>=(int) daily_festive_day.size() || daily_festive_day[k]==0) fprintf(g,"NO;");
      else fprintf(g,"YES;");
    }
    if(k<(int) i0.size() && k<(int) c.size()){
      if(k<(int) V.size()){
        if(WeeklyBiasCorrection==false) fprintf(g,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",R[k],Epi[k],i1[k],i0[k],c[k],V[k],V95[k],V90[k]);
        else fprintf(g,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",R[k],Epi[k],i1[k],i0[k],i0[k]*q[k],c[k],V[k],V95[k],V90[k]);
      }
      else{
        if(WeeklyBiasCorrection==false) fprintf(g,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",R[k],Epi[k],i1[k],i0[k],c[k],V[V.size()-1],V95[V.size()-1],V90[V.size()-1]);
        else fprintf(g,"%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",R[k],Epi[k],i1[k],i0[k],i0[k]*q[k],c[k],V[V.size()-1],V95[V.size()-1],V90[V.size()-1]);
      }
    }
    else if(k<(int) i0.size()){
      fprintf(g,";;%lf;%lf\n",i1[k],i0[k]);
    }
    else{
      if(k<(int) V.size()){
        if(WeeklyBiasCorrection==false) fprintf(g,"%lf;%lf;;;%lf;%lf\n",R[k],Epi[k],c[k],V[k]);
        else fprintf(g,"%lf;%lf;;;;%lf;%lf\n",R[k],Epi[k],c[k],V[k]);
      }
      else{
        if(WeeklyBiasCorrection==false){
          if(n<14) fprintf(g,"%lf;;;;%lf;;;;%lf;%lf;\n",R[k],c[k],c[k]*Forecat_Rt_CI95[n],c[k]*Forecat_CI90[n]);
          else fprintf(g,"%lf;;;;%lf;\n",R[k],c[k]);
          n++;
        }
        else{
          if(n<14) fprintf(g,"%lf;;;;;%lf;;;;%lf;%lf;\n",R[k],c[k],c[k]*Forecat_Rt_CI95[n],c[k]*Forecat_CI90[n]);
          else fprintf(g,"%lf;;;;;%lf\n",R[k],c[k]);
          n++;
        }
      }
    }
  }
  
  fclose(g);
  
  /// DEMO
  
  /// WE WRITE IN AN AUXILIARY FILE THE OPTIMAL SHIFT WITH EPIESTIM
  // g=fopen ("ShiftEpiEstim.txt", "w");
  // fprintf(g,"%d\n", (int) (tmin+0.75));
  //fclose(g);
}

#endif
