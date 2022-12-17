/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file utilities.cpp
 * @brief functions for the EpiIndicators method to compute the
 *        functional relationship between 2 epidemiological indicators
 *
 * @author Luis Alvarez <lalvarez.mat@gmail.com>
 */


//#include <bits/stdc++.h>
#include "EpiIndicators.h"
#include "utilities.h"

//#ifndef R_COMPILE
//#define R_COMPILE
//#endif // R_COMPILE

/// ------------------------------------------------------------------------------------------------------
/// MEAN OF THE LAST ELEMENTS OF A VECTOR
/// ------------------------------------------------------------------------------------------------------
double mean_tail(
    vector<double> &i /** input vector */,
    int Ndays /** size of the tail to compute the mean */){
  double mean=0;
  for(int k=i.size()-Ndays;k<(int) i.size();k++) mean+=i[k];
  return mean/Ndays;
}

/// ------------------------------------------------------------------------------------------------------
/// GEOMETRIC MEAN OF THE LAST ELEMENTS OF A VECTOR
/// ------------------------------------------------------------------------------------------------------
double geometric_mean_tail(
    vector<double> &i /** input vector */,
    int Ndays /** size of the tail to compute the geometric mean */){
  double mean=1;
  for(int k=i.size()-Ndays;k<(int) i.size();k++) mean*=i[k];
  return pow(mean,1./Ndays);
}

/// ------------------------------------------------------------------------------------------------------
/// JOIN INDICATOR VALUES BY DATE
/// ------------------------------------------------------------------------------------------------------
int joint_indicators_by_date(
    vector<string> &date0,
    vector<double> &i0,
    vector<string> &date1,
    vector<double> &i1,
    vector<string> &date,
    vector<double> &f,
    vector<double> &g){
  if(date0.size()<2 || date0.size()!=i0.size() ) return -1;
  if(date0.size()<2 || date0.size()!=i0.size() ) return -1;
  
  time_t date_min=string2date(date0[0].c_str());
  time_t temp=string2date(date1[0].c_str());
  if(temp<date_min) date_min=temp;
  time_t date_max=string2date(date0[date0.size()-1].c_str());
  temp=string2date(date1[date1.size()-1].c_str());
  if(temp>date_max) date_max=temp;
  
  date.clear();
  f.clear();
  g.clear();
  for(time_t t=date_min+1000;t<date_max+2000;t+=86400){
    string current_date=date2string(t);
    date.push_back(current_date);
    f.push_back(0.);
    g.push_back(0.);
    for(int k=0;k< (int) date0.size();k++){
      if(current_date.compare(date0[k])==0){
        f[f.size()-1]=i0[k];
        break;
      }
    }
    for(int k=0;k< (int) date1.size();k++){
      if(current_date.compare(date1[k])==0){
        g[g.size()-1]=i1[k];
        break;
      }
    }
  }
  
  
  return 0;
}

/// ----------------------------------------------------------------------------------------------
/** CROUT METHOD TO SOLVE TRIANGULAR LINEAR SYSTEM. IT RETURNS THE VECTOR SOLUTION */
/// ----------------------------------------------------------------------------------------------
vector <double> crout(
    vector <double> &a /** diagonal values */,
    vector <double> &b /** upper-diagonal values */,
    vector <double> &c /** lower-diagonal values */,
    vector <double> &t /** independent system vector */)
{
  b.resize(a.size()-1);
  c.resize(a.size()-1);
  
  /// CROUT DECOMPOSITION
  vector <double> l(a.size());
  vector <double> m(c.size());
  vector <double> u(b.size());
  
  l[0]=a[0];
  if(l[0]==0) return vector<double>();
  u[0]=b[0]/l[0];
  for(int i=1;i< (int) (a.size()-1);i++){
    m[i-1]=c[i-1];
    l[i]=a[i]-m[i-1]*u[i-1];
    if(l[i]==0) return vector<double>();
    u[i]=b[i]/l[i];
  }
  m[m.size()-1]=c[c.size()-1];
  l[a.size()-1]=a[a.size()-1]-m[a.size()-2]*u[a.size()-2];
  
  /// SOLVING FIRST TRIANGULAR SYSTEM
  vector <double> z(l.size()); // declaramos array de salida
  if(l[0]==0) return(vector <double>()); // comprobamos division por cero
  z[0]=t[0]/l[0];
  for(int i=1;i< (int) z.size();i++){
    if(l[i]==0) return(vector <double>()); // comprobamos division por cero
    z[i]=(t[i]-z[i-1]*m[i-1])/l[i];
  }
  
  /// SOLVING SECOND TRIANGULAR SYSTEM
  vector <double> v(z.size()); // declaramos array de salida
  v[z.size()-1]=z[z.size()-1];
  for(int i=z.size()-2;i>=0;i--){
    v[i]=(z[i]-v[i+1]*u[i]);
  }
  
  return(v);
  
}

/// ----------------------------------------------------------------------------------------------
/// BASIC LINEAR INTERPOLATION. IT RETURNS THE INTERPOLATION VALUE (USING NEUMAN BOUNDARY CONDITION)
/// ----------------------------------------------------------------------------------------------
double linear_interpolation2(
    vector<double> &N /** data vector to interpolate */,
    double t /** time (index) location to interpolate the vector */ )
{
  if(t<=0.) return N[0];
  int t1=t;
  if(t1>=(int) N.size()-1) return N[N.size()-1];
  if(N[t1]==0 || N[t1+1]==0) return 0.;
  double dt=t-t1;
  return (1.-dt)*N[t1]+dt*N[t1+1];
}

/// ----------------------------------------------------------------------------------------------
/// APPLY VECTOR SHIFT TO A VECTOR DATA. IT RETURNS THE SHIFTED VECTOR
/// ----------------------------------------------------------------------------------------------
vector<double> apply_shift(
    vector <double> &v /** data vector */,
    vector<double> &shift /** shift vector */)
{
  vector<double> u(v.size(),0.);
  for(int k=0;k< (int) u.size();k++){
    double t=k+shift[k];
    u[k]=linear_interpolation2(v,t);
  }
  return u;
}

/// ----------------------------------------------------------------------------------------------
/// FIND LOCAL MAXIMA OF A VECTOR DATA. IT RETURNS THE INDEX OF THE LOCAL MAXIMA POSITIONS
/// ----------------------------------------------------------------------------------------------
vector<int> find_maxima(
    vector<double> &c /** vector with data */,
    int radius /** radius of the interval to check the local maxima */,
    double max_ratio_left /** minimum ratio difference for the data in the left extreme of the interval */,
    double max_ratio_right /** minimum ratio difference for the data in the right extreme of the interval */){
  
  vector<int> maxima;
  
  for(int k=radius;k< (int) c.size()-radius-1;k++){
    if(maxima.size()>0 && k-maxima[maxima.size()-1]<=radius) continue;
    bool max_pos=true;
    for(int m=k-radius;m<=k+radius;m++){
      if(c[m]>c[k]){
        max_pos=false;
        break;
      }
    }
    if(max_pos==false) continue;
    if( c[k]-c[k-radius]> c[k]*max_ratio_left && c[k]-c[k+radius]> c[k]*max_ratio_right){
      maxima.push_back(k);
    }
  }
  
  return maxima;
}

/// ----------------------------------------------------------------------------------------------
/// SHIFT INITIALIZATION BETWEEN 2 INDICATOR SEQUENCES. IT RETURNS THE SHIFT
/// ----------------------------------------------------------------------------------------------
vector<double> shift_initialization(
    vector<double> &c0 /** first epidemiological indicator sequence */,
    vector<double> &d0 /** second epidemiological indicator sequence */,
    int radius /** parameter when using the local maxima inicialization */,
    double max_ratio_left /** parameter when using the local maxima inicialization */,
    double max_ratio_right /** parameter when using the local maxima inicialization */,
    int s_min /** minimum value of the shift */,
    int s_max /** maximum value of the shift */,
    int type=0 /** initialization type. If type==0 we use the local maxima. Otherwise we divide the sequence
 in "type" subintervals and we look for local correspondences beetwen subintervals */)
{
  vector<int> x;
  vector<double> y;
  if(type==0){
    vector<double> c=c0;
    vector<double> d=d0;
    gauss_conv(c,2.,1);
    gauss_conv(d,2.,1);
    vector<int> max_pos_c=find_maxima(c,radius,max_ratio_left,max_ratio_right);
    vector<int> max_pos_d=find_maxima(d,radius,max_ratio_left,max_ratio_right);
    //    for(int k=0;k<max_pos_c.size();k++){
    //      printf("max_pos_c[%d]=%d c=%lf\n",k,max_pos_c[k],c[max_pos_c[k]]);
    //    }
    //    printf("\n");
    //    for(int k=0;k<max_pos_d.size();k++){
    //      printf("max_pos_d[%d]=%d d=%lf\n",k,max_pos_d[k],d[max_pos_d[k]]);
    //    }
    
    for(int k=0;k< (int) max_pos_c.size();k++){
      for(int m=0;m< (int) max_pos_d.size();m++){
        if( max_pos_d[m]>max_pos_c[k]+s_min && max_pos_d[m]<max_pos_c[k]+s_max){
          x.push_back(max_pos_c[k]);
          y.push_back((double) max_pos_d[m]-max_pos_c[k]);
          break;
        }
      }
    }
  }
  else if(type>0){
    x.resize(type);
    y.resize(type);
    double step = (double) c0.size()/(x.size()+1.);
    for(int k=0;k< (int) x.size();k++){
      x[k]=(k+1)*step;
      double max_cor=-1;
      for(double shift=s_min;shift<=s_max;shift+=0.1){
        int tmin=x[k]-step,tmin2=shift+x[k]-step;
        int tmax=x[k]+step,tmax2=shift+x[k]+step;
        if(tmin2>tmin) tmin=tmin2;
        if(tmin<0) tmin=0;
        if(tmax2<tmax) tmax=tmax2;
        if(tmax>= (int) c0.size()) tmax=c0.size()-1;
        vector<double> c(tmax-tmin+1);
        vector<double> d(tmax-tmin+1);
        for(int m=0;m< (int) c.size();m++){
          c[m]=c0[tmin+m];
          d[m]=linear_interpolation2(d0,tmin+m+shift);
        }
        double cm,cs,dm,ds,cor=0.;
        basic_statistics(c,cm,cs);
        basic_statistics(d,dm,ds);
        for(int m=0;m< (int) c.size();m++){
          cor+=(c[m]-cm)*(d[m]-dm);
        }
        cor=cor/(cs*ds*c.size());
        if(cor>max_cor){
          max_cor=cor;
          y[k]=shift;
        }
      }
    }
  }
  
  if(x.size()==0) return vector<double>(c0.size(),0.);
  if(x.size()==1) return vector<double>(c0.size(),y[0]);
  vector<double> s(c0.size(),y[0]);
  for(int k=x[ (int) x.size()-1];k< (int) c0.size();k++) s[k]=y[y.size()-1];
  for(int k=0;k< (int) x.size()-1;k++){
    for(int m=x[k];m<=x[k+1];m++){
      s[m]=y[k]+(m-x[k])*(y[k+1]-y[k])/(x[k+1]-x[k]);
    }
  }
  
  return s;
}



/// ------------------------------------------------------------------------------------------------------
/// ESTIMATION OF THE RATIO AND SHIFT BETWEEN INDICATORS FOR A GIVEN INITIALIZATION OF s(t). IT RETURNS
/// A NORMALIZED ERROR
/// ------------------------------------------------------------------------------------------------------
double shift_and_mortality_estimation(
    vector<double> &c0 /** first epidemiological indicator sequence */,
    vector<double> &d0 /** second epidemiological indicator sequence */,
    vector<double> &r2 /** output ratio between indicators */,
    vector<double> &s2 /** input/output shift between indicators */,
    int tmin /** left time interval extreme to compare the 2 indicators */,
    int tmax /** right time interval extreme to compare the 2 indicators */,
    double s_min /** minimum value of the shift */,
    double s_max /**  maximum value of the shift */,
    double wr /** regularity weight for the ratio */,
    double ws /** regularity weight for the shift */,
    int tail /** tail of values used to fit f(t)*r(t)-g(t+s(t)) for the last values of the sequence */,
    double tail_mu /** exponential coefficient to increase thw weight of the tail values */,
    double r_init /** optional value of the ratio on the left time interval extreme */,
    double r_end /** optional value of the ratio on the right time interval extreme */,
    double s_init /** optional value of the shift on the left time interval extreme */,
    double s_end  /** optional value of the shift on the right time interval extreme */)
{
  int Nc=tmax-tmin+1;
  vector<double>  c1(c0.size()),d1(d0.size()),dp(d0.size());
  vector<double>  c(Nc),d(Nc),s(Nc),r(Nc),p(Nc,1.);
  
  double c0_median=percentil( c0.size()/2,c0);
  double d0_median=percentil( d0.size()/2,d0);
  
  c1[0]=c0[0]/c0_median;
  d1[0]=d0[0]/d0_median;
  for(int k=1;k< (int) c1.size();k++){
    c1[k]=c0[k]/c0_median;
    d1[k]=d0[k]/d0_median;
    dp[k]=d1[k]-d1[k-1];
  }
  
  gauss_conv(dp,2.,1);
  
  for(int k=0;k<Nc;k++){
    c[k]=c1[tmin+k];
    s[k]=s2[tmin+k];
  }
  
  //int tail=56;
  //double mu=0.1;
  for(int k=Nc-tail;k<Nc;k++){
    p[k]=exp((k-Nc+tail+1)*tail_mu);
  }
  
  double Error0=0,Error1=0,Error2=0,Error0N=0,Error1N=0,Error2N=0,Error=0,ErrorN=0,dif=1e20;
  double Error3N=0;
  
  r2=vector<double>(c0.size(),0.);
  double minError=1e100,minError0;
  int iter=0;
  for(iter=0;iter<1000;iter++){
    //printf("%d ",iter);
    //system("pause");
    vector<double> dn(Nc),dpn(Nc);
    for(int k=0;k<Nc;k++){
      dn[k]=linear_interpolation2(d1,tmin+k+s[k]);
      dpn[k]=linear_interpolation2(dp,tmin+k+s[k]);
    }
    
    vector<double> diag(Nc,0.),sup(Nc,-wr),inf(Nc,-wr),b(Nc,0.);
    diag[0]=wr+c[0]*c[0];
    b[0]=c[0]*dn[0];
    diag[Nc-1]=wr+c[Nc-1]*c[Nc-1]*p[Nc-1];
    b[Nc-1]=c[Nc-1]*dn[Nc-1]*p[Nc-1];
    for(int k=1;k<Nc-1;k++){
      diag[k]=2*wr+c[k]*c[k]*p[k];
      b[k]=c[k]*dn[k]*p[k];
    }
    if(r_init!=-1e6){ diag[0]=1; sup[0]=0; b[0]=r_init/d0_median*c0_median; }
    if(r_end!=-1e6){ diag[Nc-1]=1; inf[Nc-2]=0; b[Nc-1]=r_end/d0_median*c0_median; }
    
    r=crout(diag,sup,inf,b);
    
    if(r.size()==0) return -1;
    
    if(iter==0){
      for(int k=1;k<Nc-1;k++){
        double x=r[k]*c[k]-linear_interpolation2(d1,tmin+k+s[k]);
        double y=k==0?0:r[k]-r[k-1];
        double z=k==0?0:s[k]-s[k-1];
        Error0+=p[k]*x*x;
        Error1+=wr*y*y;
        Error2+=ws*z*z;
      }
      Error=Error0+Error1+Error2;
    }
    
    
    /// CROUT s
    sup=vector<double>(Nc,-ws);
    inf=vector<double>(Nc,-ws);
    diag[0]=ws+dpn[0]*dpn[0];
    b[0]=(r[0]*c[0]-dn[0])*dpn[0]-ws*(s[0]-s[1]);
    diag[Nc-1]=ws+dpn[Nc-1]*dpn[Nc-1]*p[Nc-1];
    b[Nc-1]=(r[Nc-1]*c[Nc-1]-dn[Nc-1])*dpn[Nc-1]*p[Nc-1]-ws*(s[Nc-1]-s[Nc-2]);
    for(int k=1;k<Nc-1;k++){
      diag[k]=2*ws+dpn[k]*dpn[k]*p[k];
      b[k]=(r[k]*c[k]-dn[k])*dpn[k]*p[k]-ws*(2*s[k]-s[k-1]-s[k+1]);
    }
    
    for(int k=0;k<Nc;k++){
      if(s[k]>s_max){
        b[k]+=1e10*(s_max-s[k]);
        diag[k]+=1e10;
      }
      if(s[k]<s_min){
        b[k]+=1e10*(s_min-s[k]);
        diag[k]+=1e10;
      }
    }
    
    if(s_init!=-1e6){ diag[0]=1; sup[0]=0; b[0]=s_init-s[0]; }
    if(s_end!=-1e6){ diag[Nc-1]=1; inf[Nc-2]=0; b[Nc-1]=s_end-s[Nc-1]; }
    
    vector<double> sol=crout(diag,sup,inf,b);
    
    if(sol.size()==0) return -1;
    
    for(int k=0;k<Nc;k++) s[k]+=sol[k];//*dpn[k];
    
    Error0N=Error1N=Error2N=0;
    for(int k=1;k<Nc-1;k++){
      double x=r[k]*c[k]-linear_interpolation2(d1,tmin+k+s[k]);
      double y=k==0?0:r[k]-r[k-1];
      double z=k==0?0:s[k]-s[k-1];
      Error0N+=p[k]*x*x;
      Error1N+=wr*y*y;
      Error2N+=ws*z*z;
      if(s[k]>s_max) Error3N+=1e10*(s_max-s[k])*(s_max-s[k]);
      else if(s[k]<s_min) Error3N+=1e10*(s_min-s[k])*(s_min-s[k]);
    }
    ErrorN=Error0N+Error1N+Error2N+Error3N;
    dif=(Error-ErrorN)/ErrorN;
    
    if(ErrorN<minError){
      minError=ErrorN;
      minError0=Error0N;
      //printf("iter=%d minError=%e minError0=%lf\n",iter,minError,minError0); system("pause");
      for(int k=0;k< (int) r.size();k++) {
        r2[tmin+k]=r[k]*d0_median/c0_median;
        s2[tmin+k]=s[k];
      }
      for(int k=0;k<tmin;k++) {
        r2[k]=r2[tmin];
        s2[k]=s[0];
      }
      for(int k=r.size();k< (int) r2.size();k++) {
        r2[k]=r2[tmin+r.size()-1];
        s2[k]=s[s.size()-1];
      }
    }
    
    //printf("%1.0e %e\n",dif,ErrorN); system("pause");
    if(fabs(dif)<1e-6) break;
    Error0=Error0N;
    Error1=Error1N;
    Error2=Error2N;
    //Error3=Error3N;
    Error=ErrorN;
  }
  
  return minError0;//Error0N;
  
}

/// ------------------------------------------------------------------------------------------------------
/// FORECASTING THE INDICATOR f(t) FROM THE INDICATOR g(t)
/// ------------------------------------------------------------------------------------------------------
void indicators_forecast(
    vector<string> &dates /** date for the first indicator value (input/output) */,
    vector<double> &f /** first epidemiological indicator (input/output) */,
    vector<double> &g /** second epidemiological indicator with forecast (input) */,
    vector<double> &r0 /** ratio between indicators (input) */,
    vector<double> &s0 /** shift between indicators (input) */,
    double r_forecast_end /** final value of r(t) for the Hermite interpolation (input)*/,
    double s_forecast_end /** final value of s(t) for the Hermite interpolation (input) */,
    int NforecastDays /** number of forecast days for f(t) (input) */,
    int NdaysHermite /** Number of days used for the Hermite main component (input) */,
    int NdaysHermiteDerivada /** Number of days used for the Hermite derivative component (input) */
)
{
  vector<double> r=r0,s=s0;
  
  HermiteInterpolation(r,r_forecast_end,NdaysHermite,NdaysHermiteDerivada);
  HermiteInterpolation(s,s_forecast_end,NdaysHermite,NdaysHermiteDerivada);
  int Ndays=f.size()+NforecastDays;
  for(int k=f.size();k<Ndays;k++){
    double si=linear_interpolation2(s,(double) k);
    double ri=linear_interpolation2(r,(double) k);
    double gi=linear_interpolation2(g,(double) k+si);
    f.push_back(gi/ri);
    time_t temp=string2date(dates[dates.size()-1].c_str());
    temp+=86400+1000;
    dates.push_back(date2string(temp));
  }
  
  r0=r;
  s0=s;
}


/// ------------------------------------------------------------------------------------------------------
/// OPTIMIZATION OF THE RATIO AND SHIFT BETWEEN INDICATORS USING DIFFERENT INITIALIZATIONS OF s(t)
/// ------------------------------------------------------------------------------------------------------
double shift_and_ratio_optimization(
    vector<string> &dates /** date for each indicator value */,
    vector<double> &f /** first epidemiological indicator sequence */,
    vector<double> &g /** second epidemiological indicator sequence */,
    vector<double> &r /** output ratio between indicators */,
    vector<double> &s /** output shift between indicators */,
    int s_min /** minimum value of the shift */,
    int s_max /** maximum value of the shift */,
    double wr /** regularity weight for the ratio */,
    double ws /** regularity weight for the shift */,
    int tail /** tail of the dates used to fit f(t)*r(t)-g(t+s(t)) for the last values of the sequence */,
    double tail_mu /** exponential coefficient to increase the weight of the tail values */,
    double r_init /** optional value of the ratio on the left time interval extreme */,
    double r_end /** optional value of the ratio on the right time interval extreme */,
    double s_init /** optional value of the shift on the left time interval extreme */,
    double s_end /** optional value of the shift on the right time interval extreme */)
{
  while( (-s_min>0 && -s_min< (int) f.size() && f[(int) -s_min]==0)
           || (s_max>0 && s_max< (int) g.size() && g[s_max]==0) ){
    dates.erase(dates.begin());
    f.erase(f.begin());
    g.erase(g.begin());
  }
  //printf("tail=%d tail_mu=%lf\n",tail,tail_mu);
  int Nc=f.size()-1;
  while( (Nc-s_max>0 && Nc-s_max<  (int) f.size() && f[Nc-s_max]==0) ||
         (Nc+s_min>0 && Nc+s_min< (int) g.size() && g[Nc+s_min]==0) || g[Nc]==0){
    dates.resize(Nc);
    f.resize(Nc);
    g.resize(Nc--);
  }
  
  if(f.size()<56){
    return -1.;
  }
  
  Nc=f.size()/56;
  double min_error=2e20;
  double max_ratio_left=0.1;
  double max_ratio_right=0.1;
  int radius=28;
  for(int type=0;type<=Nc;type++){
    vector<double> r2,s2=shift_initialization(f,g,radius,max_ratio_left,max_ratio_right,s_min,s_max,type);
    double error=shift_and_mortality_estimation(f,g,r2,s2,0,f.size()-1,s_min,s_max,wr,ws,tail,tail_mu,r_init,r_end,s_init,s_end);
    if(error>=0 && error<min_error){
      min_error=error;
      r=r2;
      s=s2;
#ifndef R_COMPILE
      printf("type=%d min_error=%lf\n",type,min_error);
#endif
    }
  }
  
  return min_error<1e20?min_error:-1. ;
}
