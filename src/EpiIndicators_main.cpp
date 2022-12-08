#include <Rcpp.h>
using namespace Rcpp;
using namespace std; 

 
#include "EpiIndicators.h"


     
// [[Rcpp::export]]
DataFrame EpiIndicatorsC(
    CharacterVector date,
    NumericVector f,
    NumericVector g,
    int s_min=-10,
    int s_max=25,
    double wr=5000,
    double ws=50,
    double r_init=-1e6,
    double r_end=-1e6,
    double s_init=-1e6,
    double s_end=-1e6
)
{
  if(date.size()<10){
    stop("The number of samples : %d is too small\n",(int) date.size());
  }
  
  if(date.size()!=f.size() || f.size()!=g.size()){
    stop("the vector data have different size"); 
  }
  
  vector<double> c(f.size()),d(f.size());
  vector <string> dateC(f.size()); 
  for(int k=0;k<(int) c.size();k++){
    c[k]=f[k];
    d[k]=g[k];
    dateC[k]=string(date[k]); 
    if(strlen(dateC[k].c_str())!=10 || dateC[k].c_str()[4]!='-' || dateC[k].c_str()[7]!='-'){
      stop("EpiIndicators first argument must be dates in the format YYYY-MM-DD");
    }
  } 
  
  /// OUTPUTS
  vector<double> r;
  vector<double> s;
  
  
  double result=shift_and_ratio_optimization(dateC,c,d,r,s,s_min,s_max,wr,ws,r_init,r_end,s_init,s_end);
  
  if(result<0){
    stop("Unexpected problem using EpiIndicators()");
  }
  vector<double> ds=apply_shift(d,s),c_r(c.size());
  for(int k=0;k<(int) c.size();k++){
    c_r[k]=c[k]*r[k];
  }
  
  DataFrame results = DataFrame::create(
    Named("date") = dateC ,
    Named("f") = c ,
    Named("g") = d ,
    Named("r") = r ,
    Named("s") = s ,
    Named("f.r")  = c_r,
    Named("g.s")  = ds
  );
  return results;
}


// [[Rcpp::export]]
DataFrame joint_indicators_by_dateC(
    CharacterVector date0,
    NumericVector i0,
    CharacterVector date1,
    NumericVector i1)
{
  
  vector<double> i00(i0.size()),i11(i1.size());
  vector <string> date00(date0.size()),date11(date1.size()); 
  for(int k=0;k<(int) i0.size();k++){
    i00[k]=i0[k];
    date00[k]=string(date0[k]); 
    if(strlen(date00[k].c_str())!=10 || date00[k].c_str()[4]!='-' || date00[k].c_str()[7]!='-'){
      stop("EpiIndicators first argument must be dates in the format YYYY-MM-DD");
    }
  } 
  for(int k=0;k<(int) i1.size();k++){
    i11[k]=i1[k];
    date11[k]=string(date1[k]); 
    if(strlen(date11[k].c_str())!=10 || date11[k].c_str()[4]!='-' || date11[k].c_str()[7]!='-'){
      stop("EpiIndicators first argument must be dates in the format YYYY-MM-DD");
    }
  } 
  
  vector<string> date;
  vector<double> f,g; 
  
  int res=joint_indicators_by_date(date00,i00,date11,i11,date,f,g);
  if(res<0) stop("Problems joning the indicator tables");
  
  DataFrame results = DataFrame::create(
    Named("date") = date ,
    Named("f") = f ,
    Named("g") = g
  );
  return results;
  
}

// [[Rcpp::export]]
NumericVector apply_shiftC(
    NumericVector g,
    NumericVector s
){
  
  if(g.size()!=s.size()){
    stop("vectors g and s have different size"); 
  }
  
  vector<double> g2(g.size()),s2(g.size());
  for(int k=0;k<(int) g.size();k++){
    g2[k]=g[k];
    s2[k]=s[k]; 
  }
  vector<double> v=apply_shift(g2,s2); 
  NumericVector v2( v.begin(), v.end() );
  return v2; 
  
}


