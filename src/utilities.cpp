/**
 License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
 see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
 */


/**
 * @file utilities.cpp
 * @brief auxiliary function for EpiInvert method to compute the reproduction number and a restored incidence curve
 *
 * @author Luis Alvarez <lalvarez@ulpgc.es>
 */


#include "utilities.h"

/// LOG-NORMAL EVALUATION
double log_normal(
    double x,
    double mu,
    double sigma,
    double shift)
{
  x=x-shift;
  if(x<=0) return 0.;
  double aux=log(x)-mu;
  return exp(-aux*aux/(2*sigma))/( (x)*sqrt(2*M_PI*sigma) );
}

/// FUNCTION TO COMPUTE THE LINEAR REGRESSION (y=mx+n) BETWEEN VECTORS x AND y
/// IT RETURNS THE CORRELATION FACTOR
double linear_regression(
    vector<double> &x, /// INPUT VECTOR x
    vector<double> &y, /// INPUT VECTOR y
    double &m,double &n) /// PARAMETERS OF THE ESTIMATED LINEAR REGRESSION (y=mx+n)
{
  int N=x.size();
  if( ((int) y.size())!=N || N<2) return -1e10;

  /// AUXILIARY VARIABLES
  double suma_x=0,suma_y=0,suma_xy=0,suma_x2=0;
  for(int k=0;k<N;k++){
    suma_x+=x[k];
    suma_y+=y[k];
    suma_xy+=x[k]*y[k];
    suma_x2+=x[k]*x[k];
  }
  /// COMPUTATION OF m AND n
  double denominador=N*suma_x2-suma_x*suma_x;
  if(denominador==0.) return -1e10;
  m=(N*suma_xy-suma_x*suma_y)/denominador;
  n=(suma_x2*suma_y-suma_xy*suma_x)/denominador;

  /// CORRELATION FACTOR
  double mx,my,sx,sy,r=0;
  basic_statistics(x,mx,sx);
  basic_statistics(y,my,sy);

  for(int k=0;k<N;k++){
    r+=(x[k]-mx)*(y[k]-my);
  }
  r/=(x.size()*sx*sy);
  return r;


}

/// ------------------------------------------------------------------------------------------
/// ERROR MANAGEMENT
void fprintf_demo_failure(
    char mes[300] /** char array with the explanation about the execution failure */)
{
  FILE * f;
  f = fopen ("demo_failure.txt", "w");
  if(f==NULL){
    printf("We can not open demo_failure.txt");
    return;
  }
  fprintf(f,"%s\n",mes);
  fclose(f);
  return;
}

/// ------------------------------------------------------------------------------------------
/// BASIC STATISTICS: MEAN AND STANDARD DEVIATION (THE INPUT IS A SINGLE ARRAY)
void basic_statistics(
    vector<double> &i /** INPUT DATA IN VECTOR FORMAT*/,
    double &mean /** OUTPUT MEAN */,
    double &sd /** OUTPUT STANDARD DEVIATION */){

  if(i.size()==0) return;
  mean=sd=0;
  for(int k=0;k<(int) i.size();k++) mean+=i[k];
  mean/=i.size();

  for(int k=0;k<(int) i.size();k++) sd+=(i[k]-mean)*(i[k]-mean);
  sd=sqrt(sd/i.size());

}

/// ------------------------------------------------------------------------------------------
/// BASIC STATISTICS: MEAN AND STANDARD DEVIATION (THE INPUT IS A MATRIX)
void basic_statistics(
    vector < vector<double> > &i /** INPUT DATA IN MATRIX FORMAT*/,
    double &mean /** OUTPUT MEAN */,
    double &sd /** OUTPUT STANDARD DEVIATION */){

  mean=sd=0;
  int N=0;
  for(int k=0;k<(int) i.size();k++){
    for(int m=0;m<(int) i[k].size();m++){
      mean+=i[k][m];
      N++;
    }
  }
  if(N==0) return;
  mean/=N;

  for(int k=0;k<(int) i.size();k++){
    for(int m=0;m<(int) i[k].size();m++){
      sd+=(i[k][m]-mean)*(i[k][m]-mean);
    }
  }

  sd=sqrt(sd/N);

}

/// ------------------------------------------------------------------------------------------
/// PERCENTIL COMPUTATION
double percentil(
    const int k /** PERCENTIL MULTIPLIES BY THE VECTOR SIZE (k=0.5*x.size() CORRESPONDS TO THE MEDIAN)*/,
    const vector<double> &x /** INPUT VECTOR */)
{
  int i,ir,j,l,mid;
  double a,paso;
  vector<double> y=x;
  int n=x.size();

  l=0;
  ir=n-1;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && y[ir] < y[l]) {
        paso=y[l]; y[l]=y[ir]; y[ir]=paso;
      }
      a=y[k];
      return a;
    } else {
      mid=(l+ir) >> 1;
      paso=y[mid]; y[mid]=y[l+1]; y[l+1]=paso;
      if (y[l] > y[ir]) {
        paso=y[l]; y[l]=y[ir]; y[ir]=paso;
      }
      if (y[l+1] > y[ir]) {
        paso=y[l+1]; y[l+1]=y[ir]; y[ir]=paso;
      }
      if (y[l] > y[l+1]) {
        paso=y[l]; y[l]=y[l+1]; y[l+1]=paso;

      }
      i=l+1;
      j=ir;
      a=y[l+1];
      for (;;) {
        do i++; while (y[i] < a);
        do j--; while (y[j] > a);
        if (j < i) break;
        paso=y[i]; y[i]=y[j]; y[j]=paso;
      }
      y[l+1]=y[j];
      y[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}

/// ------------------------------------------------------------------------------------------
/// FUNCTION TO READ A DATABASE OF COVID-19 SEQUENCES (USED IN THE LEARNIG FORECAST PROCEDURE)
vector< vector<double> > read_matrix(
    const char name[] /** name of the file to read*/
)
{
  /// DATA SEQUENCE
  vector< vector<double> > cV0(1);

  /// READING THE DATA SEQUENCE
  FILE *f;
  f=fopen (name , "r");

  if(f==NULL) {
    printf("\nproblems reading file %s\n",name);
    //system("pause");
    exit(-1);
  }

  char c='0';
  char s[200];
  int k=0;
  int cont=0;
  while(! feof (f) ){
    c=getc(f);
    if(c==' '){
      //if(k==0) break;
      s[k]='\0';
      if(k>0) cV0[cont].push_back(atof(s));
      k=0;
      continue;
    }
    else if(c=='\n'){
      s[k]='\0';
      if(k>0) cV0[cont].push_back(atof(s));
      cont++;
      cV0.resize(cont+1);
      //break;
    }
    else {
      s[k++]=c;
    }
  }
  if(cV0[cV0.size()-1].size()==0) cV0.resize(cV0.size()-1);


  for(int m=0;m<(int) cV0.size();m++){
    if(cV0[m][cV0[m].size()-1]==0) cV0[m].resize(cV0[m].size()-1);
  }

  return cV0;
}

/// ------------------------------------------------------------------------------------------
/// ESTIMATION OF A LINEAR REGRESSION INTERPOLATING IN THE LAST 7 DAYS
vector<double> last_week_regression_interpolation(
    const vector<double> &d /** data to interpolate */)
{
  vector<double> P;
  int m=d.size()-4;
  P.push_back((d[m-3]+d[m+3]+d[m-2]+d[m+2]+d[m-1]+d[m+1]+d[m])/7.);
  P.push_back((3.*(d[m+3]-d[m-3])+2.*(d[m+2]-d[m-2])+d[m+1]-d[m-1])/28.);

  return P;
}

/// ------------------------------------------------------------------------------------------
/// EVALUATION OF THE INTERPOLATION POLYNOMIAL CENTERED IN THE LAST WEEK
double last_week_polynomial_evaluation(
    const int pos /** time to evaluate the polynomial */,
    const vector<double> &d /** data to interpolate */,
    const vector<double> &P /** polynomial coefficients */){
  if(P.size()==0) return -1e40;
  if(P.size()==1) return P[0];
  int m=d.size()-4;
  double x=pos-m;
  //printf("P.size()=%d x=%lf P[0]+P[1]*x=%lf\n",P.size(),x,P[0]+P[1]*x); system("pause");
  if(P.size()==2) return P[0]+P[1]*x;
  return P[0]+P[1]*x+P[2]*x*x;
}



/// ------------------------------------------------------------------------------------------
/// STANDARD GAUSS METHOD FOR SOLVING A LINEAR SYSTEM Au=b
vector< double > linear_system_solution(
    vector< vector<long double> > &A  /** MATRIX */,
    vector< long double > &b) /** INDEPENDENT VECTOR */
{

  vector< double >  u(b.size());

  for(int k=0;k<(int) b.size()-1;k++){
    long double max=fabs(A[k][k]);
    if(max==0) return vector<double>();

    for(int j=k+1;j<(int) b.size();j++){
      if(A[j][k]!=0){
        long double mul=-A[j][k]/A[k][k];
        A[j][k]=0.;
        for(int n=k+1;n<(int) b.size();n++) A[j][n]+=mul*A[k][n];
        b[j]+=mul*b[k];
      }
    }
  }

  if(A[b.size()-1][b.size()-1]==0.){
    if(b[b.size()-1]!=0) return vector<double>();
    else u[b.size()-1]=1.;
  }
  else{
    u[b.size()-1]=b[b.size()-1]/A[b.size()-1][b.size()-1];
  }
  for(int k=b.size()-2;k>=0;k--)
  {
    long double sum=0;
    for(int l=k+1;l<(int) b.size();l++) sum+= A[k][l]*u[l];
    u[k]=(b[k]-sum)/A[k][k];
  }

  return(u);
}

/// ------------------------------------------------------------------------------------------
/// FUNCTION TO READ COUNTRY DATA FROM THE FILE "owid-covid-data.csv"
/// RETURN  A VECTOR WITH THE NEW CASES AND THE DATE OF THE FINAL AVAILABLE DATA
vector<double>  read_country(
    const char C[] /** INPUT. ACRONIM OF THE COUNTRY IN THE FILE owid-covid-data.csv**/,
                char date[] /** OUTPUT. LAST DATE OF AVAILABLE DATA */){
  vector<double> i;
  //printf("%s\n",C);
  FILE *f;
  f=fopen ("owid-covid-data.csv", "r");

  if(f==NULL){
    char mes[300];
    sprintf(mes,"%s not found. Maybe the data file owid-covid-data.csv is not in this repository\n or there is an error in the country acronim or your own file data is not available",C);
    printf("\n%s\n\n",mes);
    fprintf_demo_failure(mes);
    //system("pause");
    exit(0);
  }
  char c='0';
  char s[200];
  int k=0;
  c=getc(f);
  while(! feof (f)){
    /// FIND THE FIRST STRING OF THE LINE
    if(c!=','){
      s[k++]=c;
      c=getc(f);
      continue;
    }
    s[k]='\0';

    /// IF THE FIRST STRING IS NOT THE SAME THAT THE SPECIFIED COUNTRY, WE CONTINUE WITH THE NEXT LINE
    if( strcmp(s,C)!=0){
      while(c!='\n' && !feof (f)){
        if(c==EOF) return i;
        c=getc(f);
      }
      if(c==EOF) return i;
      c=getc(f);
      k=0;
      continue;
    }
    /// LOOKING FOR THE DATE POSITION TO UPDATE THE LAST AVALAIBLE DATE
    c=getc(f);
    while(c!=',')  c=getc(f);
    c=getc(f);
    while(c!=',')  c=getc(f);
    c=getc(f);
    k=0;
    while(c!=','){
      date[k++]=c;
      c=getc(f);
    }
    date[k]='\0';


    /// LOOKING FOR THE NEW DAILY CASE POSITION
    c=getc(f);
    while(c!=',')  c=getc(f);
    c=getc(f);
    k=0;
    while(c!=','){
      s[k++]=c;
      c=getc(f);
    }
    s[k]='\0';

    /// ADDING THE NEW DAILY VALUE TO THE VECTOR
    //printf("%s\n",s);
    if(k>0) i.push_back(atof(s));
    else i.push_back(0.);
    k=0;

  }
  printf("last date : %s, last incidence : %1.0lf\n",date,i[i.size()-1]);

  //system("pause");
  if(i.size()==0){
    char mes[300];
    sprintf(mes,"%s not found. It can be a problem with the country acronym\n (if you are taking the data from owid-covid-data.csv) \n or the file does not exist",C);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
  }
  return i;
}


/// --------------------------------------------------------------------------------------------------------------
/// FUNCTION TO READ THE INCIDENCE DATA. THE DATA ARE DETERMINED BY THE CHAR ARRAY name[] AND CAN BE ORGANIZED AS :
///  (1) THE ACRONIM OF THE COUNTRY IN THE FILE owid-covid-data.csv
///  (2) THE NAME OF A FILE WITH THE INCIDENCE DATA OF A SINGLE REGION WITH NO DATE INFORMATION
///  (3) THE NAME OF A FILE WITH THE INCIDENCE DATA OF A SEVERAL REGIONS WITH NO DATE INFORMATION (EACH ROW REPRESENTS A REGION). THE DATA ARE
///      SEPARATED BY BLANK SPACES. EpiInvert USES THE CUMULATIVE VALUE OF THE INCIDENCE FOR ALL REGIONS.
vector< vector<double> > read_data_multiple(
    const char name[] /** FILE NAME OR ACRONIM OF THE COUNTRY IN THE FILE owid-covid-data.csv*/,
                   time_t &current_day /** OUTPUT LAST DATE OF AVAILABLE DATA */)
{

  vector< vector<double> > cV /** OUTPUT MATRIX WITH THE DATA INFORMATION*/;

  /// WE TRY TO OPEN THE FILE WITH THE DATA
  FILE *f;
  f=fopen (name , "r");
  /// IF THE FILE DOES NOT EXIST WE TRY TO GET THE DATA FROM THE FILE owid-covid-data.csv
  if(f==NULL){
    char date[300];
    vector<double> c=read_country(name,date);
    current_day = string2date(date);
    //printf("%s %d \n",date,(int)current_time ); system("pause");
    /// WE CHECK THE SIZE OF THE INCIDENCE DATA
    if(c.size()>0){
      cV.push_back(c);
      return cV;
    }
    char mes[300];
    printf("Problems reading %s\n",name);
    sprintf(mes,"Problems reading %s\n",name);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
    //return vector< vector<double> >();
  }

  /// WE CHECK IF THE FILE CONTAINS DATE INFORMATION (WE LOOK FOR THE CHARACTER ';' USED TO SEPARATE DATA FROM DATES)
  bool data_with_dates=false;
  char c=getc(f);
  while ( ! feof (f) ){
    if(c==';'){
      data_with_dates=true;
      break;
    }
    c=getc(f);
  }
  fclose(f);

  f=fopen (name , "r");


  /// WE READ THE DATA IF IT INCLUDES DATE INFORMATION OR THE FILE DOES NOT EXIST
  if(data_with_dates==true){
    cV.push_back(read_data_single_and_date(name,current_day));
    return cV;
  }


  time(&current_day);
  current_day-=3600*24;
  current_day=((int) current_day/86400)*86400;


  /// WE READ THE DATA FROM THE FILE
  c='0';
  char s[200];
  int k=0;
  while(true){
    c=getc(f);
    if(c==' '){
      if(k==0) break;
      s[k]='\0';

      cV.push_back(vector<double>(1,atof(s)));
      k=0;
      continue;
    }
    else if(c=='\n'){
      s[k]='\0';
      cV.push_back(vector<double>(1,atof(s)));
      break;
    }
    else {
      s[k++]=c;
      //      if(c!='0' && c!='1' && c!='2' && c!='3' && c!='4' && c!='5' && c!='6' && c!='7' && c!='8' && c!='9' && c!='.' && c!='-'){
      //        return vector< vector<double> >();
      //      }
    }
  }

  k=0;
  while ( ! feof (f) ){
    double x;
    fscanf(f,"%lf\n",&x);
    //printf("%1.0lf ",x); system("pause");
    cV[k%cV.size()].push_back(x);
    k++;
  }
  fclose(f);

  if(cV.size()==0){
    char mes[300];
    printf("Problems reading %s\n",name);
    sprintf(mes,"Problems reading the data in %s\n",name);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
  }

  if(cV[0].size()<20){
    char mes[300];
    sprintf(mes,"Sorry lacunary data. The country  has only %d samples\n",(int) cV[0].size());
    fprintf_demo_failure(mes);
    exit(0);
  }

  for(int k=1;k<(int) cV.size();k++){
    if(cV[0].size()!=cV[k].size()){
      cV[k].push_back(0.);
      if(cV[0].size()!=cV[k].size()){
        char mes[300];
        printf("Problems with the data file in %s\n",name);
        sprintf(mes,"Problems with the data file in %s\n",name);
        printf("%s\n",mes);
        fprintf_demo_failure(mes);
        exit(0);
      }
    }
  }
  return cV;
}

///----------------------------------------------------------------------------------------------------
/// BASIC LINEAR INTERPOLATION
double linear_interpolation(vector<double> &N,double t){
  if(t<=0.) return N[0];
  int t1=t;
  if(t1>=(int) N.size()-1) return N[N.size()-1];
  double dt=t-t1;
  return (1.-dt)*N[t1]+dt*N[t1+1];
}

///----------------------------------------------------------------------------------------------------
/// L2 ERROR (USED TO COMPUTE THE BEST MATCHING BETWEEN EPIESTIM AND EPIINVERT)
double L2(
    vector<double> &c1 /** first incidence data */,
    double t1 /** time to shift the first incidence data */,
    vector<double> &c2 /** second incidence data*/,
    double t2 /** time to shift the second incidence data */,
    int kmin /** min time to evaluate the L2 error */,
    int kmax /** max time to evaluate the L2 error */)
{
  double sum=0;
  for(int k=kmin;k<=kmax;k++){
    double aux1=linear_interpolation(c1,t1+k);
    double aux2=linear_interpolation(c2,t2+k);
    sum+=(aux1-aux2)*(aux1-aux2);
  }
  return sqrt(sum/(kmax-kmin+1));
}


/// ------------------------------------------------------------------------------------------
/// FUNCTION TO READ THE DATA IN A SINGLE SIMPLE FILE.
void read_data_single(char name[],vector<double> &c){
  FILE *f;
  f=fopen (name , "r");
  c.clear();
  if(f==NULL){
    return;
  }

  while ( ! feof (f) ){
    double x;
    fscanf(f,"%lf\n",&x);
    c.push_back(x);
    //printf("%1.0lf\n",c[c.size()-1]);
  }
  fclose(f);
}

/// ------------------------------------------------------------------------------------------
/// FUNCTION TO COMPUTE THE PERCENTIL IN THE LAST DAYS
vector<double> back_percentil(vector<double> &i,int radius){
  vector<double> nf(i.size());
  for(int k=0;k<(int) i.size();k++){
    vector<double> v;
    for(int n=k;(int) v.size()<radius && n>=0;n--){  v.push_back(i[n]);  }
    if((int) v.size()<radius){
      for(int n=k+1;(int) v.size()<radius && n<(int) i.size();n++){  v.push_back(i[n]);  }
    }
    nf[k]=fabs(percentil(0.5*v.size(),v))+1.;
  }
  return nf;
}

///--------------------------------------------------------------------------------
/// FUNCTION TO CONVERT A CHAR ARRAY WITH THE FORMAT "YYYY-MM-DD" TO time_t
time_t string2date(const char *date){

  struct tm t;
  t.tm_year=0;
  t.tm_mon=0;
  t.tm_mday=0;
  t.tm_sec=0;
  t.tm_min=0;
  t.tm_hour=0;
  t.tm_isdst=-1;

  time_t time=mktime(&t);


  /// YEAR
  char c[10];
  int k=0,n=0;
  while(date[k]!='/' && date[k]!='-' && k<(int) strlen(date)){
    c[n++]=date[k++];
  }
  c[n]='\0';
  //printf("year=%s\n",c); //system("pause");
  t.tm_year=atoi(c)-1900;
  if(k==(int) strlen(date) || t.tm_year<0 || t.tm_year>200){
    return time;
  }

  /// MONTH
  k++;
  n=0;
  while(date[k]!='\\' && date[k]!='-' && k<(int) strlen(date)){
    c[n++]=date[k++];
  }
  c[n]='\0';
  //printf("month=%s\n",c);
  t.tm_mon=atoi(c)-1;
  if(k==(int) strlen(date) || t.tm_mon<0 || t.tm_mon>11){
    return time;
  }

  /// DAY
  n=0;
  k++;
  while(date[k]!='\\' && date[k]!='-' && k<(int) strlen(date)){
    c[n++]=date[k++];
  }
  c[n]='\0';
  //printf("day=%s\n",c);
  t.tm_mday=atoi(c);
  if(t.tm_mday<0 || t.tm_mday>31){
    return time;
  }

  return mktime(&t);

}

///--------------------------------------------------------------------------------
/// FUNCTION TO READ DATA CONTAINING THE DATE AND THE INCIDENCE VALUE
vector<double>  read_data_single_and_date(
    const char filename[] /** INPUT. DATA FILENAME **/,
                       time_t &current_time /** OUTPUT. LAST DATE OF AVAILABLE DATA */)
{
  char date[200];
  vector<double> i;
  FILE *f;
  f=fopen (filename, "r");
  if(f==NULL){
    char mes[300];
    printf("Problems reading %s\n",filename);
    sprintf(mes,"Problems reading %s\n",filename);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
  }

  char c='0';
  char s[200];
  int k=0;
  while(! feof (f)){
    c=getc(f);
    /// FIND THE FIRST STRING OF THE LINE
    while(c!=';' && c!='\n' && ! feof (f)){
      s[k++]=c;
      c=getc(f);
      continue;
    }
    s[k]='\0';
    //printf("%s\n",s);
    if(c==';') strcpy(date,s);
    else if(c=='\n') i.push_back(atof(s));
    k=0;
  }
  i.push_back(atof(s));
  while(i[i.size()-1]<=0 && i.size()>0) i.resize(i.size()-1);
  //printf("%s\n",date);
  current_time=string2date(date);
  //system("pause");
  return i;
}

///--------------------------------------------------------------------------------
/// FUNCTION TO READ THE FILE WITH THE FESTIVE DAYS
void read_festive_days(
    char name[] /** INPUT FILE NAME WITH THE FESTIVE DAYS IN THE FORMAT YYYY-MM-DD */,
             vector<string> &festive_dates /** OUTPUT VECTOR WITH THE STRING OF THE FESTIVE DAYS*/)
{
  festive_dates.clear();
  FILE *f;
  f=fopen (name , "r");
  if(f==NULL){
    printf("Problems reading file %s\n",name);
    return;
  }

  while ( ! feof (f) ){
    char date[200];
    fscanf(f,"%s\n",date);
    //printf("%s\n",date);
    string s(date);
    if(s.length()>7) festive_dates.push_back(s);
  }

  //for(int k=0;k<festive_dates.size();k++){  printf("%s\n",festive_dates[k].c_str());   }
  fclose(f);
}

///-------------------------------------------------------------------------------------
/// FUNCTION TO INITIALIZE THE VECTOR TO DETERMINE THE FESTIVE DAYS. IT RETURNS A VECTOR
/// WITH THE VALUE 1 FOR A FESTIVE DAY AN 0 OTHERWISE
vector<int> daily_festive_day_initialization(
    time_t current_time /** TIME OF THE LAST AVAILABLE DATE */,
    int i_size /** SIZE OF THE INCIDENCE DATA (THE LASTA POSITION CORRESPONDS TO THE CURRENT TIME */,
                                               vector<string> &festive_dates /** INPUT VECTOR WITH THE STRING OF THE FESTIVE DAYS*/)
{
  //printf("festive_dates.size()=%d, current_time=%d, i_size=%d\n",festive_dates.size(),(int) current_time,i_size);
  if(festive_dates.size()==0 || current_time<=0 || i_size==0) return vector<int> ();
  vector<int> h(i_size,0);

  for(int k=0;k<(int) festive_dates.size();k++){
    time_t t2=current_time-string2date(festive_dates[k].c_str());
    if(t2<0) continue;
    int n = i_size-1-round(t2/86400.);
    //printf("%s n=%d  i_size=%d\n",festive_dates[k],n,i_size);
    if(n>=0 && n<i_size) h[n]=1;
  }
  //  for(int k=0;k<h.size();k++) printf("%d\n",h[k]);
  //  system("pause");
  return h;

}

/// FUNCTION TO GET ONE OF THE STORED FESTIVE DAY SEQUENCES
vector<string> get_stored_festive_days(vector<double> &i){
  vector<string> fd;

  /// WE COMPUTE THE MEAN AND SD OF THE DALY INCIDENCE THE FIRST YEAR
  if(i.size()<365) return fd;
  vector<double> i2=i;
  i2.resize(365);
  double mean,sd;
  basic_statistics(i2,mean,sd);
  //printf("mean=%lf, sd=%lf\n",mean,sd);


  /// DETECTING USA
  if( fabs(mean-67356)<1000 && fabs(sd-72005)<1000){
    printf("-> Festive days : USA detected\n"); //system("pause");
    fd.push_back(string("2020-01-01"));
    fd.push_back(string("2020-01-20"));
    fd.push_back(string("2020-02-17"));
    fd.push_back(string("2020-05-25"));
    fd.push_back(string("2020-06-21"));
    fd.push_back(string("2020-07-03"));
    fd.push_back(string("2020-09-07"));
    fd.push_back(string("2020-10-12"));
    fd.push_back(string("2020-11-11"));
    fd.push_back(string("2020-11-26"));
    fd.push_back(string("2020-12-25"));
    fd.push_back(string("2020-12-31"));
    fd.push_back(string("2021-01-01"));
    fd.push_back(string("2021-01-18"));
    fd.push_back(string("2021-02-15"));
    fd.push_back(string("2021-05-31"));
    fd.push_back(string("2021-06-18"));
    fd.push_back(string("2021-07-05"));
    fd.push_back(string("2021-09-06"));
    fd.push_back(string("2021-10-11"));
    fd.push_back(string("2021-11-11"));
    fd.push_back(string("2021-11-25"));
    fd.push_back(string("2021-11-28"));
    fd.push_back(string("2021-12-24"));
    fd.push_back(string("2021-12-31"));
    fd.push_back(string("2022-01-01"));
    fd.push_back(string("2022-01-17"));
    fd.push_back(string("2022-02-21"));
    fd.push_back(string("2022-05-30"));
    fd.push_back(string("2022-06-18"));
    fd.push_back(string("2022-07-04"));
    fd.push_back(string("2022-09-05"));
    fd.push_back(string("2022-10-10"));
    fd.push_back(string("2022-11-11"));
    fd.push_back(string("2022-11-24"));
    fd.push_back(string("2022-12-26"));
    fd.push_back(string("2022-12-31"));
    return fd;
  }
  /// DETECTING FRANCE
  if( (fabs(mean-8250)<100 && fabs(sd-12188)<100) ||  (fabs(mean-8486)<100 && fabs(sd-13664)<100) ){
    printf("-> Festive days : FRANCE detected\n"); //system("pause");
    fd.push_back(string("2020-01-01"));
    fd.push_back(string("2020-04-10"));
    fd.push_back(string("2020-04-13"));
    fd.push_back(string("2020-05-01"));
    fd.push_back(string("2020-05-08"));
    fd.push_back(string("2020-05-21"));
    fd.push_back(string("2020-06-01"));
    fd.push_back(string("2020-07-14"));
    fd.push_back(string("2020-08-15"));
    fd.push_back(string("2020-11-01"));
    fd.push_back(string("2020-11-11"));
    fd.push_back(string("2020-12-25"));

    fd.push_back(string("2021-01-01"));
    fd.push_back(string("2021-04-02"));
    fd.push_back(string("2021-04-05"));
    fd.push_back(string("2021-05-01"));
    fd.push_back(string("2021-05-08"));
    fd.push_back(string("2021-05-13"));
    fd.push_back(string("2021-05-24"));
    fd.push_back(string("2021-07-14"));
    fd.push_back(string("2021-08-15"));
    fd.push_back(string("2021-11-01"));
    fd.push_back(string("2021-11-11"));
    fd.push_back(string("2021-12-25"));

    fd.push_back(string("2022-01-01"));
    fd.push_back(string("2022-04-15"));
    fd.push_back(string("2022-04-18"));
    fd.push_back(string("2022-05-01"));
    fd.push_back(string("2022-05-08"));
    fd.push_back(string("2022-05-26"));
    fd.push_back(string("2022-06-06"));
    fd.push_back(string("2022-07-14"));
    fd.push_back(string("2022-08-15"));
    fd.push_back(string("2022-11-01"));
    fd.push_back(string("2022-11-11"));
    fd.push_back(string("2022-12-25"));
    return fd;
  }
  /// DETECTING GERMANY
  if( (fabs(mean-5885)<100 && fabs(sd-8199)<100) || (fabs(mean-5903)<100 && fabs(sd-8836)<100) ){
    printf("-> Festive days : GERMANY detected\n"); //system("pause");
    fd.push_back(string("2020-01-01"));
    fd.push_back(string("2020-04-10"));
    fd.push_back(string("2020-04-13"));
    fd.push_back(string("2020-05-01"));
    fd.push_back(string("2020-05-21"));
    fd.push_back(string("2020-06-01"));
    fd.push_back(string("2020-10-03"));
    fd.push_back(string("2020-12-25"));

    fd.push_back(string("2021-01-01"));
    fd.push_back(string("2021-04-02"));
    fd.push_back(string("2021-04-05"));
    fd.push_back(string("2021-05-01"));
    fd.push_back(string("2021-05-13"));
    fd.push_back(string("2021-05-24"));
    fd.push_back(string("2021-10-03"));
    fd.push_back(string("2021-12-25"));

    fd.push_back(string("2022-01-01"));
    fd.push_back(string("2022-04-15"));
    fd.push_back(string("2022-04-18"));
    fd.push_back(string("2022-05-01"));
    fd.push_back(string("2022-05-26"));
    fd.push_back(string("2022-06-06"));
    fd.push_back(string("2022-10-03"));
    fd.push_back(string("2022-12-25"));
    return fd;
  }
  /// DETECTING SPAIN
  if( (fabs(mean-13032)<100 && fabs(sd-22625)<100) || (fabs(mean-7515)<100 && fabs(sd-12590)<100)  ){
    printf("-> Festive days : SPAIN detected\n"); //system("pause");
    fd.push_back(string("2020-01-01"));
    fd.push_back(string("2020-01-06"));
    fd.push_back(string("2020-04-10"));
    fd.push_back(string("2020-04-13"));
    fd.push_back(string("2020-05-01"));
    fd.push_back(string("2020-08-15"));
    fd.push_back(string("2020-10-12"));
    fd.push_back(string("2020-11-01"));
    fd.push_back(string("2020-12-06"));
    fd.push_back(string("2020-12-08"));
    fd.push_back(string("2020-12-25"));

    fd.push_back(string("2021-01-01"));
    fd.push_back(string("2021-01-06"));
    fd.push_back(string("2021-04-01"));
    fd.push_back(string("2021-05-01"));
    fd.push_back(string("2021-08-13"));
    fd.push_back(string("2021-10-12"));
    fd.push_back(string("2021-11-01"));
    fd.push_back(string("2021-12-06"));
    fd.push_back(string("2021-12-08"));
    fd.push_back(string("2021-12-25"));

    fd.push_back(string("2022-01-01"));
    fd.push_back(string("2022-01-06"));
    fd.push_back(string("2022-04-15"));
    fd.push_back(string("2021-05-01"));
    fd.push_back(string("2022-08-15"));
    fd.push_back(string("2022-10-12"));
    fd.push_back(string("2022-11-01"));
    fd.push_back(string("2022-12-06"));
    fd.push_back(string("2022-12-08"));
    fd.push_back(string("2022-12-25"));

    return fd;
  }


  return fd;

}


double evaluation_init_extrapolation_14(int t,vector <double> &C)
{
  /// LINEAR REGRESSION
  if(C.size()==2){
    return C[0]+C[1]*t;
  }
  /// EXPONENTIAL APPROXIMATION
  else if(C.size()==3){
    return C[0]*exp(C[1]*t)+C[2];
  }
  return 1e20;
}

double linear_regression_14(
    vector<double> &i /** daily infected */,
    vector <double> &C /** linear regression interpolation coefficients */)
{
  if(i.size()<14) return -1.;

  double suma_x=0,suma_y=0,suma_xy=0,suma_x2=0;
  int N=14;
  for(int k=0;k<N;k++){
    suma_x+=k;
    suma_y+=i[k];
    suma_xy+=k*i[k];
    suma_x2+=k*k;
  }

  double denominador=N*suma_x2-suma_x*suma_x;
  if(denominador==0.) return -1.;

  C.clear();
  C.push_back((suma_x2*suma_y-suma_xy*suma_x)/denominador);
  C.push_back((N*suma_xy-suma_x*suma_y)/denominador);

  //printf("C[0]=%lf C[1]=%lf\n",C[0],C[1]);

  double error=0;
  for(int k=0;k<14;k++){
    double y = evaluation_init_extrapolation_14(k,C);
    //printf("i[%d]=%lf y=%lf\n",k,i[k],y);
    error+=(y-i[k])*(y-i[k]);
  }
  //printf("error=%e\n",error);
  return sqrt(error/14.);

}

double exponential_approximation_14(
    vector<double> &i /** daily infected */,
    vector <double> &C /** exponential approximation interpolation coefficients */)
{
  if(i.size()<14) return -1.;

  double a,b,suma_x=0,suma_y=0,suma_xy=0,suma_x2=0;
  int N=14;
  for(int k=0;k<N;k++){
    suma_x+=k;
    double y=i[k]>0?log(i[k]):0;
    suma_y+=y;
    suma_xy+=k*y;
    suma_x2+=k*k;
  }

  double denominador=N*suma_x2-suma_x*suma_x;
  if(denominador==0.) return -1.;

  C.clear();
  a=(N*suma_xy-suma_x*suma_y)/denominador;
  b=(suma_x2*suma_y-suma_xy*suma_x)/denominador;

  C.push_back(exp(b));
  C.push_back(a);


  double c=0;
  for(int k=0;k<14;k++){
    c+=i[k]-C[0]*exp(C[1]*k);
  }
  C.push_back(c/14.);

  //printf("C[0]=%lf C[1]=%lf C[2]=%lf\n",C[0],C[1],C[2]);


  double error=0;
  for(int k=0;k<14;k++){
    double y = evaluation_init_extrapolation_14(k,C);
    //printf("i[%d]=%lf y=%lf\n",k,i[k],y);
    error+=(y-i[k])*(y-i[k]);
  }
  //printf("error=%e\n",error);

  return sqrt(error/14.);

}

///----------------------------------------------------------------------------------------------------
/// BUILD A PARAMETRIC SERIAL INTERVAL FROM A SHIFTED LOG-NORMAL
int parametric_si_distr(
    double mean,
    double sd,
    double shift,
    vector<double> &si_distr /** output vector with the serial interval*/)
{

  int k0=round(shift);
  double variance = sd*sd;

  double sigma=log(variance/(mean*mean)+1);
  double mu=log(mean)-sigma/2.;
  //printf("mu=%lf sigma=%lf\n",mu,sigma);
  si_distr.clear();
  double aux1,aux2=0,sum=0;
  for(double x=k0;x<1000;x++){
    aux1=aux2;
    aux2=(4.*log_normal(x,mu,sigma,shift)+log_normal(x-0.5,mu,sigma,shift)+log_normal(x+0.5,mu,sigma,shift))/6.;
    if(aux2<aux1 && aux2<1e-5) break;
    if(aux2<1e-6){
      k0++;
      continue;
    }
    //printf("si[%d]=%lf\n",si_distr.size(),aux2); system("pause");
    si_distr.push_back(aux2);
    sum+=aux2;
  }

  for(int k=0;k< (int) si_distr.size();k++){
    si_distr[k]/=sum;
    //printf("%lf\n",si_distr[k]);
  }
  //printf("k0=%d\n",k0);
  //system("pause");
  return k0;

}

///----------------------------------------------------------------------------------------------------
/// READING THE SERIAL INTERVAL FROM A FILE. IT RETURNS f0 (THE POSITION OF ZERO IN THE SERIAL INTERVAL)
int read_si_distr(
    const char name[] /** file name */,
                   vector<double> &si_distr /** output vector with the serial interval*/)
{
  FILE *f;
  int k1=0,f0=-100;
  f=fopen (name , "r");
  if(f==NULL){
    printf("Problems reading the serial interval file : %s\n",name);
    char mes[300];
    sprintf(mes,"Problems reading the serial interval file : %s\n",name);
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
  }
  si_distr.clear();
  while ( ! feof (f) ){
    double x;
    int k;
    fscanf(f,"%d  %lf\n",&k,&x);
    if(f0==-100){
      f0=-k;
      k1=k;
    }
    else{
      if(k!=k1+1) break;
      k1=k;
    }
    si_distr.push_back(x);
  }
  fclose(f);

  if(si_distr.size()<5){
    printf("The size of the serial interval is : %d (too small)\n",(int) si_distr.size());
    char mes[300];
    sprintf(mes,"The size of the serial interval is : %d (too small) \n",(int) si_distr.size());
    printf("%s\n",mes);
    fprintf_demo_failure(mes);
    exit(0);
  }

  /// NORMALIZATION OF THE SERIAL INTERVAL
  double sum=0;
  for(int k=0;k<(int) si_distr.size();k++) sum+=si_distr[k];
  //printf("sum=%lf\n",sum);
  for(int k=0;k<(int) si_distr.size();k++) si_distr[k]/=sum;

  return f0;
}
