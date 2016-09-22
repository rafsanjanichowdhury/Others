/*
  Muti Dimension integration using Monte Carlo method
  Integration of f(x1, x2, ... xn) on (a[i],b[i]) (i=1,n) intervals
  AG: February 2007
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double f(double[], int);
double int_mcnd(double(*)(double[],int), double[], double[], int, int); 

int main()
{
    const int n = 6;       /* define how many integrals */
//    const int m = 1000000; /* define how many points */
    double a[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /* left end-points */
    double b[n] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; /* right end-points */
    double result;
    int i, m;
    int ntimes;
    
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint); 
// current time in seconds (begin calculations)
    time_t seconds_i;
    seconds_i = time (NULL);

    m = 2;                // initial number of intervals
    ntimes = 20;          // number of interval doublings with nmax=2^ntimes
    cout << setw(12) << n <<"D Integral" << endl;
    for (i=0; i <=ntimes; i=i+1)
    {
        result = int_mcnd(f, a, b, n, m);
        cout << setw(10)  << m << setw(12) << result <<endl;
        m = m*2;
    }
    
//    cout << " " << n <<"D Integral = " << result << endl;

// current time in seconds (end of calculations)
    time_t seconds_f;
    seconds_f = time (NULL);
    cout << endl << "total elapsed time = " << seconds_f - seconds_i << " seconds" << endl << endl;
 

    system ("pause");
    return 0;
}

/*
 *  Function f(x1, x2, ... xk)
*/             
    double f(double x[], int n)
{
    double y;
    int j;
    y = 0.0;
/* a user may define the function explicitly like below */       
//    y = pow((x[0]+x[1]+x[2]+x[3]+x[4]),2);
/* or for a specific function like (x1 + x2 + ... xn)^2
   use a loop over n */
    for (j = 0; j < n; j = j+1)
      {
         y = y + x[j];
      }     
    y = pow(y,2);
     
    return y;
}

/*==============================================================
    Muti Dimension integration of f(x1,x2,...xn) 
    method: Monte-Carlo method
    written by: Alex Godunov (February 2007)
----------------------------------------------------------------
input:
    fn  - a multiple argument real function (supplied by the user)
    a[] - left end-points of the interval of integration
    b[] - right end-points of the interval of integration
    n  - dimension of integral
    m  - number of random points
output:
    r      - result of integration
Comments:
    be sure that following headers are included
    #include <cstdlib>
    #include <ctime>     
================================================================*/

    double int_mcnd(double(*fn)(double[],int),double a[], double b[], int n, int m)

{
    double r, x[n], p;
    int i, j;

    srand(time(NULL));  /* initial seed value (use system time) */

    r = 0.0;
    p = 1.0;

// step 1: calculate the common factor p
    for (j = 0; j < n; j = j+1)
      {
         p = p*(b[j]-a[j]);
      } 

// step 2: integration
    for (i = 1; i <= m; i=i+1)
    {
//      calculate random x[] points
        for (j = 0; j < n; j = j+1)
        {
            x[j] = a[j] + (b[j]-a[j])*rand()/(RAND_MAX+1);
        }         
        r = r + fn(x,n);
    }
    r = r*p/m;
    return r;
}
