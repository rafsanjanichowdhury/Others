/*
 Solver for a system of n first order Ordinary Differential Equations 
 (initial value problem)
 Method: calls 4th-order Runge-Kutta
 Can be used for solving a system of n/2 second order ODE

 Demo version for 2D projectile motion
 AG: Last revision October 2009
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

/* function prototypes */
double dnx(double, double [], double [], int );
double rk4_dn1(double(*)(double, double [], double [], int ),
               double, double, double [], double [], int);

/* global variables */
const double g = 9.81;               // free fall acceleration in m/s^2
const double m = 1.0;                // mass of a projectile in kg
const double rad = 3.1415926/180.0;  // radians

int main()
{
    const int n=4;                   // number of first-order equations 
    double ti, tf, dt, tmax;
    double xi[n], xf[n];
    double v0, a0;
    int i, key;


/* output: file and formats */
    ofstream file;
    file.open ("table01c.dat");         // write results to this file
/* output format */
    file.precision(3);
    file.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);

/* initial information */
    ti = 0.0;                // initial value for variable t
    v0 = 180.0;              // initial speed (m/s)
    a0 =  45.0;              // initial angle (degrees)
    xi[0] = 0.0;             // initial position in x (m)
    xi[1] = 0.0;             // initial position in y (m)
    xi[2] = v0*cos(a0*rad);  // initial speed in x direction (m.s)
    xi[3] = v0*sin(a0*rad);  // initial speed in y direction (m/s)
    
    dt = 0.2;             // step size for integration (s)
    tmax = 60.0;          // integrate till tmax (s)
/* end of initial information */

    file << setw(12) << "t"  << setw(12) << "x"   << setw(12) << "y"
         << setw(12) << "x'" << setw(12) << "y'"  << endl; 

/* integration of ODE */
    while (ti <= tmax)
    {

     file << setw(12) << ti   << setw(12) << xi[0] << setw(12) << xi[1]
          << setw(12) << xi[2]<< setw(12) << xi[3] << endl;

     if (xi[1] < 0.0) break;     // stop if the projectile is below the ground 

     tf = ti + dt;
/*=================================*/
     rk4_dn1(dnx, ti, tf, xi, xf, n);
/*=================================*/

// prepare for the next step
        ti = tf;
        for (i = 0; i<=n-1; i = i+1)
        { 
           xi[i] = xf[i];
        }
   }

    system ("pause");
    return 0;
}

/*============================================================== 
  System of first order differential equations for the RK solver
  
  For a system of n first-order ODEs
  x [] array - x values
  dx[] array - dx/dt values
  
  For a system of n/2 second order ODEs follow the agreement
  In:  x[] array 
  # first n/2 elements are x
  # last  n/2 elements are dx/dt
  Out: dx[] array
  # first n/2 elements are first order derivatives  (or dx/dt)
  # last  n/2 elements are second order derivatives (or d2x/dt2)
  example: 2D projectile motion in (x,y) plane
  In           Out
  x[0] - x     dx[0] - x'
  x[1] - y     dx[1] - y'
  x[2] - x'    dx[2] - x"
  x[3] - y'    dx[3] - y"
==============================================================*/

    double dnx(double t, double x[], double dx[], int n)
{
/* first order */
    dx[0] = x[2];
    dx[1] = x[3];
/* second order */
    dx[2] = 0.0;
    dx[3] = (-1.0)*g;
}    

/*==========================================================
 rk4_dn1.cpp: Solution of a system of n first-order ODE
 method:      Runge-Kutta 4th-order
 written by: Alex Godunov
 last revision: 7 October 2009
----------------------------------------------------------
 call ...
 dnx(t,x[],dx[],n)- functions dx/dt   (supplied by a user)
 input ...
 ti    - initial time
 tf    - solution time
 xi[]  - initial values 
 n     - number of first order equations
 output ...
 xf[]  - solutions
==========================================================*/
double rk4_dn1(double(dnx)(double, double [], double [], int), 
               double ti, double tf, double xi[], double xf[], int n)
{
      double h, t, x[n], dx[n];
      double k1[n],k2[n],k3[n],k4[n];
      int j;

      h = tf-ti;
      t = ti;
//k1
      dnx(t, xi, dx, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k1[j] = h*dx[j];
          x[j]  = xi[j] + k1[j]/2.0;  
        }      
//k2
      dnx(t+h/2.0, x, dx, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k2[j] = h*dx[j];
          x[j]  = xi[j] + k2[j]/2.0;  
        }
//k3
      dnx(t+h/2.0, x, dx, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k3[j] = h*dx[j];
          x[j]  = xi[j] + k3[j];  
        }      
//k4 and result      
      dnx(t+h, x, dx, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k4[j] = h*dx[j];
          xf[j] = xi[j] + k1[j]/6.0+k2[j]/3.0+k3[j]/3.0+k4[j]/6.0;
        }      
    return 0.0;
}


