/*
 Solver for second order Ordinary Differential Equations

 x"(t) = f(t,x,x')    equation
 x(ti) = xi           initial value 1
 x'(ti)= vi           initial value 2

 the second order ODE is solved as a set of two first order ODEs
 x'(t) = f1(t,x)
 x"(t) = f2(t,x,x')

 Methods (select one by a key)
 key = 0; simple Euler
 key = 1; modified Euler (predictor-corrector)
 key = 2; 4-th order Runge-Kutta

 Alex Godunov: Last revision March 2007
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

/* function prototypes */
double f1(double, double, double);
double f2(double, double, double);
double euler2d(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);
double euler2m(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);
double rk4_2nd(double(*)(double, double, double),
               double(*)(double, double, double),
               double, double, double, double,
               double&, double&);

int main()
{
    double ti, xi, vi, tf, xf, vf, dt, tmax;
    double energy;
    int key;
    const string method[3] = {"simple Euler","modified Euler","4th order Runge-Kutta"};

/* output: file and formats */
    ofstream file;
    file.open ("ode02.dat");
    file.precision(6);
    file.setf(ios::fixed | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);

/* initial information */
    key =  2;             // select a method (key = 0, 1, 2)
    ti = 0.0;             // initial value for variable
    xi = 0.0;             // initial value for function x(t)
    vi = 1.0;             // initial
    dt = 0.1;             // step size for integration
    tmax = 12.0;          // integrate from ti till tmax
/* end of initial information */

    file << setw(30) << method[key] << endl;
    file << setw(12) << "t" << setw(12) <<"x"<< setw(12) << "x'" << endl;
    file << setw(12) << ti << setw(12) << xi << setw(12) << vi   << endl;

/* integration of ODE */
    while (ti <= tmax)
    {
        tf = ti + dt;
        if(key == 0) euler2d(f1,f2,ti,xi,vi,tf,xf,vf);
        if(key == 1) euler2m(f1,f2,ti,xi,vi,tf,xf,vf);
        if(key == 2) rk4_2nd(f1,f2,ti,xi,vi,tf,xf,vf);

        file << setw(12) << tf << setw(12) << xf  << setw(12) << vf << endl;
        ti = tf;
        xi = xf;
        vi = vf;
    }
    system ("pause");
    return 0;
}

/*
  Definition of the x'(t) = f1(t,x,x') = x' by the definition
*/
    double f1(double t, double x, double v)
{
    double d1x;
    d1x = v;
    return d1x;
}
/*
 *  Definition of the x"(t) = f2(t,x,x')
*/
    double f2(double t, double x, double v)
{
    double d2x;
    d2x = (-1.0)*x;  //simple harmonic oscillator
    return d2x;
}

/*-----------------------------------------------------------------------
 Program to solve 1st order Ordinary Differential Equations
 x'(t) = d1x(t,x)
 x"(t) = d2x(t,x)
 method:    simple Euler method
 input ...
 d1x(t,x) - first order derivative: supplied by a user
 d2x(t,x) - second order derivative: supplied by a user
 ti  - initial value for an independent variable (t)
 xi  - initial value for x(t)
 vi  - initial value for x'(t)
 tf  - find solution for this point t
 output ...
 xf and vf - solutions at point tf, i.e. x(tf) and x'(tf)
-----------------------------------------------------------------------*/
double euler2d(double(*d1x)(double, double, double),
               double(*d2x)(double, double, double),
               double ti, double xi, double vi, double tf,
               double& xf, double& vf)
{
    xf = xi + d1x(ti,xi,vi)*(tf-ti);
    vf = vi + d2x(ti,xi,vi)*(tf-ti);
   return 0.0;
}

/*-----------------------------------------------------------------------
 Program to solve 1st order Ordinary Differential Equations
 x'(t) = d1x(t,x)
 x"(t) = d2x(t,x)
 method: modified Euler method (predictor-corrector)
 input ...
 d1x(t,x) - first order derivative: supplied by a user
 d2x(t,x) - second order derivative: supplied by a user
 ti  - initial value for an independent variable (t)
 xi  - initial value for x(t)
 vi  - initial value for x'(t)
 tf  - find solution for this point t
 output ...
 xf and vf - solutions at point tf, i.e. x(tf) and x'(tf)
-----------------------------------------------------------------------*/
double euler2m(double(*d1x)(double, double, double),
               double(*d2x)(double, double, double),
               double ti, double xi, double vi, double tf,
               double& xf, double& vf)
{
    xf = xi + d1x(ti,xi,vi)*(tf-ti);
    vf = vi + d2x(ti,xi,vi)*(tf-ti);
/* correction */
    xf = xi + (d1x(ti,xi,vi)+d1x(ti,xf,vf))*0.5*(tf-ti);
    vf = vi + (d2x(ti,xi,vi)+d2x(ti,xf,vf))*0.5*(tf-ti);
   return 0.0;
}

/*-----------------------------------------------------------------------
 Program to solve 1st order Ordinary Differential Equations
 x'(t) = d1x(t,x)
 x"(t) = d2x(t,x)
 method: 4th-order Runge-Kutta method
 input ...
 d1x(t,x) - first order derivative: supplied by a user
 d2x(t,x) - second order derivative: supplied by a user
 ti  - initial value for an independent variable (t)
 xi  - initial value for x(t)
 vi  - initial value for x'(t)
 tf  - find solution for this point t
 output ...
 xf and vf - solutions at point tf, i.e. x(tf) and x'(tf)
-----------------------------------------------------------------------*/
double rk4_2nd(double(*d1x)(double, double, double),
               double(*d2x)(double, double, double),
               double ti, double xi, double vi, double tf,
               double& xf, double& vf)
{
      double h,t,k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v;

      h = tf-ti;
      t = ti;

      k1x = h*d1x(t,xi,vi);
      k1v = h*d2x(t,xi,vi);

      k2x = h*d1x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);
      k2v = h*d2x(t+h/2.0,xi+k1x/2.0,vi+k1v/2.0);

      k3x = h*d1x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);
      k3v = h*d2x(t+h/2.0,xi+k2x/2.0,vi+k2v/2.0);

      k4x = h*d1x(t+h,xi+k3x,vi+k3v);
      k4v = h*d2x(t+h,xi+k3x,vi+k3v);

      xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)/6.0;
      vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)/6.0;

      return 0.0;
}
