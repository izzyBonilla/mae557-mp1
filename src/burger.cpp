#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "burger.hpp"

#define EXPLICIT 1
#define IMPLICIT 2
#define NGHOST 2
#define L 2*M_PI

int main(int argc, char *argv[]) {

    // struct to hold all PDE integration parameters
    struct pdeParams params;

    int scheme;

    double tf;

    if(argc==6) { // unpacking command line arguments
        params.nx = atoi(argv[1]);
        params.nwg = params.nx + NGHOST;
        params.nsteps = atoi(argv[2]);
        tf = atof(argv[3]);
        params.v = atof(argv[4]);
        scheme = atoi(argv[5]);

    } else { // error handling for wrong number of arguments
        std::cout << "Please supply 5 arguments specifying nx, nsteps, tf, v,";
        std::cout << "and which scheme: 1 for explicit, 2 for implicit \n";
        exit(1);
    }

    params.dx = L/params.nx;
    params.dt = tf/params.nsteps;

    // state vector
    Vec u(params.nwg);

    // apply initial conditions
    init(u.data(),params);
    
    if(scheme == EXPLICIT) {
        // explicit, upwind scheme
        expl(u.data(),params);
    } else if(scheme == IMPLICIT) {
        // implicit, centered scheme
        implicit(u.data(),params);
    }
    
    
    // write to file
    std::ofstream arr_file("u.dat");
    if(arr_file.is_open()) {
        for(int i = 0; i < params.nx; ++i) {
            arr_file << params.dx*i << "\t" << u[i] << "\n";
        }
    } else {
        std::cout << "Can't write to file u.dat \n";
        exit(1);
    }

    return 0;
     
}

int init(double* u, struct pdeParams params) {
    // Writes initial condition to state vector u
    int nwg = params.nwg;
    double dx = params.dx;

    for(int i = 0; i < nwg; ++i) {
        double x = dx*i;
        u[i] = sin(x)*exp(-(pow(x-M_PI,2)));
    }

    return 0;
}

int expl(double* u, struct pdeParams params) {
    // implementation of explicit scheme using 
    // 2-point upwinded differencing and forward euler

    double v = params.v;
    double dx = params.dx;
    double dt = params.dt;

    int i;
    int j;

    double f;
    int coeff;

    // main loop
    // could parallelize this easily but it's not worth it
    for(i = 0; i < params.nsteps; ++i){
        // perodic boundary conditions with ghost points
        u[0] = u[params.nwg-4];
        u[1] = u[params.nwg-3];
        u[params.nwg-2] = u[2];
        u[params.nwg-1] = u[3];
        for (j = 1; j < params.nwg-2; ++j) {
            coeff = (u[j]>0); // if u[j] > 0, then use backward difference else forward difference
            f = v/(dx*dx)*(u[j+1]-2*u[j]+u[j-1])-u[j]/dx*(u[j-coeff+1]-u[j-coeff]);
            u[j] += dt*f;
        }
    }
    return 0;
}

int implicit(double* u, struct pdeParams params) {

    // unpack struct for convenient naming
    int nwg = params.nwg;
    int nsteps = params.nsteps;
    double v = params.v;
    double dx = params.dx;
    double dt = params.dt;
    
    Vec jac(nwg*nwg);
    Vec u0(nwg);

    // rhs value placeholder
    double f;

    // gameplan: do an fe step as a first guess, then do newton iterations to find next timestep
    for(int i = 0; i < nsteps; ++i) {

        // FE Step
        for(int j = 1; j < nwg-1; ++j) {
            f = v/(dx*dx)*(u[j+1]-2*u[j]+u[j-1])-u[j]/(2*dx)*(u[j+1]-u[j-1]);
            u0[j] = u[j] + dt*f;
        }

        // Periodic Boundary Conditions
        u[0] = u[nwg-4];
        u[1] = u[nwg-3];
        u[nwg-2] = u[2];
        u[nwg-1] = u[3];
        // newton iterations
        newtonMethod(jac.data(), u, u0.data(), params);
        
        // now assign u to guess u0
        for(int i = 0; i < nwg; ++i) {
            u[i] = u0[i];
        }
    }    

    return 0;
}

int jacobian(double* jac, const double* u, struct pdeParams params) {
    // given state vector u and PDE parameters params, determine
    // jacobian of nonlinear system

    int nx = params.nx;
    int nwg = params.nwg;
    double dx = params.dx;
    double dt = params.dt;
    double v = params.v;

    for(int i = 1; i < nwg - 1; ++i) {
        jac[(i*nwg)+i] = 1 + dt*(1/(2*dx)*(u[i+1]-u[i-1])+2*v/(dx*dx));   // main diagonal
        jac[(i*nwg)+i-1] = dt*(1/(2*dx)*u[i] - v/(dx*dx));               // sub diagonal
        jac[(i*nwg)+i+1] = dt*(1/(2*dx)*(-u[i]) - v/(dx*dx));           // super diagonal
    }

    // hard coding periodicity like a chump, exclude ghost points
    // 1st row 1st column element move to 1st row nth column
    jac[nwg+nx] = jac[nwg];
    jac[nwg] = 0;
    // nth row nth column element move to nth row 1st column
    jac[nwg*(nwg-2)+1] = jac[nwg*(nwg-1)-1];
    jac[nwg*(nwg-1)-1] = 0; 

    return 0;
}

int newtonMethod(double* jac, double* u, double* u0, struct pdeParams params) {

    // given a jacobian jac, state vector u, and initial guess u0,
    // iterate a pre-specified number of times to solve for u_n+1

    // unpack params
    int nwg = params.nwg;
    int nx = params.nx;
    double dx = params.dx;
    double dt = params.dt;
    double v = params.v;

    int iter = 10;

    Vec a(nwg);

    // main loop
    for(int k = 0; k < iter; ++k) {
        // evaluate nonlinear function a, negate here for ease
        // please note if you ever reuse this, you will probably have to rewrite what a is

        // calculate jacobian based on guess u0
        jacobian(jac,u0,params);

        for(int i = 1; i < nwg-1; ++i) {
            a[i] = -(u0[i]+dt*(u0[i]/(2*dx)*(u0[i+1]-u0[i-1])-v/(dx*dx)*(u0[i+1]-2*u0[i]+u0[i-1]))-u0[i]);
        }

        // solve cyclic will modify u0 in place
        solveCyclic(jac, u0, a.data(),params);

        // now have u0 = -Ja(u) = u_n+1 - u_n
        for(int i = 1; i < nwg-1; ++i) {
            u0[i] = u0[i] + u[i];
        }
    }

    return 0;
}

int solveCyclic(double* tc, double* u0, double* a, struct pdeParams params) {
    // solves tridiagonal linear system governed as tc*u0 = a

    // inputs:
    // tc: size nwg*nwg double, jacobian
    // a: size nwg double, nonlinear input
    // params: struct with PDE parameters

    // outputs:
    // u0: size nwg double, represents next guess for newton iteration

    // unpack useful dimensions from params
    int nx = params.nx;
    int nwg = params.nwg;

    int nc = nx-1; // dimension of condensed matrix

    Vec tc_cond(nc*nc);

    condense(tc,tc_cond.data(),params);

    // LU Factorization
    Vec l(nc-1);
    Vec s(nc-1);
    Vec d(nc);

    d[0] = tc_cond[0];
    s[0] = tc_cond[1];

    for(int i = 0; i < nc-2; ++i) {
        l[i] = tc_cond[(i+1)*nc+i] / d[i];
        d[i+1] = tc_cond[(i+1)*nc+i+1] - l[i]*s[i];
        s[i+1] = tc_cond[(i+1)*nc+i+2];
    }

    l[nc-2] = tc_cond[nc*nc-2]/d[nc-2];
    d[nc-1] = tc_cond[nc*nc-1] - l[nc-2]*s[nc-2];

    // Construct known vectors q1 and q2
    Vec q1(nc); // q1 will represent knowns from 0 to nx-1
    Vec q2(nc); // q2 will capture the condensed information

        // direct copy of a vector
    for(int i = 0; i < nc; ++i) {
        q1[i] = a[i+1]; // skip ghost point in a
    }

    q2[0] = -tc[nwg+nx];
    q2[nc-1] = -tc[(nwg-3)*nwg+nx]; // access nth element of n-1th row with 2 ghost points

    // Partial Solutions y1 and y2 for forward substitution step
    Vec y1(nc);
    Vec y2(nc);

    // Full Solutions x1 and x2
    Vec x1(nc);
    Vec x2(nc);

    // Forward Substitution
    y1[0] = q1[0];
    y2[0] = q2[0];

    for(int i = 1; i < nc; ++i) {
        y1[i] = q1[i] - l[i-1]*y1[i-1];
        y2[i] = q2[i] - l[i-1]*y2[i-1];
    }

    // Backward Substitution

    x1[nc-1] = y1[nc-1]/d[nc-1];
    x2[nc-1] = y2[nc-1]/d[nc-1];

    for(int i = nc-2; i >= 0; --i) {
        x1[i] = (y1[i] - s[i]*x1[i+1])/d[i];
        x2[i] = (y2[i] - s[i]*x2[i+1])/d[i];
    }

    // calculate u at the last point
    double num = (a[nx]-tc[nwg*(nwg-2)+1]*x1[0]-tc[nwg*(nwg-1)-3]*x1[nc-1]);
    double denom = (tc[nwg*(nwg-1)-2]+tc[nwg*(nwg-2)+1]*x2[0]+tc[nwg*(nwg-1)-3]*x2[nc-1]);
    u0[nx] = num/denom;

    for(int i = 1; i < nx; ++i) {
        u0[i] = x1[i-1] + x2[i-1]*u0[nx];
    }

    return 0;
}
int condense(double* a, double* ac, struct pdeParams params) {
    // given cyclic tridiagonal matrix a whose dimensions are defined in params,
    // condense a to ac by removing last row and last column

    // a is nwg*nwg, thus the first and last row and columns are ghost points
    int nx = params.nx;
    int nwg = params.nwg;
    int nc = params.nx-1;

    for(int i = 0; i < nc; ++i) {
        for(int j = 0; j < nc; ++j) {
            ac[i*nc+j] = a[nwg*(i+1)+j+1];
        }
    }

    return 0;
}