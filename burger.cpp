#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "burger.hpp"

#define EXPLICIT 1
#define IMPLICIT 2
#define NGHOST 2

int main(int argc, char *argv[]) {

    // struct to hold all PDE integration parameters
    struct pdeParams params;

    int scheme;

    if(argc==5) {
        params.nx = atoi(argv[1]);
        params.nwg = params.nx + NGHOST;
        params.nsteps = atoi(argv[2]);
        params.v = atof(argv[3]);
        scheme = atoi(argv[4]);
    } else {
        std::cout << "Please supply 4 arguments specifying nx, nsteps, v,";
        std::cout << "and which scheme: 1 for explicit, 2 for implicit \n";
        exit(1);
    }

    // specify length of domain and the final time
    double L = 2*M_PI;
    double tf = 5;

    params.dx = L/params.nx;
    params.dt = tf/params.nsteps;


    Vec u(params.nwg);

    // apply initial conditions
    init(u.data(),params);

    if(scheme == EXPLICIT) {
        // explicit, upwind scheme
        expl(u.data(),params);
    } else if(scheme == IMPLICIT) {
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
    int nx = params.nx;
    int nsteps = params.nsteps;
    double v = params.v;
    double dx = params.dx;
    double dt = params.dt;
    
    Vec jac(nwg*nwg);
    Vec u0(nwg);
    Vec a(nwg);

    // rhs value placeholder
    double f;

    // number of newton iterations per step
    int iter = 10;

    // gameplan: do an fe step as a first guess, then do newton iterations to find next timestep
    for(int i = 0; i < nsteps; ++i) {
        // Periodic Boundary Conditions
        u[0] = u[nwg-4];
        u[1] = u[nwg-3];
        u[nwg-2] = u[2];
        u[nwg-1] = u[3];

        // FE Step
        for(int j = 1; j < nwg-1; ++j) {
            f = v/(dx*dx)*(u[j+1]-2*u[j]+u[j-1])-u[j]/(2*dx)*(u[j+1]-u[j-1]);
            u0[j] = u[j] + dt*f;
        }

        // calculate jacobian based on guess u0
        jacobian(jac.data(),u0.data(),params);

        // newton iterations
        newtonMethod(jac.data(), u, u0.data(), params);
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

    for(int i = 1; i < nwg-1; ++i) {
        jac[nwg*i+i] = 1 + dt*(1/(2*dx)*(u[i+1]-u[i-1])+2*v/(dx*dx));   // main diagonal
        jac[nwg*i+i-1] = dt*(1/(2*dx)*(-u[i])-v/(dx*dx));              // sub diagonal
        jac[nwg*i+i+1] = dt*(1/(2*dx)*(u[i])-v/(dx*dx));              // super diagonal
    }

    // hard coding periodicity like a chump, exclude ghost points
    // 1st row 1st column element move to 1st row nth column
    jac[nwg+nx] = jac[nwg];
    // nth row nth column element move to nth row 1st column
    jac[nwg*(nwg-2)+1] = jac[nwg*(nwg-1)-1]; 

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

    int iter = 10;

    Vec a(nwg);

    // main loop
    for(int k = 0; k < iter; ++k) {
        // evaluate nonlinear function a, negate here for ease
        // please note if you ever reuse this, you will probably have to rewrite what a is
        for(int i = 1; i < nwg-1; ++i) {
            a[i] = -(u0[i]+dt*(u[i]/(2*dx)*(u[i+1]-u[i-1])-v/(dx*dx)*(u[i+1]-2*u[i]+u[i-1]))-u[i]);
        }

        // solve cyclic will modify u0 in place
        solveCyclic(jac, u0, a.data(),params);
        for(int i = 0; i < nwg; ++i) {
            u0[i] = u0[i] + u[i];
        }
    }

    u = u0;

    return 0;
}

int solveCyclic(double* tc, double* u0, double* a, struct pdeParams params) {
    // solves tridiagonal linear system governed as tc*u0 = a

    // unpack useful dimensions from params
    int nx = params.nx;
    int nwg = params.nwg;

    int nc = nx-1; // dimension of condensed matrix

    return 0;
}
int condense(double* a, double* ac, struct pdeParams params) {
    // given cyclic tridiagonal matrix a whose dimensions are defined in params,
    // condense a to ac by ignoring first and last row and columns

    // a is nwg*nwg, thus the first and last row and columns are ghost points
    int nwg = params.nwg;

    for(int i = 2; i < nwg-2; ++i) {
        ac[i-2] = a[i];
    }

    return 0;
}