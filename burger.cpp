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
    int iter = 100;

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
        
        // Newton iteration
        for(int k = 1; k < iter; ++k) {
            // gameplan: evaluate a
            // invert jac
            //
            for(int j = 1; j < nwg-1; ++j) {
            }
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

    // periodicity automatically enforced with ghost points
    for(int i = 1; i < nwg-1; ++i) {
        jac[nwg*i+i] = 1 + dt*(1/(2*dx)*(u[i+1]-u[i-1])+2*v/(dx*dx));   // main diagonal
        jac[nwg*i+i-1] = dt*(1/(2*dx)*(-u[i])-v/(dx*dx));              // sub diagonal
        jac[nwg*i+i+1] = dt*(1/(2*dx)*(u[i])-v/(dx*dx));              // super diagonal
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