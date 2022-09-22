#include <iostream>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "burger.h"

int main(int argc, char *argv[]) {

    // struct to hold all PDE integration parameters
    struct pdeParams params;

    int scheme;

    if(argc==5) {
        params.nx = atoi(argv[1]);
        params.nsteps = atoi(argv[2]);
        params.v = atof(argv[3]);
        scheme = atoi(argv[4]);
    } else {
        std::cout << "Please supply 4 arguments specifying nx, nsteps, v, 
                      and which scheme: 1 for explicit, 2 for implicit \n";
        exit(1);
    }

    // specify length of domain and the final time
    double L = 2*M_PI;
    double tf = 5;

    params.dx = L/params.nx;
    params.dt = tf/params.nsteps;


    Vec u(params.nx);
    
    // apply initial conditions
    init(u.data(),params);

    // explicit, upwind scheme
    expl(u.data(),params);

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
    int nx = params.nx;
    double dx = params.dx;

    for(int i = 0; i < nx; ++i) {
        double x = dx*i;
        u[i] = sin(x)*exp(-(pow(x-M_PI,2)));
    }

    return 0;
}

int implicit(double* u, struct pdeParams params) {
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
        for (j = 1; j < params.nx-1; ++j) {
            coeff = (u[j]>0); // if u[j] < 0, then use forward difference
            f = v/(dx*dx)*(u[j+1]-2*u[j]+u[j-1])-u[j]/dx*(u[j-coeff+1]-u[j-coeff]);
            u[j] += dt*f;
        }
    }
    return 0;
}