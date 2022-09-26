// burger.h
#ifndef BURGER_HPP_
#define BURGER_HPP_

#include <vector>

typedef std::vector<double> Vec;

struct pdeParams {
    int nx;
    int nwg;
    int nsteps;
    double v;
    double dx;
    double dt;
};

int init(double* u, struct pdeParams params);

int implicit(double* u, struct pdeParams params);

int jacobian(double* jac, const double* u, struct pdeParams params);

int expl(double* u, struct pdeParams params);

#endif // BURGER_HPP_