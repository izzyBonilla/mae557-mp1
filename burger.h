// burger.h
#ifndef BURGER_H_
#define BURGER_H_

#include <vector>

typedef std::vector<double> Vec;

struct pdeParams {
    int nx;
    int nsteps;
    double v;
    double dx;
    double dt;
};

int init(double* u, struct pdeParams params);

int implicit(double* u, struct pdeParams params);

int expl(double* u, struct pdeParams params);

#endif // BURGER_H_