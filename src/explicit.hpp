#ifndef EXPLICIT_HPP_
#define EXPLICIT_HPP_

#include "integrator.hpp"

class Explicit: public Integrator
{
    public:
        Explicit(double* x, double dt, double dx, double v, int nx, int nwg)
        ~Explict();
}

#endif // EXPLICIT_HPP_
