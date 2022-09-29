#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

class Integrator
{
    public:
        virtual ~Integrator() {}

        virtual int Step(double t, double *x) = 0;
};

#endif // INTEGRATOR_HPP_