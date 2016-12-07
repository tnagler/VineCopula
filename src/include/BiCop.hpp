#ifndef BICOP_HPP
#define BICOP_HPP

#include <math.h>

extern "C" {
#include "hfunc.h"
}

extern "C" {
#include "likelihood.h"
}

class BiCop
{
public:
    // constructors --------------------------------

    // default constructor (independence copula)
    BiCop();
    // construct with family and parameters
    BiCop(const int &family, const double &par, const double &par2);
    // calculate number of parameters
    int calculateNPars(const int &family) ;

    // getters and setters --------------------------------

    int getFamily() const;
    double getPar() const;
    double getPar2() const;
    double getNPars() const;
    void setPar(const double &par);
    void setFamily(const int &family);
    void setPar2(const double &par2);

    // Family-specific functions --------------------------------

    // hfunctions: the conditioning variable is put second
    void hFunc1(double* u1, double* u2, double* out, int* n);
    void hFunc2(double* u1, double* u2, double* out, int* n);

    // PDF
    void PDF(double* u1, double* u2, double* out, int* n);

    // fit statistics
    void logLik(double* u1, double* u2, double* out, int* n);
    void AIC(double* u1, double* u2, double* out, int* n);
    void BIC(double* u1, double* u2, double* out, int* n);

private:
    int _family;     // copula family
    double _par;     // first copula parameter
    double _par2;    // second copula parameter
    int _npars;      // number of parameters
};

#endif
