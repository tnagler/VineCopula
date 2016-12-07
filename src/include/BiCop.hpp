#ifndef BICOP_HPP
#define BICOP_HPP

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
    BiCop()
    {
        _family = 0;
        _par = 0;
        _par2 = 0;
        _npars = 0;
    }
    // construct with family and parameters
    BiCop(const int &family, const double &par, const double &par2)
    {
        _family = family;
        _par = par;
        _par2 = par2;
        _npars = calculateNPars(family);
    }
    // calculate number of parameters
    int calculateNPars(const int &family) {
        // indepence copula has no parameters
        if (family == 0)
            return(0);
        // check if it is a one-parameter family
        int onepars[] = {1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36};
        bool isOnePar = false;
        for (int k = 0; k < sizeof(onepars) / sizeof(onepars[0]); k++) {
            if (family == onepars[k])
                return(1);
        }
        // if not, it must be a two-parameter family
        return(2);
    }

    // getters and setters --------------------------------

    int getFamily() const {return _family;}
    double getPar() const {return _par;}
    double getPar2() const {return _par2;}
    double getNPars() const {return _npars;}

    void setFamily(const int &family) {
        _family = family;
        _npars = calculateNPars(family);
    }
    void setPar(const double &par) {_par = par;}
    void setPar2(const double &par2) {_par2 = par2;}

    // Family-specific functions --------------------------------

    // hfunctions: the conditioning variable is put second
    void hFunc1(double* u1, double* u2, double* out, int* n) {
        Hfunc1(&_family, n, u2, u1, &_par, &_par2, out);
    };
    void hFunc2(double* u1, double* u2, double* out, int* n) {
        Hfunc2(&_family, n, u1, u2, &_par, &_par2, out);
    };

    // PDF and log likelihood
    void PDF(double* u1, double* u2, double* out, int* n) {
        PDF_seperate(&_family, n, u1, u2, &_par, &_par2, out);
    };
    void logLik(double* u1, double* u2, double* out, int* n) {
        LL_mod2(&_family, n, u1, u2, &_par, &_par2, out);
    };

private:
    int _family;     // copula family
    double _par;     // first copula parameter
    double _par2;    // second copula parameter
    int _npars;      // number of parameters
};

#endif
