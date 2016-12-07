#include "include/BiCop.hpp"

// constructors --------------------------------

// default constructor (independence copula)
BiCop::BiCop()
{
    _family = 0;
    _par = 0;
    _par2 = 0;
    _npars = 0;
}
// construct with family and parameters
BiCop::BiCop(const int &family, const double &par, const double &par2) {
    _family = family;
    _par = par;
    _par2 = par2;
    _npars = calculateNPars(family);
}
// calculate number of parameters
int BiCop::calculateNPars(const int &family) {
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

int BiCop::getFamily() const {return _family;}
double BiCop::getPar() const {return _par;}
double BiCop::getPar2() const {return _par2;}
double BiCop::getNPars() const {return _npars;}

void BiCop::setFamily(const int &family) {
    _family = family;
    _npars = calculateNPars(family);
}
void BiCop::setPar(const double &par) {_par = par;}
void BiCop::setPar2(const double &par2) {_par2 = par2;}

// Family-specific functions --------------------------------

// hfunctions: the conditioning variable is put second
void BiCop::hFunc1(double* u1, double* u2, double* out, int* n) {
    Hfunc1(&_family, n, u2, u1, &_par, &_par2, out);
};
void BiCop::hFunc2(double* u1, double* u2, double* out, int* n) {
    Hfunc2(&_family, n, u1, u2, &_par, &_par2, out);
};

// PDF
void BiCop::PDF(double* u1, double* u2, double* out, int* n) {
    PDF_seperate(&_family, n, u1, u2, &_par, &_par2, out);
};

// fit statistics
void BiCop::logLik(double* u1, double* u2, double* out, int* n) {
    LL_mod2(&_family, n, u1, u2, &_par, &_par2, out);
};
void BiCop::AIC(double* u1, double* u2, double* out, int* n) {
    LL_mod2(&_family, n, u1, u2, &_par, &_par2, out);
    *out = -2 * (*out) + 2 * _npars;
};
void BiCop::BIC(double* u1, double* u2, double* out, int* n) {
    LL_mod2(&_family, n, u1, u2, &_par, &_par2, out);
    *out = -2 * (*out) + log(*n) * _npars;
};
