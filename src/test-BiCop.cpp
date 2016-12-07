#include "include/BiCop.hpp"

// functions for unit testing in R
extern "C" {
    void test_bicop(int* family, double* par, double* par2) {
        BiCop obj(*family, *par, *par2);
        if (*family == -1)
            obj = BiCop();
        *family = obj.getFamily();
        *par = obj.getPar();
        *par2 = obj.getPar2();
    }
}
extern "C" {
    void test_bicop_set(int* family, double* par, double* par2) {
        BiCop obj;
        obj.setFamily(*family);
        obj.setPar(*par);
        obj.setPar2(*par2);
    }
}
extern "C" {
    void test_bicop_hfunc1(double* u1, double* u2, double* out, int* n,
                           int* family, double* par, double* par2) {
        BiCop obj(*family, *par, *par2);
        obj.hFunc1(u1, u2, out, n);
    }
}
extern "C" {
    void test_bicop_hfunc2(double* u1, double* u2, double* out, int* n,
                           int* family, double* par, double* par2) {
        BiCop obj(*family, *par, *par2);
        obj.hFunc2(u1, u2, out, n);
    }
}
extern "C" {
    void test_bicop_pdf(double* u1, double* u2, double* out, int* n,
                           int* family, double* par, double* par2) {
        BiCop obj(*family, *par, *par2);
        obj.PDF(u1, u2, out, n);
    }
}
