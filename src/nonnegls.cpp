#include <iostream>
//#include "PatWie/meta.h"
//#include "PatWie/boundedproblem.h"
//#include "PatWie/solver/lbfgsbsolver.h"

#include "lbgfsb.hpp"

using namespace lbfgsb;

// we will solve ||Xb-y|| s.t. b>=0
class NonNegativeLeastSquares : public BoundedProblem<double> {
  public:
    typedef BoundedProblem<double> Superclass;
    typedef Superclass::TVector TVector;
    typedef Superclass::THessian TMatrix;

    const TMatrix X;
    const TVector y;

  public:
    NonNegativeLeastSquares(const TMatrix &X_, const TVector y_) :
        Superclass(X_.rows()),
        X(X_), y(y_) {}

    double value(const TVector &beta) {
        return (X*beta-y).dot(X*beta-y);
    }


    void gradient(const TVector &beta, TVector &grad) {
        grad = X.transpose()*2*(X*beta-y);
    }

};

int main(int argc, char const *argv[]) {

    const size_t DIM = 4;
    const size_t NUM = 10;
    typedef NonNegativeLeastSquares TNNLS;
    typedef TNNLS::TVector TVector;
    typedef TNNLS::TMatrix TMatrix;

    // create model X*b for arbitrary b
    TMatrix X         = TMatrix::Random(NUM, DIM);
    TVector true_beta = TVector::Random(DIM);
    TMatrix y         = X*true_beta;

    // perform non-negative least squares
    TNNLS f(X, y);
    f.setLowerBound(TVector::Zero(DIM));
    // create initial guess (make sure it's valid >= 0)
    TVector beta = TVector::Random(DIM);
    beta = (beta.array() < 0).select(-beta, beta);
    std::cout << "true b  = " << true_beta.transpose() << "\tloss:" << f(true_beta) << std::endl;
    std::cout << "start b = " << beta.transpose() << "\tloss:" << f(beta) << std::endl;
    // init L-BFGS-B for box-constrained solving
    LbfgsbSolver<TNNLS> solver;
    solver.minimize(f, beta);
    std::cout << "final b = " << beta.transpose() << "\tloss:" << f(beta) << std::endl;

    return 0;
}
