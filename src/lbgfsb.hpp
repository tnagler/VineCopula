/* Eigen-based L-BFGS-C sover adapted to C++-03/Boost from C++11/non-Boost code at
   https://github.com/PatWie/CppNumericalSolvers
   (by Patrick Wieschollek, MIT license, original header below)
    
   * merged into one file with multi-include guards left intentionally
     to depict the original division into multiple files
   * intended as a temporary solution before we can use the original C++11 code
*/
/*
 *  This file is part of CppNumericalSolvers
 *
 *  Copyright (C) Tobias Wood @spinicist 2016
 *
 *  This Source Code form is subject to the terms of the MIT license.
 *  Please see the LICENSE file.
 *
 */

#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <list>

#include <boost/array.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

namespace lbfgsb 
{

#ifndef META_H
#define META_H

enum DebugLevelEnum { None = 0, Low, High };
struct DebugLevel { static const int None = 0, Low = 1, High = 2; };
enum StatusEnum {
    NotStarted = -1,
    Continue = 0,
    IterationLimit,
    XDeltaTolerance,
    FDeltaTolerance,
    GradNormTolerance,
    Condition
};
struct Status { static const int
    NotStarted = -1,
    Continue = 0,
    IterationLimit = 1,
    XDeltaTolerance = 2,
    FDeltaTolerance = 3,
    GradNormTolerance = 4,
    Condition = 5;
};

template<typename T>
class Criteria {
public:
    size_t iterations; //!< Maximum number of iterations
    T xDelta;          //!< Minimum change in parameter vector
    T fDelta;          //!< Minimum change in cost function
    T gradNorm;        //!< Minimum norm of gradient vector
    T condition;       //!< Maximum condition number of Hessian

    Criteria() {
        reset();
    }

    static Criteria defaults() {
        Criteria d;
        d.iterations = 10000;
        d.xDelta = 0;
        d.fDelta = 0;
        d.gradNorm = 1e-4;
        d.condition = 0;
        return d;
    }

    void reset() {
        iterations = 0;
        xDelta = 0;
        fDelta = 0;
        gradNorm = 0;
        condition = 0;
    }

    void print(std::ostream &os) const {
        os << "Iterations: " << iterations << std::endl;
        os << "xDelta:     " << xDelta << std::endl;
        os << "fDelta:     " << fDelta << std::endl;
        os << "GradNorm:   " << gradNorm << std::endl;
        os << "Condition:  " << condition << std::endl;
    }
};

template<typename T>
int checkConvergence(const Criteria<T> &stop, const Criteria<T> &current) {

    if ((stop.iterations > 0) && (current.iterations > stop.iterations)) {
        return Status::IterationLimit;
    }
    if ((stop.xDelta > 0) && (current.xDelta < stop.xDelta)) {
        return Status::XDeltaTolerance;
    }
    if ((stop.fDelta > 0) && (current.fDelta < stop.fDelta)) {
        return Status::FDeltaTolerance;
    }
    if ((stop.gradNorm > 0) && (current.gradNorm < stop.gradNorm)) {
        return Status::GradNormTolerance;
    }
    if ((stop.condition > 0) && (current.condition > stop.condition)) {
        return Status::Condition;
    }
    return Status::Continue;
}

inline std::ostream &operator<<(std::ostream &os, const int &s) {
    switch (s) {
        case Status::NotStarted: os << "Solver not started."; break;
        case Status::Continue:   os << "Convergence criteria not reached."; break;
        case Status::IterationLimit: os << "Iteration limit reached."; break;
        case Status::XDeltaTolerance: os << "Change in parameter vector too small."; break;
        case Status::FDeltaTolerance: os << "Change in cost function value too small."; break;
        case Status::GradNormTolerance: os << "Gradient vector norm too small."; break;
        case Status::Condition: os << "Condition of Hessian/Covariance matrix too large."; break;
    }
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Criteria<T> &c) {
    c.print(os);
    return os;
}

#endif /* META_H */


#ifndef PROBLEM_H
#define PROBLEM_H

template<typename Scalar_, int Dim_ = Eigen::Dynamic>
class Problem {
 public:
  static const int Dim = Dim_;
  typedef Scalar_ Scalar;
  typedef Eigen::Matrix<Scalar, Dim, 1> TVector;
  typedef Eigen::Matrix<Scalar, Dim, Dim> THessian;
  typedef Criteria<Scalar> TCriteria;
  typedef typename TVector::Index TIndex;

 public:
  Problem() {}
  virtual ~Problem() {}

  virtual bool callback(const Criteria<Scalar> &state, const TVector &x) {
    return true;
  }

  /**
   * @brief returns objective value in x
   * @details [long description]
   *
   * @param x [description]
   * @return [description]
   */
  virtual Scalar value(const  TVector &x) = 0;
  /**
   * @brief overload value for nice syntax
   * @details [long description]
   *
   * @param x [description]
   * @return [description]
   */
  Scalar operator()(const  TVector &x) {
    return value(x);
  }
  /**
   * @brief returns gradient in x as reference parameter
   * @details should be overwritten by symbolic gradient
   *
   * @param grad [description]
   */
  virtual void gradient(const  TVector &x,  TVector &grad) {
    finiteGradient(x, grad);
  }

  /**
   * @brief This computes the hessian
   * @details should be overwritten by symbolic hessian, if solver relies on hessian
   */
  virtual void hessian(const TVector &x, THessian &hessian) {
    finiteHessian(x, hessian);
  }

  virtual bool checkGradient(const TVector &x, int accuracy = 3) {
    // TODO: check if derived class exists:
    // int(typeid(&Rosenbrock<double>::gradient) == typeid(&Problem<double>::gradient)) == 1 --> overwritten
    const TIndex D = x.rows();
    TVector actual_grad(D);
    TVector expected_grad(D);
    gradient(x, actual_grad);
    finiteGradient(x, expected_grad, accuracy);
    for (TIndex d = 0; d < D; ++d) {
      Scalar scale = std::max((std::max(fabs(actual_grad[d]), fabs(expected_grad[d]))), 1.);
      if(fabs(actual_grad[d]-expected_grad[d])>1e-2 * scale)
        return false;
    }
    return true;

  }

  virtual bool checkHessian(const TVector &x, int accuracy = 3) {
    // TODO: check if derived class exists:
    // int(typeid(&Rosenbrock<double>::gradient) == typeid(&Problem<double>::gradient)) == 1 --> overwritten
    const TIndex D = x.rows();

    THessian actual_hessian = THessian::Zero(D, D);
    THessian expected_hessian = THessian::Zero(D, D);
    hessian(x, actual_hessian);
    finiteHessian(x, expected_hessian, accuracy);
    for (TIndex d = 0; d < D; ++d) {
      for (TIndex e = 0; e < D; ++e) {
        Scalar scale = std::max(static_cast<Scalar>(std::max(fabs(actual_hessian(d, e)), fabs(expected_hessian(d, e)))), Scalar(1.));
        if(fabs(actual_hessian(d, e)- expected_hessian(d, e))>1e-1 * scale)
          return false;
      }
    }
    return true;
  }

  void finiteGradient(const  TVector &x, TVector &grad, int accuracy = 0) {
    // accuracy can be 0, 1, 2, 3
    const Scalar eps = 2.2204e-6;
    static boost::array<std::vector<Scalar>, 4> coeff;
    static boost::array<std::vector<Scalar>, 4> coeff2;
    {
      using namespace boost::assign;

      coeff[0] += 1., -1; 
      coeff[1] += 1., -8, 8, -1; 
      coeff[2] += -1., 9, -45, 45, -9, 1;
      coeff[3] += 3., -32, 168, -672, 672, -168, 32, -3;

      coeff2[0] += 1., -1;
      coeff2[1] += -2., -1, 1, 2;
      coeff2[2] += -3., -2, -1, 1, 2, 3;
      coeff2[3] += -4., -3, -2, -1, 1, 2, 3, 4; 
    }
    static const boost::array<Scalar, 4> dd = {2, 12, 60, 840};

    grad.resize(x.rows());
    TVector& xx = const_cast<TVector&>(x);

    const int innerSteps = 2*(accuracy+1);
    const Scalar ddVal = dd[accuracy]*eps;

    for (TIndex d = 0; d < x.rows(); d++) {
      grad[d] = 0;
      for (int s = 0; s < innerSteps; ++s)
      {
        Scalar tmp = xx[d];
        xx[d] += coeff2[accuracy][s]*eps;
        grad[d] += coeff[accuracy][s]*value(xx);
        xx[d] = tmp;
      }
      grad[d] /= ddVal;
    }
  }

  void finiteHessian(const TVector &x, THessian &hessian, int accuracy = 0) {
    const Scalar eps = std::numeric_limits<Scalar>::epsilon()*10e7;

    hessian.resize(x.rows(), x.rows());
    TVector& xx = const_cast<TVector&>(x);

    if(accuracy == 0) {
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar f4 = value(xx);
          xx[i] += eps;
          xx[j] += eps;
          Scalar f1 = value(xx);
          xx[j] -= eps;
          Scalar f2 = value(xx);
          xx[j] += eps;
          xx[i] -= eps;
          Scalar f3 = value(xx);
          hessian(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);

          xx[i] = tmpi;
          xx[j] = tmpj;
        }
      }
    } else {
      /*
        \displaystyle{{\frac{\partial^2{f}}{\partial{x}\partial{y}}}\approx
        \frac{1}{600\,h^2} \left[\begin{matrix}
          -63(f_{1,-2}+f_{2,-1}+f_{-2,1}+f_{-1,2})+\\
          63(f_{-1,-2}+f_{-2,-1}+f_{1,2}+f_{2,1})+\\
          44(f_{2,-2}+f_{-2,2}-f_{-2,-2}-f_{2,2})+\\
          74(f_{-1,-1}+f_{1,1}-f_{1,-1}-f_{-1,1})
        \end{matrix}\right] }
      */
      for (TIndex i = 0; i < x.rows(); i++) {
        for (TIndex j = 0; j < x.rows(); j++) {
          Scalar tmpi = xx[i];
          Scalar tmpj = xx[j];

          Scalar term_1 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -2*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -1*eps;  term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 1*eps;   term_1 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 2*eps;   term_1 += value(xx);

          Scalar term_2 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -2*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -1*eps;  term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 2*eps;   term_2 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 1*eps;   term_2 += value(xx);

          Scalar term_3 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += -2*eps;  term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += 2*eps;   term_3 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2*eps; xx[j] += -2*eps;  term_3 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2*eps;  xx[j] += 2*eps;   term_3 -= value(xx);

          Scalar term_4 = 0;
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += -1*eps;  term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += 1*eps;   term_4 += value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1*eps;  xx[j] += -1*eps;  term_4 -= value(xx);
          xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1*eps; xx[j] += 1*eps;   term_4 -= value(xx);

          xx[i] = tmpi;
          xx[j] = tmpj;

          hessian(i, j) = (-63 * term_1+63 * term_2+44 * term_3+74 * term_4)/(600.0 * eps * eps);
        }
      }
    }

  }

};

#endif /* PROBLEM_H */


#ifndef BOUNDEDPROBLEM_H
#define BOUNDEDPROBLEM_H

template<typename Scalar_, int CompileDim_ = Eigen::Dynamic>
class BoundedProblem : public Problem<Scalar_, CompileDim_> {
public:
    typedef Problem<Scalar_, CompileDim_> Superclass;
    typedef typename Superclass::Scalar Scalar;
    typedef typename Superclass::TVector TVector;

protected:
    TVector m_lowerBound;
    TVector m_upperBound;

public:
    BoundedProblem(int RunDim = CompileDim_) : Superclass() {
        TVector infBound(RunDim);
        infBound.setConstant(std::numeric_limits<Scalar>::infinity());
        m_lowerBound = -infBound;
        m_upperBound = infBound;
    }

    BoundedProblem(const TVector &l, const TVector &u) :
        Superclass(),
        m_lowerBound(l),
        m_upperBound(u)
    {}

    const TVector &lowerBound() const { return m_lowerBound; }
    void setLowerBound(const TVector &lb) { m_lowerBound = lb; }
    const TVector &upperBound() const { return m_upperBound; }
    void setUpperBound(const TVector &ub) { m_upperBound = ub; }

    void setBoxConstraint(TVector  lb, TVector  ub) {
        setLowerBound(lb);
        setUpperBound(ub);
    }
};

#endif // BOUNDEDPROBLEM_H


#ifndef ISOLVER_H_
#define ISOLVER_H_

template<typename ProblemType, int Ord>
class ISolver {
public:
    typedef typename ProblemType::Scalar Scalar;
    typedef typename ProblemType::TVector TVector;
    typedef typename ProblemType::THessian THessian;
    typedef typename ProblemType::TCriteria TCriteria;
protected:
    const int order_;// = Ord;
    TCriteria m_stop, m_current;
    int m_status;// = Status::NotStarted;
    int m_debug;// = DebugLevel::None;

public:
    virtual ~ISolver() {}
    ISolver() : order_(Ord), m_status(Status::NotStarted), m_debug(DebugLevel::None) {
        m_stop = TCriteria::defaults();
        m_current.reset();
    }

    ISolver(const TCriteria &s) : order_(Ord), m_status(Status::NotStarted), m_debug(DebugLevel::None) {
        m_stop = s;
        m_current.reset();
    }

    void setStopCriteria(const TCriteria &s) { m_stop = s; }
    const TCriteria &criteria() { return m_current; }
    const int &status() { return m_status; }
    void setDebug(const int &d) { m_debug = d; }

    /**
     * @brief minimize an objective function given a gradient (and optinal a hessian)
     * @details this is just the abstract interface
     *
     * @param x0 starting point
     * @param funObjective objective function
     * @param funGradient gradient function
     * @param funcHession hessian function
     */
    virtual void minimize(ProblemType &objFunc, TVector &x0) = 0;

};

#endif /* ISOLVER_H_ */

#ifndef MORETHUENTE_H_
#define MORETHUENTE_H_

template<typename ProblemType, int Ord>
class MoreThuente {

 public:
  typedef typename ProblemType::Scalar Scalar;
  typedef typename ProblemType::TVector TVector;

  /**
   * @brief use MoreThuente Rule for (strong) Wolfe conditiions
   * @details [long description]
   *
   * @param searchDir search direction for next update step
   * @param objFunc handle to problem
   *
   * @return step-width
   */

  static Scalar linesearch(const TVector &x, const TVector &searchDir, ProblemType &objFunc, const  Scalar alpha_init = 1.0) {
    // assume step width
    Scalar ak = alpha_init;

    Scalar fval = objFunc.value(x);
    TVector  g  = x.eval();
    objFunc.gradient(x, g);

    TVector s = searchDir.eval();
    TVector xx = x.eval();

    cvsrch(objFunc, xx, fval, g, ak, s);

    return ak;
  }

  static int cvsrch(ProblemType &objFunc, TVector &x, Scalar f, TVector &g, Scalar &stp, TVector &s) {
    // we rewrite this from MIN-LAPACK and some MATLAB code
    int info           = 0;
    int infoc          = 1;
    const Scalar xtol   = 1e-15;
    const Scalar ftol   = 1e-4;
    const Scalar gtol   = 1e-2;
    const Scalar stpmin = 1e-15;
    const Scalar stpmax = 1e15;
    const Scalar xtrapf = 4;
    const int maxfev   = 20;
    int nfev           = 0;

    Scalar dginit = g.dot(s);
    if (dginit >= 0.0) {
      // no descent direction
      // TODO: handle this case
      return -1;
    }

    bool brackt      = false;
    bool stage1      = true;

    Scalar finit      = f;
    Scalar dgtest     = ftol * dginit;
    Scalar width      = stpmax - stpmin;
    Scalar width1     = 2 * width;
    TVector wa = x.eval();

    Scalar stx        = 0.0;
    Scalar fx         = finit;
    Scalar dgx        = dginit;
    Scalar sty        = 0.0;
    Scalar fy         = finit;
    Scalar dgy        = dginit;

    Scalar stmin;
    Scalar stmax;

    while (true) {

      // make sure we stay in the interval when setting min/max-step-width
      if (brackt) {
        stmin = std::min(stx, sty);
        stmax = std::max(stx, sty);
      } else {
        stmin = stx;
        stmax = stp + xtrapf * (stp - stx);
      }

      // Force the step to be within the bounds stpmax and stpmin.
      stp = std::max(stp, stpmin);
      stp = std::min(stp, stpmax);

      // Oops, let us return the last reliable values
      if (
      (brackt && ((stp <= stmin) || (stp >= stmax)))
      || (nfev >= maxfev - 1 ) || (infoc == 0)
      || (brackt && ((stmax - stmin) <= (xtol * stmax)))) {
        stp = stx;
      }

      // test new point
      x = wa + stp * s;
      f = objFunc.value(x);
      objFunc.gradient(x, g);
      nfev++;
      Scalar dg = g.dot(s);
      Scalar ftest1 = finit + stp * dgtest;

      // all possible convergence tests
      if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0))
        info = 6;

      if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest))
        info = 5;

      if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest)))
        info = 4;

      if (nfev >= maxfev)
        info = 3;

      if (brackt & (stmax - stmin <= xtol * stmax))
        info = 2;

      if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))
        info = 1;

      // terminate when convergence reached
      if (info != 0)
        return -1;

      if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol)*dginit))
        stage1 = false;

      if (stage1 & (f <= fx) & (f > ftest1)) {
        Scalar fm = f - stp * dgtest;
        Scalar fxm = fx - stx * dgtest;
        Scalar fym = fy - sty * dgtest;
        Scalar dgm = dg - dgtest;
        Scalar dgxm = dgx - dgtest;
        Scalar dgym = dgy - dgtest;

        cstep( stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
      } else {
        // this is ugly and some variables should be moved to the class scope
        cstep( stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
      }

      if (brackt) {
        if (fabs(sty - stx) >= 0.66 * width1)
          stp = stx + 0.5 * (sty - stx);
        width1 = width;
        width = fabs(sty - stx);
      }
    }

    return 0;
  }

  static int cstep(Scalar& stx, Scalar& fx, Scalar& dx, Scalar& sty, Scalar& fy, Scalar& dy, Scalar& stp,
  Scalar& fp, Scalar& dp, bool& brackt, Scalar& stpmin, Scalar& stpmax, int& info) {
    info = 0;
    bool bound = false;

    // Check the input parameters for errors.
    if ((brackt & ((stp <= std::min(stx, sty) ) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
    | (stpmax < stpmin)) {
      return -1;
    }

    Scalar sgnd = dp * (dx / fabs(dx));

    Scalar stpf = 0;
    Scalar stpc = 0;
    Scalar stpq = 0;

    if (fp > fx) {
      info = 1;
      bound = true;
      Scalar theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max(theta, std::max(dx, dp));
      Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
      if (stp < stx)
        gamma = -gamma;
      Scalar p = (gamma - dx) + theta;
      Scalar q = ((gamma - dx) + gamma) + dp;
      Scalar r = p / q;
      stpc = stx + r * (stp - stx);
      stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
      if (fabs(stpc - stx) < fabs(stpq - stx))
        stpf = stpc;
      else
        stpf = stpc + (stpq - stpc) / 2;
      brackt = true;
    } else if (sgnd < 0.0) {
      info = 2;
      bound = false;
      Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max(theta, std::max(dx, dp));
      Scalar gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
      if (stp > stx)
        gamma = -gamma;

      Scalar p = (gamma - dp) + theta;
      Scalar q = ((gamma - dp) + gamma) + dx;
      Scalar r = p / q;
      stpc = stp + r * (stx - stp);
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (fabs(stpc - stp) > fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      brackt = true;
    } else if (fabs(dp) < fabs(dx)) {
      info = 3;
      bound = 1;
      Scalar theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      Scalar s = std::max(theta, std::max( dx, dp));
      Scalar gamma = s * sqrt(std::max(static_cast<Scalar>(0.), (theta / s) * (theta / s) - (dx / s) * (dp / s)));
      if (stp > stx)
        gamma = -gamma;
      Scalar p = (gamma - dp) + theta;
      Scalar q = (gamma + (dx - dp)) + gamma;
      Scalar r = p / q;
      if ((r < 0.0) & (gamma != 0.0)) {
        stpc = stp + r * (stx - stp);
      } else if (stp > stx) {
        stpc = stpmax;
      } else {
        stpc = stpmin;
      }
      stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (brackt) {
        if (fabs(stp - stpc) < fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }
      } else {
        if (fabs(stp - stpc) > fabs(stp - stpq)) {
          stpf = stpc;
        } else {
          stpf = stpq;
        }

      }
    } else {
      info = 4;
      bound = false;
      if (brackt) {
        Scalar theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
        Scalar s = std::max(theta, std::max(dy, dp));
        Scalar gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
        if (stp > sty)
          gamma = -gamma;

        Scalar p = (gamma - dp) + theta;
        Scalar q = ((gamma - dp) + gamma) + dy;
        Scalar r = p / q;
        stpc = stp + r * (sty - stp);
        stpf = stpc;
      } else if (stp > stx)
        stpf = stpmax;
      else {
        stpf = stpmin;
      }
    }

    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    } else {
      if (sgnd < 0.0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }

      stx = stp;
      fx = fp;
      dx = dp;
    }

    stpf = std::min(stpmax, stpf);
    stpf = std::max(stpmin, stpf);
    stp = stpf;

    if (brackt & bound) {
      if (sty > stx) {
        stp = std::min(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
      } else {
        stp = std::max(stx + static_cast<Scalar>(0.66) * (sty - stx), stp);
      }
    }

    return 0;

  }

};

#endif /* MORETHUENTE_H_ */

#ifndef LBFGSBSOLVER_H
#define LBFGSBSOLVER_H

template<typename TProblem>
class LbfgsbSolver : public ISolver<TProblem, 1> {
  public:
    typedef ISolver<TProblem, 1> Superclass;
    typedef typename Superclass::Scalar Scalar;
    typedef typename Superclass::TVector TVector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VariableTVector;
  protected:
  // last updates
  std::list<TVector> xHistory;
  // workspace matrices
  MatrixType W, M;
  Scalar theta;
  int DIM;
  int m_historySize; // = 5;

    struct sort_pred 
    {
      const std::vector< std::pair<int, Scalar> > &v;
      sort_pred(const std::vector< std::pair<int, Scalar> > &v) : v(v) {}
      bool operator()(size_t i1, size_t i2) { return v[i1].second < v[i2].second; }
    };

  /**
   * @brief sort pairs (k,v) according v ascending
   * @details [long description]
   *
   * @param v [description]
   * @return [description]
   */
  std::vector<int> sort_indexes(const std::vector< std::pair<int, Scalar> > &v) {
    std::vector<int> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
      idx[i] = v[i].first;
    sort(idx.begin(), idx.end(), sort_pred(v)
//[&v](size_t i1, size_t i2) {
//      return v[i1].second < v[i2].second;
//    }
    );
    return idx;
  }
  /**
   * @brief Algorithm CP: Computation of the generalized Cauchy point
   * @details PAGE 8
   *
   * @param c [description]
   */
  void getGeneralizedCauchyPoint(const TProblem &problem, TVector &x, TVector &g, TVector &x_cauchy, VariableTVector &c) {
    const int DIM = x.rows();
    // Given x,l,u,g, and B = \theta I-WMW
    // {all t_i} = { (idx,value), ... }
    // TODO: use "std::set" ?
    std::vector<std::pair<int, Scalar> > SetOfT;
    // the feasible set is implicitly given by "SetOfT - {t_i==0}"
    TVector d = -g;
    // n operations
    for (int j = 0; j < DIM; j++) {
      if (g(j) == 0) {
        SetOfT.push_back(std::make_pair(j, std::numeric_limits<Scalar>::max()));
      } else {
        Scalar tmp = 0;
        if (g(j) < 0) {
          tmp = (x(j) - problem.upperBound()(j)) / g(j);
        } else {
          tmp = (x(j) - problem.lowerBound()(j)) / g(j);
        }
        SetOfT.push_back(std::make_pair(j, tmp));
      }
    }
    // sortedindices [1,0,2] means the minimal element is on the 1-st entry
    std::vector<int> sortedIndices = sort_indexes(SetOfT);
    x_cauchy = x;
    // Initialize
    // p :=     W^Scalar*p
    VariableTVector p = (W.transpose() * d);                     // (2mn operations)
    // c :=     0
    c = VariableTVector::Zero(W.cols());
    // f' :=    g^Scalar*d = -d^Td
    Scalar f_prime = -d.dot(d);                         // (n operations)
    // f'' :=   \theta*d^Scalar*d-d^Scalar*W*M*W^Scalar*d = -\theta*f' - p^Scalar*M*p
    Scalar f_doubleprime = (Scalar)(-1.0 * theta) * f_prime - p.dot(M * p); // (O(m^2) operations)
    // \delta t_min :=  -f'/f''
    Scalar dt_min = -f_prime / f_doubleprime;
    // t_old :=     0
    Scalar t_old = 0;
    // b :=     argmin {t_i , t_i >0}
    int i = 0;
    for (int j = 0; j < DIM; j++) {
      i = j;
      if (SetOfT[sortedIndices[j]].second > 0)
        break;
    }
    int b = sortedIndices[i];
    // see below
    // t                    :=  min{t_i : i in F}
    Scalar t = SetOfT[b].second;
    // \delta Scalar             :=  t - 0
    Scalar dt = t ;
    // examination of subsequent segments
    while ((dt_min >= dt) && (i < DIM)) {
      if (d(b) > 0)
        x_cauchy(b) = problem.upperBound()(b);
      else if (d(b) < 0)
        x_cauchy(b) = problem.lowerBound()(b);
      // z_b = x_p^{cp} - x_b
      Scalar zb = x_cauchy(b) - x(b);
      // c   :=  c +\delta t*p
      c += dt * p;
      // cache
      VariableTVector wbt = W.row(b);
      f_prime += dt * f_doubleprime + (Scalar) g(b) * g(b) + (Scalar) theta * g(b) * zb - (Scalar) g(b) *
      wbt.transpose() * (M * c);
      f_doubleprime += (Scalar) - 1.0 * theta * g(b) * g(b)
                       - (Scalar) 2.0 * (g(b) * (wbt.dot(M * p)))
                       - (Scalar) g(b) * g(b) * wbt.transpose() * (M * wbt);
      p += g(b) * wbt.transpose();
      d(b) = 0;
      dt_min = -f_prime / f_doubleprime;
      t_old = t;
      ++i;
      if (i < DIM) {
        b = sortedIndices[i];
        t = SetOfT[b].second;
        dt = t - t_old;
      }
    }
    dt_min = std::max(dt_min, (Scalar)0.0);
    t_old += dt_min;
    #pragma omp parallel for
    for (int ii = i; ii < x_cauchy.rows(); ii++) {
      x_cauchy(sortedIndices[ii]) = x(sortedIndices[ii]) + t_old * d(sortedIndices[ii]);
    }
    c += dt_min * p;
  }
  /**
   * @brief find alpha* = max {a : a <= 1 and  l_i-xc_i <= a*d_i <= u_i-xc_i}
   * @details [long description]
   *
   * @param FreeVariables [description]
   * @return [description]
   */
  Scalar findAlpha(const TProblem &problem, TVector &x_cp, VariableTVector &du, std::vector<int> &FreeVariables) {
    Scalar alphastar = 1;
    const unsigned int n = FreeVariables.size();
    assert(du.rows() == n);
    for (unsigned int i = 0; i < n; i++) {
      if (du(i) > 0) {
        alphastar = std::min(alphastar, (problem.upperBound()(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
      } else {
        alphastar = std::min(alphastar, (problem.lowerBound()(FreeVariables[i]) - x_cp(FreeVariables[i])) / du(i));
      }
    }
    return alphastar;
  }
  /**
   * @brief solving unbounded probelm
   * @details [long description]
   *
   * @param SubspaceMin [description]
   */
  void SubspaceMinimization(const TProblem &problem, TVector &x_cauchy, TVector &x, VariableTVector &c, TVector &g,
  TVector &SubspaceMin) {
    Scalar theta_inverse = 1 / theta;
    std::vector<int> FreeVariablesIndex;
    for (int i = 0; i < x_cauchy.rows(); i++) {
      if ((x_cauchy(i) != problem.upperBound()(i)) && (x_cauchy(i) != problem.lowerBound()(i))) {
        FreeVariablesIndex.push_back(i);
      }
    }
    const int FreeVarCount = FreeVariablesIndex.size();
    MatrixType WZ = MatrixType::Zero(W.cols(), FreeVarCount);
    for (int i = 0; i < FreeVarCount; i++)
      WZ.col(i) = W.row(FreeVariablesIndex[i]);
    TVector rr = (g + theta * (x_cauchy - x) - W * (M * c));
    // r=r(FreeVariables);
    MatrixType r = MatrixType::Zero(FreeVarCount, 1);
    for (int i = 0; i < FreeVarCount; i++)
      r.row(i) = rr.row(FreeVariablesIndex[i]);
    // STEP 2: "v = w^T*Z*r" and STEP 3: "v = M*v"
    VariableTVector v = M * (WZ * r);
    // STEP 4: N = 1/theta*W^T*Z*(W^T*Z)^T
    MatrixType N = theta_inverse * WZ * WZ.transpose();
    // N = I - MN
    N = MatrixType::Identity(N.rows(), N.rows()) - M * N;
    // STEP: 5
    // v = N^{-1}*v
    v = N.lu().solve(v);
    // STEP: 6
    // HERE IS A MISTAKE IN THE ORIGINAL PAPER!
    VariableTVector du = -theta_inverse * r - theta_inverse * theta_inverse * WZ.transpose() * v;
    // STEP: 7
    Scalar alpha_star = findAlpha(problem, x_cauchy, du, FreeVariablesIndex);
    // STEP: 8
    VariableTVector dStar = alpha_star * du;
    SubspaceMin = x_cauchy;
    for (int i = 0; i < FreeVarCount; i++) {
      SubspaceMin(FreeVariablesIndex[i]) = SubspaceMin(FreeVariablesIndex[i]) + dStar(i);
    }
  }
 public:
  LbfgsbSolver() : m_historySize(5) {}

  void setHistorySize(const int hs) { m_historySize = hs; }

  void minimize(TProblem &problem, TVector &x0) {
    DIM = x0.rows();
    theta = 1.0;
    W = MatrixType::Zero(DIM, 0);
    M = MatrixType::Zero(0, 0);
    xHistory.push_back(x0);
    MatrixType yHistory = MatrixType::Zero(DIM, 0);
    MatrixType sHistory = MatrixType::Zero(DIM, 0);
    TVector x = x0, g = x0;
    Scalar f = problem.value(x);
    problem.gradient(x, g);
    // conv. crit.
    this->m_current.reset();
    this->m_status = Status::Continue;
std::cerr << "x.rows(): " << x.rows() << std::endl;
std::cerr << "x.cols(): " << x.cols() << std::endl;
std::cerr << "g.rows(): " << g.rows() << std::endl;
std::cerr << "g.cols(): " << g.cols() << std::endl;
std::cerr << problem.lowerBound().rows() << std::endl;
std::cerr << problem.lowerBound().cols() << std::endl;
std::cerr << problem.upperBound().rows() << std::endl;
std::cerr << problem.upperBound().cols() << std::endl;
    while (
      problem.callback(this->m_current, x) 
      && (((x - g).cwiseMax(problem.lowerBound()).cwiseMin(problem.upperBound()) - x).template lpNorm<Eigen::Infinity>() >= 1e-4) 
      && (this->m_status == Status::Continue)
    ) {
      Scalar f_old = f;
      TVector x_old = x;
      TVector g_old = g;
      // STEP 2: compute the cauchy point
      TVector CauchyPoint = TVector::Zero(DIM);
      VariableTVector c = VariableTVector::Zero(W.cols());
      getGeneralizedCauchyPoint(problem, x, g, CauchyPoint, c);
      // STEP 3: compute a search direction d_k by the primal method for the sub-problem
      TVector SubspaceMin;
      SubspaceMinimization(problem, CauchyPoint, x, c, g, SubspaceMin);
      // STEP 4: perform linesearch and STEP 5: compute gradient
      Scalar alpha_init = 1.0;
      const Scalar rate = MoreThuente<TProblem, 1>::linesearch(x,  SubspaceMin-x ,  problem, alpha_init);
      // update current guess and function information
      x = x - rate*(x-SubspaceMin);
      f = problem.value(x);
      problem.gradient(x, g);
      xHistory.push_back(x);
      // prepare for next iteration
      TVector newY = g - g_old;
      TVector newS = x - x_old;
      // STEP 6:
      Scalar test = newS.dot(newY);
      test = (test < 0) ? -1.0 * test : test;
      if (test > 1e-7 * newY.squaredNorm()) {
        if (yHistory.cols() < m_historySize) {
          yHistory.conservativeResize(DIM, this->m_current.iterations + 1);
          sHistory.conservativeResize(DIM, this->m_current.iterations + 1);
        } else {
          yHistory.leftCols(m_historySize - 1) = yHistory.rightCols(m_historySize - 1).eval();
          sHistory.leftCols(m_historySize - 1) = sHistory.rightCols(m_historySize - 1).eval();
        }
        yHistory.rightCols(1) = newY;
        sHistory.rightCols(1) = newS;
        // STEP 7:
        theta = (Scalar)(newY.transpose() * newY) / (newY.transpose() * newS);
        W = MatrixType::Zero(yHistory.rows(), yHistory.cols() + sHistory.cols());
        W << yHistory, (theta * sHistory);
        MatrixType A = sHistory.transpose() * yHistory;
        MatrixType L = A.template triangularView<Eigen::StrictlyLower>();
        MatrixType MM(A.rows() + L.rows(), A.rows() + L.cols());
        MatrixType D = -1 * A.diagonal().asDiagonal();
        MM << D, L.transpose(), L, ((sHistory.transpose() * sHistory) * theta);
        M = MM.inverse();
      }
      if (fabs(f_old - f) < 1e-8) {
        // successive function values too similar
        break;
      }
      ++this->m_current.iterations;
      this->m_current.gradNorm = g.norm();
      this->m_status = checkConvergence(this->m_stop, this->m_current);
    }
    x0 = x;
    if (this->m_debug > DebugLevel::None) {
        std::cout << "Stop status was: " << boost::lexical_cast<std::string>(this->m_status) << std::endl;
        std::cout << "Stop criteria were: " << std::endl << this->m_stop << std::endl;
        std::cout << "Current values are: " << std::endl << this->m_current << std::endl;
    }
  }
};
#endif /* LBFGSBSOLVER_H_ */
}
