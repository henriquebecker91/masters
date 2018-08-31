#ifndef HBM_CPLEX_UKP_MODEL_HPP
#define HBM_CPLEX_UKP_MODEL_HPP

#include <cassert>
#include <cmath>
#include <ilcplex/ilocplex.h>
#include "ukp_common.hpp"

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << (#var ": ") << var << std::endl
#endif

#define HBM_TOLERANCE 0.00001

namespace hbm {
  template<typename W, typename P, typename I>
  void cplex_solve_ukp(
    instance_t<W, P> &ukpi,
    solution_t<W, P, I> &sol,
    int argc,
    argv_t argv
  ) {
    using namespace std;
    auto &out = cout;

    const IloInt n = static_cast<IloInt>(ukpi.items.size());
    const IloInt c = static_cast<IloInt>(ukpi.c);

    HBM_PRINT_VAR(n);
    HBM_PRINT_VAR(c);

    IloEnv env;

    IloIntArray w(env, n);
    IloIntArray p(env, n);

    for (IloInt i = 0; i < n; ++i) {
      w[i] = static_cast<IloInt>(ukpi.items[i].w);
      p[i] = static_cast<IloInt>(ukpi.items[i].p);
    }

    IloIntVarArray x(env, n);

    IloModel model(env);
    cout << "before maximize" << endl;
    model.add(IloMaximize(env, IloScalProd(p, x)));
    cout << "before IloScalProd" << endl;
    model.add(IloScalProd(w, x) <= c);

    IloCplex cplex(model);

    // configure cplex solver
    cplex.setOut(env.getNullStream()); // disable output
    cplex.setParam(IloCplex::RandomSeed, 0);
    cplex.setParam(IloCplex::EpGap,      0);
    cplex.setParam(IloCplex::Threads,    1);

    cout << "before solve" << endl;
    cplex.solve();

    // The variables are from an IloIntVarArray, and even so their values need
    // to be captured in IloNumArray variable (i.e., a floating point variable)
    // and IloRound be used to get the real integer values.
    IloNumArray xv(env, n);
    cplex.getValues(xv, x);

    sol.opt = static_cast<P>(IloRound(cplex.getObjValue()));
    for (int i = 0; i < n; ++i) {
      if (IloRound(xv[i]) >= 1) sol.used_items.emplace_back(
        ukpi.items[i], static_cast<W>(xv[i]), i
      );
    }

    //model.end();
    //cplex.end();
    //x.end();
    //w.end();
    //p.end();
    env.end();
  }
}

#endif //HBM_CPLEX_UKP_MODEL_HPP

