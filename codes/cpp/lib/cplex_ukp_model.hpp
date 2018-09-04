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

    IloExpr max_profit(env), constraint_capacity(env);
    model.add(IloMaximize(env, IloScalProd(p, x)));
    cout << "before IloScalProd" << endl;
    model.add(IloScalProd(w, x) <= c);

    IloCplex cplex(model);

    // configure cplex solver
    cplex.setOut(env.getNullStream()); // disable output
    cplex.setParam(IloCplex::Param::RandomSeed, 0);
    // The AbsMIPGap is the absolute gap between best solution found and
    // the optimistic guess (upper or lower bound). The default is 1e-06
    // and we prefer to not change it because it can mess the CPLEX performance
    // without giving any precision gain. The MIPGap is a relative value
    // with default 1e-04, so any model with an objective value over 0.01
    // would stop first because of this relative gape than because of
    // the absolute one.
    //cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.0);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);

    // The cplex guarantees that setting the threads to one gives a
    // completely deterministic execution, so the Parallel param does
    // not need to be used to guarantee determinism.
    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setParam(IloCplex::Param::Parallel, 1);

    // Maybe define an internal time limit?
    cplex.setParam(IloCplex::Param::TimeLimit, 1000);
    cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 7*1024);
    // Guarantee CPLEX is using the wall clock time (default).
    cplex.setParam(IloCplex::Param::ClockType, 2);
    // TODO: make some preliminary tests with the parameters below to decide
    // if one of them should be used.
    cplex.setParam(IloCplex::Param::Emphasis::MIP, 0);
    // Values for IloCplex::Param::Emphasis::MIP
    // 0 balance optimality and feasibility
    // 1 feasibility over optimality
    // 2 optimality over feasibility
    // 3 "even greater emphasis is placed on proving optimality"
    // 4 "consider this setting when the FEASIBILITY setting has
    //    difficulty finding solutions of acceptable quality."

    cplex.setParam(IloCplex::Param::Simplex::Display, 0);

    // TODO: check if "numerical precision emphasis" should be set to one,
    // the default is zero that is to not worry much about numerical precision.
    // TODO: "presolve switch" maybe there is no reason to try to simplify
    // the model
    // TODO: check if "CPU mask to bind threads to cores" should be used
    // instead of taskset
    // TODO: check if CPLEX "integrality tolerance" has to be configured too.

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
        ukpi.items[i], static_cast<W>(IloRound(xv[i])), i
      );
    }

    env.end();
  }
}

#endif //HBM_CPLEX_UKP_MODEL_HPP

