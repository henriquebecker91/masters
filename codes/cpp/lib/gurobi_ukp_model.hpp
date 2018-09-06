#ifndef HBM_GUROBI_UKP_MODEL_HPP
#define HBM_GUROBI_UKP_MODEL_HPP

#include <memory>
#include <cassert>
#include <cmath>
#include "gurobi_c++.h"
#include "ukp_common.hpp"

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

#define HBM_TOLERANCE 0.00001

namespace hbm {
  template<typename W, typename P, typename I>
  void gurobi_solve_ukp(
    instance_t<W, P> &ukpi,
    solution_t<W, P, I> &sol,
    int argc,
    argv_t argv
  ) {
    using namespace std;

    const int n = static_cast<int>(ukpi.items.size());
    const double c = static_cast<double>(ukpi.c);

    unique_ptr<double[]> w(new double[n]);
    unique_ptr<double[]> p(new double[n]);
    // NOTE: I (Henrique Becker) tried adding the trivial variable bounds in
    // Gurobi 8.0.1 (Thu Sep  6 13:40:08 -03 2018) but this made the code
    // take 0.6s for instance corepb to 97s. Startling how the code could
    // be much slower if I had the bounds from start and never tested
    // without them.
    //unique_ptr<double[]> ub(new double[n]);
    //unique_ptr<char[]> type(new char[n]);

    for (int i = 0; i < n; ++i) {
      w[i] = static_cast<double>(ukpi.items[i].w);
      p[i] = static_cast<double>(ukpi.items[i].p);
      //ub[i] = static_cast<double>(ukpi.c / ukpi.items[i].w);
      //type[i] = GRB_INTEGER;
    }

    try {
      GRBEnv env = GRBEnv();
      GRBModel model = GRBModel(env);
      unique_ptr<GRBVar[]> x(model.addVars(n, GRB_INTEGER));
      //unique_ptr<GRBVar[]> x(model.addVars(nullptr, ub.get(), nullptr, type.get(), 0, n));
      GRBLinExpr objective_expr, capacity_lhs_expr;

      objective_expr.addTerms(p.get(), x.get(), n);
      model.setObjective(objective_expr, GRB_MAXIMIZE);
      capacity_lhs_expr.addTerms(w.get(), x.get(), n);
      model.addConstr(capacity_lhs_expr, GRB_LESS_EQUAL, c);
      
      model.set(GRB_IntParam_Seed, 0);
      model.set(GRB_IntParam_Threads, 1);
      model.set(GRB_IntParam_DisplayInterval, std::numeric_limits<int>::max());
      // This is the relative MIPGap. The absolute MIPGap (MIPGapAbs) is 1e-10
      // and will not be messed with.
      model.set(GRB_DoubleParam_MIPGap, 0.0);
      model.set(GRB_DoubleParam_TimeLimit, 1800.0);
      // From the manual: "By default, the Gurobi MIP solver strikes a balance
      // between finding new feasible solutions and proving that the current
      // solution is optimal. If you are more interested in finding feasible
      // solutions quickly, you can select MIPFocus=1. If you believe the
      // solver is having no trouble finding good quality solutions, and wish
      // to focus more attention on proving optimality, select MIPFocus=2. If
      // the best objective bound is moving very slowly (or not at all), you
      // may want to try MIPFocus=3 to focus on the bound.
      model.set(GRB_IntParam_MIPFocus, 0);
      model.optimize();

      sol.opt = static_cast<P>(model.get(GRB_DoubleAttr_ObjVal) + HBM_TOLERANCE);
      for (int i = 0; i < n; ++i) {
        double xi = x[i].get(GRB_DoubleAttr_X);
        assert(abs(xi - static_cast<int>(xi + HBM_TOLERANCE)) < HBM_TOLERANCE);
        if (xi >= (1.0 - HBM_TOLERANCE)) sol.used_items.emplace_back(
          ukpi.items[i], static_cast<I>(xi + HBM_TOLERANCE), i
        );
      }
    } catch(GRBException e) {
      cout << "Error code = " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch(...) {
      cout << "Exception during optimization" << endl;
    }
  }
}

#endif //HBM_GUROBI_UKP_MODEL_HPP

