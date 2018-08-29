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

    const I n = static_cast<I>(ukpi.items.size());
    const double c = ukpi.c;

    unique_ptr<double[]> w(new double[n]);
    unique_ptr<double[]> p(new double[n]);

    for (I i = 0; i < n; ++i) {
      w[i] = static_cast<double>(ukpi.items[i].w);
      p[i] = static_cast<double>(ukpi.items[i].p);
    }

    try {
      GRBEnv env = GRBEnv();
      GRBModel model = GRBModel(env);
      unique_ptr<GRBVar[]> x(model.addVars(n, GRB_INTEGER));
      GRBLinExpr objective_expr, capacity_lhs_expr;

      objective_expr.setTerms(p, x, n);
      model.setObjective(objective_expr, GRB_MAXIMIZE);
      capacity_lhs_expr.setTerms(w, x, n);
      model.addConstr(capacity_lhs_expr, GRB_LESS_EQUAL, c);
      
      model.optimize();

      sol.opt = static_cast<P>(model.get(GRB_DoubleAttr_ObjVal) + HBM_TOLERANCE);
      for (I i = 0; i < n; ++i) {
        double xi = x[i].get(GRB_DoubleAttr_X);
        assert(abs(xi - static_cast<int>(xi + HBM_TOLERANCE)) < HBM_TOLERANCE);
        if (xi >= (1.0 - HBM_TOLERANCE)) sol.used_items.emplace_back(
          ukpi[i], static_cast<I>(xi + HBM_TOLERANCE), i
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

