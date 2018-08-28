#ifndef HBM_GUROBI_UKP_MODEL_HPP
#define HBM_GUROBI_UKP_MODEL_HPP

#include <memory> // for unique_ptr
#include "gurobi_c++.h"
#include "ukp_common.hpp"

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

namespace hbm {
  template<typename W, typename P, typename I>
  void gurobi_solve_ukp(
    instance_t<W, P> &ukpi,
    solution_t<W, P, I> &sol,
    int argc,
    argv_t argv
  ) {
    using namespace std;

    try {
      GRBEnv env = GRBEnv();
      GRBModel model = GRBModel(env);
      // Creates capacity constraint with only the first variable.
      const w0 = static_cast<double>(ukpi.items[0].w);
      const p0 = static_cast<double>(ukpi.items[0].p);
      GRBConstr capacity = model.addConstr(
        w0 * model.addVar(0.0, ukpi.c/w0, p0, GRB_INTEGER),
        GRB_LESS_EQUAL,
        static_cast<double>(ukpi.c)
      );
      // CREATE THE FIRST VARIABLE FIRST, AND THEN OBJECTIVE, AND THEN CONSTRAINT,
      // AN THEN ALL OTHER VARIABLES LINKED TO OBJECTIVE AND CONSTRAINT
      model.setObjective(, GRB_MAXIMIZE);
      
      //std::unique_ptr< GRBVar[] > vars = new GRBVar[ukpi.n];
      double wi, pi;
      for (I i = 1; i < n; ++i) {
        wi = static_cast<double>(ukpi.items[i].w);
        pi = static_cast<double>(ukpi.items[i].p);
        static_cast<void>(
          model.addVar(0.0, ukpi.c/wi, pi, GRB_INTEGER, 1, &capacity, &wi)
        );
      }

      model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

      // Add constraint: x + 2 y + 3 z <= 4

      model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

      // Add constraint: x + y >= 1

      model.addConstr(x + y >= 1, "c1");

      // Optimize model

      model.optimize();

      cout << x.get(GRB_StringAttr_VarName) << " "
           << x.get(GRB_DoubleAttr_X) << endl;
      cout << y.get(GRB_StringAttr_VarName) << " "
           << y.get(GRB_DoubleAttr_X) << endl;
      cout << z.get(GRB_StringAttr_VarName) << " "
           << z.get(GRB_DoubleAttr_X) << endl;

      cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    } catch(GRBException e) {
      cout << "Error code = " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch(...) {
      cout << "Exception during optimization" << endl;
    }
  }
}

#endif //HBM_GUROBI_UKP_MODEL_HPP

