// Code based on the cutstock.cpp sample of the CPLEX licensed materials (property of IBM); and changed by Henrique Becker (by the year of 2016).

#define CPLEX       1
#define UKP5_FP     2
#define UKP5_FP_NS  3
#define UKP5_INT    4
#define UKP5_INT_NS 5
#define MTU1        6
#define MTU2        7
#define MGREENDP    8
#define MGREENDP1   9

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <chrono>
#include <cstdint>

#if KNAPSACK_SOLVER == UKP5_FP
#include <vector>
#include <numeric>

#include "ukp5.hpp"
#endif

#define WEIGHT uint_fast32_t
#define INT_P  uint_fast64_t
// IloNum is another name for double.
#define FP_P   IloNum
#define IX_TYPE uint_fast32_t

#ifndef KNAPSACK_SOLVER
  #error KNAPSACK_SOLVER macro not defined
#endif

ILOSTLBEGIN

static void readData (const char* filename, IloNum& rollWidth,
               IloNumArray& size, IloNumArray& amount);
static void report1 (IloCplex& cutSolver, IloNumVarArray Cut, IloRangeArray Fill);

using std::cout;
using std::endl;
using namespace std::chrono;

/// MAIN PROGRAM ///

int main(int argc, char **argv)
{
  IloEnv env;

  string name = "undefined";
  #if KNAPSACK_SOLVER == CPLEX
  name = "cplex";
  #elif KNAPSACK_SOLVER == UKP5_FP
  name = "ukp5_fp";
  #else
    #error VALUE FOR MACRO KNAPSACK_SOLVER WAS NOT EXPECTED
  #endif
  cout << "algorithm_name: " << name << endl;

  IloNum      rollWidth;
  IloNumArray amount(env);
  IloNumArray size(env);

  if (argc > 1) {
    readData(argv[1], rollWidth, size, amount);
  } else {
    /* Example Input file:
    115
    [25, 40, 50, 55, 70]
    [50, 36, 24, 8, 30]
    * Where the first line has the master roll width, the second line
    * a list of the available sizes, and the third line, the respective
    * demand for each size on the second line.
    */
    cout << "usage: " << argv[0] << " <filename in csp format>" << endl;
    return EXIT_FAILURE;
  }

  /// MASTER PROBLEM
  IloModel cutOpt(env);

  IloObjective   RollsUsed = IloAdd(cutOpt, IloMinimize(env));
  IloRangeArray  Fill = IloAdd(cutOpt, IloRangeArray(env, amount, rollWidth));
  IloNumVarArray Cut(env);

  IloInt nWdth = size.getSize();

  // Create initial set of columns/patterns. Each initial pattern
  // is comprised of only one size, cut as many times as possible
  // over a roll of size rollWidth.
  for (IloInt j = 0; j < nWdth; j++) {
    int max_qt_j_in_roll = static_cast<int>(rollWidth / size[j]);
    Cut.add(IloNumVar(RollsUsed(1) + Fill[j](max_qt_j_in_roll)));
  }

  IloCplex cutSolver(cutOpt);
  // configure cplex solver
  cutSolver.setOut(env.getNullStream()); // disable output
  cutSolver.setParam(IloCplex::RandomSeed, 42);
  cutSolver.setParam(IloCplex::EpGap,       0);
  cutSolver.setParam(IloCplex::Threads,     1);

  /// PRICING PROBLEM (PATTERN-GENERATION/KNAPSACK)

  #if KNAPSACK_SOLVER == CPLEX
  IloModel patGen (env);

  IloObjective ReducedCost = IloAdd(patGen, IloMaximize(env));
  IloNumVarArray Use(env, nWdth, 0.0, IloInfinity, ILOINT);
  patGen.add(IloScalProd(size, Use) <= rollWidth);

  IloCplex patSolver(patGen);
  IloNumArray price(env, nWdth);
  // configure cplex solver
  patSolver.setOut(env.getNullStream()); // disable output
  patSolver.setParam(IloCplex::RandomSeed, 42);
  patSolver.setParam(IloCplex::EpGap,       0);
  patSolver.setParam(IloCplex::Threads,     1);
  #elif KNAPSACK_SOLVER == UKP5_FP
  hbm::instance_t<size_t, IloNum> ukpi;
  ukpi.c = static_cast<size_t>(rollWidth);
  ukpi.items.resize(nWdth);
  std::vector<size_t> idx(nWdth, 0);
  hbm::solution_t<size_t, IloNum, size_t> sol;
  hbm::ukp5_conf_t<size_t> conf;
  #endif

  /// COLUMN-GENERATION PROCEDURE ///
  size_t num_iter = 0;
  bool last_iter  = false;
  steady_clock::time_point before, after;
  double curr_knapsack_time = 0, all_knapsacks_time = 0;
  double curr_master_prob_time= 0, all_master_prob_time = 0;
  IloNumArray newPatt(env, nWdth);

  while (!last_iter) {
    /// OPTIMIZE OVER CURRENT PATTERNS ///
    before = steady_clock::now();
    cutSolver.solve();
    after = steady_clock::now();
    curr_master_prob_time = duration_cast<duration<double>>(after - before).count();
    all_master_prob_time += curr_master_prob_time;

    /// FIND AND ADD A NEW PATTERN ///
    
    // Initialize knapsack instance
    for (IloInt i = 0; i < nWdth; i++) {
      #if KNAPSACK_SOLVER == CPLEX
      price[i] = cutSolver.getDual(Fill[i]);
      #elif KNAPSACK_SOLVER == UKP5_FP
      ukpi.items[i].w = static_cast<size_t>(size[i]);
      ukpi.items[i].p = cutSolver.getDual(Fill[i]);
      #endif
    }

    #if KNAPSACK_SOLVER == CPLEX
    ReducedCost.setLinearCoefs(Use, price);
    #endif

    // WE COUNT THE TIME USED TO SOLVE THE KNAPSACK AND RECOVER THE SOLUTION
    before = steady_clock::now();
    #if KNAPSACK_SOLVER == UKP5_FP
    std::iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(),
      [&ukpi](size_t a, size_t b) { return ukpi.items[a] < ukpi.items[b]; });

    hbm::ukp5(ukpi, sol, conf);
    for (IloInt i = 0; i < nWdth; i++) newPatt[i] = 0;
    for (auto &it : sol.used_items) {
      newPatt[idx[it.ix]] = static_cast<IloNum>(it.qt);
    }
    sol.used_items.clear();
    #elif KNAPSACK_SOLVER == CPLEX
    patSolver.solve();
    patSolver.getValues(newPatt, Use);
    #endif
    after = steady_clock::now();
    curr_knapsack_time = duration_cast<duration<double>>(after - before).count();
    all_knapsacks_time += curr_knapsack_time;

    double sigma = ldexp(1, -40); // == 1*pow(2, -40) ~= pow(10, 12)

    cout << "num_iter: " << ++num_iter << endl;
    #if KNAPSACK_SOLVER == CPLEX
    if (patSolver.getValue(ReducedCost) <= 1.0 + sigma) last_iter = true;
    cout << "hex_opt: " << hexfloat << patSolver.getValue(ReducedCost) << endl;
    cout << "dec_opt: " << defaultfloat << patSolver.getValue(ReducedCost) << endl;
    #elif KNAPSACK_SOLVER == UKP5_FP
    if (sol.opt <= 1.0 + sigma) last_iter = true;

    cout << "hex_opt: " << hexfloat << sol.opt << endl;
    cout << "dec_opt: " << defaultfloat << sol.opt << endl;
    #endif
    cout << "hex_knapsack_time: " << hexfloat << curr_knapsack_time << endl;
    cout << "dec_knapsack_time: " << defaultfloat << curr_knapsack_time << endl;

    for (IloInt i = 0; i < nWdth; ++i) {
      if (newPatt[i] > 0) {
        cout << "ix: " << i << " qt: " << static_cast<size_t>(newPatt[i]) << " w: " << static_cast<size_t>(size[i]) << " hex_p: " << hexfloat << cutSolver.getDual(Fill[i]) << " dec_p: " << defaultfloat << cutSolver.getDual(Fill[i]) << endl;
      }
    }
    cout << endl;

    if (last_iter) {
      cout << "hex_sum_master_prob_time: " << hexfloat << all_master_prob_time << endl;
      cout << "dec_sum_master_prob_time: " << defaultfloat << all_master_prob_time << endl;
      cout << endl;
      cout << "hex_sum_knapsack_time: " << hexfloat << all_knapsacks_time << endl;
      cout << "dec_sum_knapsack_time: " << defaultfloat << all_knapsacks_time << endl;
      cout << endl;
      cout << "total_iter: " << num_iter << endl;
      report1(cutSolver, Cut, Fill);
    } else {
      Cut.add(IloNumVar(RollsUsed(1) + Fill(newPatt)));
    }
  }

  env.end();

  return EXIT_SUCCESS;
}


static void readData (const char* filename, IloNum& rollWidth,
               IloNumArray& size, IloNumArray& amount)
{
  ifstream in(filename);
  if (in) {
    in >> rollWidth;
    in >> size;
    in >> amount;
  }
  else {
    cout << "No such file: " << filename << endl;
    throw(1);
  }
}

static void report1 (IloCplex& cutSolver, IloNumVarArray Cut,
              IloRangeArray Fill)
{
  cout << "dec_rolls: " << defaultfloat << cutSolver.getObjValue() << endl;
  cout << "hex_rolls: " << hexfloat << cutSolver.getObjValue() << endl;
  for (IloInt j = 0; j < Cut.getSize(); j++) {
    if (cutSolver.getValue(Cut[j]) > 0) {
      cout << defaultfloat << "dec cut " << j << " " << cutSolver.getValue(Cut[j]) << endl;
    }
  }
  for (IloInt j = 0; j < Cut.getSize(); j++) {
    if (cutSolver.getValue(Cut[j]) > 0) {
      cout << hexfloat << "hex cut " << j << " " << cutSolver.getValue(Cut[j]) << endl;
    }
  }
  cout << endl;
  for (IloInt i = 0; i < Fill.getSize(); i++) {
    if (cutSolver.getDual(Fill[i])) {
      cout << defaultfloat << "dec fill " << i << " " << cutSolver.getDual(Fill[i]) << endl;
    }
  }
  for (IloInt i = 0; i < Fill.getSize(); i++) {
    if (cutSolver.getDual(Fill[i])) {
      cout << hexfloat << "hex fill " << i << " " << cutSolver.getDual(Fill[i]) << endl;
    }
  }
}

