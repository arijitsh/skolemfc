/******************************************
 SkolemFC

 Copyright (C) 2024, Arijit Shaw, Brendan Juba, and Kuldeep S. Meel.

 All rights reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
***********************************************/

#include "skolemfc.h"

#include <unigen/unigen.h>

#include <iomanip>

#include "GitSHA1.h"
#include "skolemfc-int.h"
#include "time_mem.h"

using namespace SkolemFCInt;
using CMSat::lbool;
using CMSat::Lit;

struct SkolemFC::SklFCPrivate
{
  SklFCPrivate(SkolemFCInt::SklFCInt* _p) : p(_p) {}
  ~SklFCPrivate() { delete p; }
  SkolemFCInt::SklFCInt* p = NULL;
};

SkolemFC::SklFC::SklFC(const double epsilon,
                       const double delta,
                       const uint32_t seed,
                       const uint32_t verbosity)
{
  skolemfc = new SklFCPrivate(
      new SkolemFCInt::SklFCInt(epsilon, delta, seed, verbosity));
}
SkolemFC::SklFC::~SklFC() { delete skolemfc; }

uint32_t SkolemFC::SklFC::nVars() { return skolemfc->p->nvars; }
void SkolemFC::SklFC::new_vars(uint32_t num) { skolemfc->p->nvars += num; }

bool SkolemFC::SklFC::add_clause(const std::vector<Lit>& cl)
{
  return skolemfc->p->add_clause(cl);
}

bool SkolemFC::SklFC::add_exists_var(uint32_t var)
{
  return skolemfc->p->add_exists_var(var);
}

bool SkolemFC::SklFC::add_forall_var(uint32_t var)
{
  return skolemfc->p->add_forall_var(var);
}

void SkolemFC::SklFC::check_ready() { skolemfc->p->check_ready(); }

void SkolemFC::SklFC::set_constants()
{
  start_time_skolemfc = cpuTime();

  epsilon = skolemfc->p->epsilon;
  delta = skolemfc->p->delta;

  cout << "c [sklfc] running with epsilon: " << epsilon
       << " and delta: " << delta << endl;

  epsilon_f = 0.6 * epsilon;
  delta_f = 0.5 * delta;

  thresh = 4.0 * log(2 / delta_f) * (1 + epsilon_f) / (epsilon_f * epsilon_f);
  thresh *= (double)skolemfc->p->exists_vars.size();

  cout << "c [sklfc] threshold (x |Y|) is set to: " << thresh << endl;
}

void SkolemFC::SklFC::get_est0()
{
  ApproxMC::AppMC appmc;
  appmc.new_vars(skolemfc->p->nVars());
  for (auto& clause : skolemfc->p->clauses)
  {
    appmc.add_clause(clause);
  }
  appmc.set_projection_set(skolemfc->p->forall_vars);
  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc) << "] counting for F formula"
       << endl;
  ApproxMC::SolCount c = appmc.count();
  uint64_t f_cnt = std::pow(2, c.hashCount) * c.cellSolCount;
  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc)
       << "]  F formula has (projected) count: " << f_cnt << endl;

  if (f_cnt == 0)
  {
    cout << "c [sklfc] F is UNSAT. Est1 = 0" << endl;
    okay = false;
  }
  value_est0 = skolemfc->p->exists_vars.size();
  uint64_t pow_n = pow(2, skolemfc->p->forall_vars.size());

  if (pow_n > f_cnt)
    value_est0 *= log(pow_n - f_cnt);
  else
    value_est0 = 0;

  log_skolemcount = value_est0;
  cout << "c [sklfc] Est0 = " << value_est0 << endl;
}

void SkolemFC::SklFC::get_g_count()
{
  skolemfc->p->create_g_formula();

  appmc_g.new_vars(skolemfc->p->nGVars());
  for (auto& clause : skolemfc->p->g_formula_clauses)
  {
    appmc_g.add_clause(clause);
  }
  appmc_g.set_projection_set(skolemfc->p->forall_vars);

  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc) << "] counting for G formula"
       << endl;

  ApproxMC::SolCount c = appmc_g.count();
  count_g_formula = std::pow(2, c.hashCount) * c.cellSolCount;
  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc)
       << "] G formula has (projected) count: " << count_g_formula << endl;
  if (count_g_formula == 0)
  {
    cout << "c [sklfc] G is UNSAT. Est1 = 0" << endl;
    okay = false;
  }
}

void SkolemFC::SklFC::get_sample_num_est()
{
  CMSat::SATSolver cms;
  ApproxMC::AppMC apmc;
  vector<Lit> new_clause;
  cms.new_vars(skolemfc->p->nGVars());
  apmc.new_vars(skolemfc->p->nVars());
  for (auto& clause : skolemfc->p->g_formula_clauses)
  {
    cms.add_clause(clause);
  }
  auto res = cms.solve();
  if (res == CMSat::l_False)
  {
    cout << "c Unsat G" << endl;
  }

  for (auto& clause : skolemfc->p->clauses)
  {
    apmc.add_clause(clause);
  }
  vector<lbool> model = cms.get_model();
  for (auto var : skolemfc->p->forall_vars)
  {
    new_clause.clear();
    if (model[var] == CMSat::l_True)
    {
      new_clause.push_back(Lit(var, false));
    }
    else
    {
      assert(model[var] == CMSat::l_False);
      new_clause.push_back(Lit(var, true));
    }
    apmc.add_clause(new_clause);
  }

  ApproxMC::SolCount c = apmc.count();
  uint64_t f_cnt = std::pow(2, c.hashCount) * c.cellSolCount;
  cout << "c [sklfc] Estimated count from each it: " << f_cnt << endl;
  sample_num_est = 2 * thresh / (c.hashCount + log2(c.cellSolCount));
  cout << "c [sklfc]  [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc)
       << "] approximated number of iterations: " << sample_num_est << endl;
}

void SkolemFC::SklFC::unigen_callback(const vector<int>& solution, void*)
{
  for (uint32_t i = 0; i < solution.size(); i++)
  {
    samples_from_unisamp.push_back(solution);
  }
}

void SkolemFC::SklFC::get_samples()
{
  get_g_count();
  get_sample_num_est();
  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc) << "] starting to get "
       << sample_num_est << " samples" << endl;

  auto ug_appmc = new ApproxMC::AppMC;
  auto unigen = new UniGen::UniG(ug_appmc);
  ug_appmc->set_verbosity(0);
  unigen->set_verbosity(0);
  //   unigen->set_callback(unigen_callback,NULL);
  unigen->set_callback([this](const vector<int>& solution,
                              void*) { this->unigen_callback(solution, NULL); },
                       NULL);

  ug_appmc->new_vars(skolemfc->p->nGVars());
  for (auto& clause : skolemfc->p->g_formula_clauses)
  {
    ug_appmc->add_clause(clause);
  }
  ug_appmc->set_projection_set(skolemfc->p->forall_vars);
  ApproxMC::SolCount c = ug_appmc->count();
  unigen->set_full_sampling_vars(skolemfc->p->forall_vars);
  unigen->sample(&c, sample_num_est);
  cout << "c [sklfc] [" << std::setprecision(2) << std::fixed
       << (cpuTime() - start_time_skolemfc) << "] generated " << sample_num_est
       << " samples" << endl;
}

double SkolemFC::SklFC::get_final_count()
{
  return (thresh / (double)iteration) * (double)count_g_formula;
}

vector<vector<Lit>> SkolemFC::SklFC::create_formula_from_sample(int sample_num)
{
  // TODO if already made samples are not enough, then create a new set
  vector<vector<Lit>> formula = skolemfc->p->clauses;
  vector<Lit> new_clause;
  for (auto int_lit : samples_from_unisamp[sample_num])
  {
    bool isNegated = int_lit < 0;
    uint32_t varIndex =
        isNegated ? -int_lit - 1 : int_lit - 1;  // Convert to 0-based index
    Lit literal = isNegated ? ~Lit(varIndex, false) : Lit(varIndex, false);
    new_clause.clear();
    new_clause.push_back(literal);
    formula.push_back(new_clause);
  }
  if (skolemfc->p->verbosity > 2) skolemfc->p->print_formula(formula);
  return formula;
}

void SkolemFC::SklFC::get_and_add_count_for_a_sample()
{
  assert(iteration < sample_num_est);
  vector<vector<Lit>> sampling_formula = create_formula_from_sample(iteration);
  ApproxMC::AppMC appmc;
  appmc.new_vars(skolemfc->p->nVars());
  for (auto& clause : sampling_formula)
  {
    appmc.add_clause(clause);
  }
  ApproxMC::SolCount c = appmc.count();
  // TODO c.hashCount + 1 seems a bad hack
  double logcount_this_it = c.hashCount + log2(c.cellSolCount);

  iteration++;
  log_skolemcount += logcount_this_it;

  if (skolemfc->p->verbosity > 1)
  {
    cout << "c [sklfc] logcount at iteration " << iteration << ": "
         << logcount_this_it << " log_skolemcount: " << log_skolemcount << endl;
  }
}

double SkolemFC::SklFC::count()
{
  set_constants();

  get_est0();

  if (okay) get_samples();

  while ((log_skolemcount <= thresh) && okay)
  {
    get_and_add_count_for_a_sample();
  }

  if (okay) log_skolemcount = (double)value_est0 + get_final_count();

  return log_skolemcount;
}

const char* SkolemFC::SklFC::get_version_info()
{
  return SkolemFCInt::get_version_sha1();
}
const char* SkolemFC::SklFC::get_compilation_env()
{
  return SkolemFCInt::get_compilation_env();
}

void SkolemFC::SklFC::set_verbosity(uint32_t verb) {}

void SkolemFC::SklFC::set_seed(uint32_t seed) {}
