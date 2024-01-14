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

#include <approxmc/approxmc.h>

#include "GitSHA1.h"
#include "skolemfc-int.h"

using namespace SkolemFCInt;
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

uint64_t SkolemFC::SklFC::count()
{
  uint64_t logcount = 0;

  ApproxMC::AppMC appmc;
  appmc.new_vars(skolemfc->p->nVars());
  for (auto& clause : skolemfc->p->clauses)
  {
    appmc.add_clause(clause);
  }
  appmc.set_projection_set(skolemfc->p->forall_vars);
  ApproxMC::SolCount c = appmc.count();
  uint64_t f_cnt = std::pow(2, c.hashCount) * c.cellSolCount;
  cout << "c [sklfc] F formula has (projected) count: " << f_cnt << endl;

  skolemfc->p->create_g_formula();

  ApproxMC::AppMC appmc_g;
  appmc_g.new_vars(skolemfc->p->nGVars());
  for (auto& clause : skolemfc->p->g_formula_clauses)
  {
    appmc_g.add_clause(clause);
  }
  appmc_g.set_projection_set(skolemfc->p->forall_vars);
  c = appmc_g.count();
  uint64_t g_cnt = std::pow(2, c.hashCount) * c.cellSolCount;
  cout << "c [sklfc] G formula has (projected) count: " << g_cnt << endl;

  return logcount;
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
