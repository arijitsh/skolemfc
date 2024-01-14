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

#include "skolemfc-int.h"

#include <math.h>
#include <string.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include "GitSHA1.h"
#include "time_mem.h"

using std::cout;
using std::endl;

using namespace SkolemFCInt;

SklFCInt::SklFCInt(const double _epsilon,
                   const double _delta,
                   const uint32_t seed,
                   const uint32_t _verbosity)
{
  epsilon = _epsilon;
  delta = _delta;
  verbosity = _verbosity;
}

SklFCInt::~SklFCInt() {}

const char* SklFCInt::get_version_info() const
{
  return SkolemFCInt::get_version_sha1();
}

const char* SklFCInt::get_compilation_env() const
{
  return SkolemFCInt::get_compilation_env();
}

bool SklFCInt::add_clause(const vector<Lit>& cl)
{
  clauses.push_back(cl);
  return false;
}

bool SkolemFCInt::SklFCInt::add_exists_var(uint32_t e_var)
{
  exists_vars.push_back(e_var);
  return true;
}

void printFormula(const vector<vector<Lit>>& formula) {
  cout << "c Below is G formula" << endl;
    for (const auto& clause : formula) {
      cout << "c ";
        for (const Lit& lit : clause) {
           cout << lit << " ";
        }
       cout << endl;
    }
  cout << "c Finished printing G formula" << endl;


}


bool SkolemFCInt::SklFCInt::create_g_formula()
{
  // Add F(X, Y') to g_formula_clauses
  uint32_t numvars = nVars();
  cout << "c nvars = " << numvars << endl;
  for (const auto& clause : clauses)
  {
    g_formula_clauses.push_back(clause);
    newClause.clear();
    for (const Lit& lit : clause)
    {
      uint32_t var = lit.var();
      // Check if the variable is in exists_vars (Y)
      if (std::find(exists_vars.begin(), exists_vars.end(), var)
          != exists_vars.end())
      {
        // Replace Y variable with Y' variable (Y' = Y + offset)
        newClause.push_back(Lit(var + exists_vars.size(), lit.sign()));
      }
      else
      {
        newClause.push_back(lit);
      }
    }
    g_formula_clauses.push_back(newClause);
  }

  // Add (Y ≠ Y') to g_formula_clauses
  for (uint32_t y : exists_vars)
  {
    diffClause.clear();
    diffClause.push_back(Lit(y, false));           // Y
    diffClause.push_back(~Lit(y + nvars, false));  // Not Y'
    g_formula_clauses.push_back(diffClause);

    diffClause.clear();
    diffClause.push_back(~Lit(y, false));         // Not Y
    diffClause.push_back(Lit(y + nvars, false));  // Y'
    g_formula_clauses.push_back(diffClause);
  }
  cout << "c [sklfc] G formula created with " << g_formula_clauses.size()
       << " clauses" << endl;
  if (verbosity > 2 ) printFormula(g_formula_clauses);
}


bool SkolemFCInt::SklFCInt::add_forall_var(uint32_t a_var)
{
  forall_vars.push_back(a_var);
  return true;
}

void SklFCInt::check_ready() const
{
  cout << "c solver got clauses: " << clauses.size()
       << " e vars: " << exists_vars.size() << " a vars: " << forall_vars.size()
       << endl;
}
