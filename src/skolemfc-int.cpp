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
                   const uint32_t _seed,
                   const uint32_t _verbosity)
{
  epsilon = _epsilon;
  delta = _delta;
  verbosity = _verbosity;
  seed = _seed;
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

void SkolemFCInt::SklFCInt::print_formula(const vector<vector<Lit>>& formula)
{
  cout << "c Below is created formula" << endl;
  for (const auto& clause : formula)
  {
    cout << "c ";
    for (const Lit& lit : clause)
    {
      cout << lit << " ";
    }
    cout << endl;
  }
  cout << "c Finished printing G formula" << endl;
}

void SkolemFCInt::SklFCInt::create_g_formula()
{
  g_formula_clauses.clear();
  std::vector<uint32_t> mapped_exists_vars(exists_vars.size());
  // Create a mapping for exists_vars to new variables
  for (size_t i = 0; i < exists_vars.size(); ++i)
  {
    mapped_exists_vars[i] =
        nVars() + i;  // New variable starting from nVars + 1
  }

  // Add F(X, Y') to g_formula_clauses
  for (const auto& clause : clauses)
  {
    g_formula_clauses.push_back(clause);
    new_clause.clear();
    for (const Lit& lit : clause)
    {
      uint32_t var = lit.var();
      auto it = std::find(exists_vars.begin(), exists_vars.end(), var);
      if (it != exists_vars.end())
      {
        // Map Y variable to corresponding Y' variable
        size_t index = std::distance(exists_vars.begin(), it);
        new_clause.push_back(Lit(mapped_exists_vars[index], lit.sign()));
      }
      else
      {
        new_clause.push_back(lit);
      }
    }
    g_formula_clauses.push_back(new_clause);
  }

  // Add (Y ≠ Y') to g_formula_clauses

  diff_clause.clear();
  for (size_t i = 0; i < exists_vars.size(); ++i)
  {
    uint32_t y = exists_vars[i];
    uint32_t y_prime = nVars() + i;
    uint32_t aux_y = nVars() + exists_vars.size()
                     + i;  // Auxiliary variable for each pair y, y'
    if (verbosity > 2)
    {
      cout << "c y: " << y + 1 << ", y': " << y_prime + 1
           << ", aux_y: " << aux_y + 1 << endl;
    }

    // Clauses (y, y', -aux_y) and (-y, -y', -aux_y)
    g_formula_clauses.push_back(
        {Lit(y, false), Lit(y_prime, false), ~Lit(aux_y, false)});
    g_formula_clauses.push_back(
        {~Lit(y, false), ~Lit(y_prime, false), ~Lit(aux_y, false)});

    // Collect aux_y for the final clause
    diff_clause.push_back(Lit(aux_y, false));
  }

  g_formula_clauses.push_back(diff_clause);
  n_g_vars = nVars() + 2 * exists_vars.size();

  cout << "c [sklfc] G formula created with " << g_formula_clauses.size()
       << " clauses and " << nGVars() << " variables." << endl;
  if (verbosity > 2) print_formula(g_formula_clauses);
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
