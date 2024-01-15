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

#ifndef SKOLEMFC_H__
#define SKOLEMFC_H__

#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#ifdef CMS_LOCAL_BUILD
#include "cryptominisat.h"
#else
#include <cryptominisat5/cryptominisat.h>
#endif

#include <approxmc/approxmc.h>

using std::vector;

namespace SkolemFC {

struct SklFCPrivate;

struct SklFC
{
 public:
  SklFC(const double epsilon = 0.5,
        const double delta = 0.5,
        const uint32_t seed = 1,
        const uint32_t verbosity = 0);
  ~SklFC();
  static const char* get_version_info();
  static const char* get_compilation_env();

  // Adding CNF
  uint32_t nVars();
  void new_var();
  void new_vars(uint32_t num);
  bool add_clause(const std::vector<CMSat::Lit>& lits);
  bool add_forall_var(uint32_t var);
  bool add_exists_var(uint32_t var);

  void check_ready();

  void set_constants();
  void get_est0();
  void get_g_count();
  void get_samples();
  void get_and_add_count_for_a_sample();
  double get_final_count();
  void get_sample_num_est();

  uint64_t count();

  // Set config
  void set_seed(uint32_t seed);
  void set_verbosity(uint32_t verb);

 private:
  SklFCPrivate* skolemfc = NULL;
  double log_skolemcount;
  double sum_logcount;
  double thresh;
  uint64_t sample_num_est;
  ApproxMC::AppMC appmc_g;
  void unigen_callback(const vector<int>& solution, void*);
};

}  // namespace SkolemFC

#endif
