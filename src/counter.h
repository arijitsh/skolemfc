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

#ifndef COUNTER_H_
#define COUNTER_H_

#include <cstdint>
#include <fstream>
#include <map>
#include <mutex>
#include <random>

#include "config.h"
#include "constants.h"
#include "skolemfc.h"

using std::map;
using std::string;
using std::vector;
using namespace CMSat;

namespace SkolemFC {

class Counter
{
 public:
  uint64_t solve(Config _conf);
  string get_version_info() const;

 private:
  Config conf;
  uint64_t count();
  void add_appmc_options();

  ////////////////
  // Helper functions
  ////////////////

  void write_log(bool sampling,
                 int iter,
                 uint32_t hashCount,
                 int found_full,
                 uint32_t num_sols,
                 uint32_t repeat_sols,
                 double used_time);
  void openLogFile();
  void call_after_parse();
  void ban_one(const uint32_t act_var, const vector<lbool>& model);

  void readInAFile(SATSolver* solver2, const string& filename);
  void readInStandardInput(SATSolver* solver2);
  // Data so we can output temporary count when catching the signal
  vector<uint64_t> numHashList;
  vector<int64_t> numCountList;

  ////////////////
  // internal data
  ////////////////
  double startTime;
  std::ofstream logfile;
  std::mt19937 randomEngine;
  uint32_t orig_num_vars;
  double total_inter_simp_time = 0;
  uint32_t threshold;  // precision, it's computed
  uint32_t cnf_dump_no = 0;

  int argc;
  char** argv;
};

}  // namespace SkolemFC

#endif  // COUNTER_H_
