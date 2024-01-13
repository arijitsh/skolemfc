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

#include <errno.h>
#include <string.h>
#include <sys/stat.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
//#include <coz.h>

#include "GitSHA1.h"
#include "counter.h"
#include "time_mem.h"

using std::cout;
using std::endl;
using std::map;
using namespace SkolemFC;

uint64_t Counter::solve(Config _conf)
{
  conf = _conf;
  startTime = cpuTimeTotal();

  openLogFile();
  randomEngine.seed(conf.seed);

  uint64_t solCount = count();

  if (conf.verb)
  {
    cout << "c [sklfc] SkolemFC T: " << (cpuTimeTotal() - startTime) << " s"
         << endl;
  }
  return solCount;
}

uint64_t Counter::count()
{
  // TODO Do actual counting here
  return 42;
}

void Counter::openLogFile()
{
  if (!conf.logfilename.empty())
  {
    logfile.open(conf.logfilename.c_str());
    if (!logfile.is_open())
    {
      cout << "[appmc] Cannot open Counter log file '" << conf.logfilename
           << "' for writing." << endl;
      exit(1);
    }

    logfile << std::left << std::setw(5) << "sampl"
            << " " << std::setw(4) << "iter"
            << " " << std::setw(4) << "hash"
            << " " << std::setw(4) << "full"
            << " " << std::setw(4) << "sols"
            << " " << std::setw(4) << "rep"
            << " " << std::setw(7) << "T"
            << " " << std::setw(7) << "total T" << endl;
  }
}

void Counter::write_log(bool sampling,
                        int iter,
                        uint32_t hashCount,
                        int found_full,
                        uint32_t num_sols,
                        uint32_t repeat_sols,
                        double used_time)
{
  if (!conf.logfilename.empty())
  {
    logfile << std::left << std::setw(5) << (int)sampling << " " << std::setw(4)
            << iter << " " << std::setw(4) << hashCount << " " << std::setw(4)
            << found_full << " " << std::setw(4) << num_sols << " "
            << std::setw(4) << repeat_sols << " " << std::setw(7) << std::fixed
            << std::setprecision(2) << used_time << " " << std::setw(7)
            << std::fixed << std::setprecision(2)
            << (cpuTimeTotal() - startTime) << endl;
  }
}
