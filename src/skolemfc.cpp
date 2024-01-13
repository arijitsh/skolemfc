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

#include "GitSHA1.h"

using namespace SkolemFC;
using CMSat::Lit;

struct SkolemFC::SklPrivateData
{
  SklPrivateData(SkolemFC::SklFC* _p) : p(_p) {}
  ~SklPrivateData() { delete p; }
  SkolemFC::SklFC* p = NULL;
};

SkolemFC::SklFC::SklFC(const double epsilon,
                       const double delta,
                       const uint32_t seed,
                       const uint32_t verbosity)
{
  skolemfc =
      new SklPrivateData(new SkolemFC::SklFC(epsilon, delta, seed, verbosity));
}
SkolemFC::SklFC::~SklFC() { delete skolemfc; }

uint32_t SkolemFC::SklFC::nVars() { return nvars; }
void SkolemFC::SklFC::new_vars(uint32_t num) { nvars += num; }

bool SkolemFC::SklFC::add_clause(const std::vector<Lit>& cl)
{
  return skolemfc->p->add_clause(cl);
}

const char* SkolemFC::SklFC::get_version_info()
{
  return SkolemFC::get_version_sha1();
}
const char* SkolemFC::SklFC::get_compilation_env()
{
  return SkolemFC::get_compilation_env();
}

void SklFC::set_verbosity(uint32_t verb) {}

void SklFC::set_seed(uint32_t seed) {}
