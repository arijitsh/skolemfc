/*
 SkolemFC

 Copyright (c) 2024, Arijit Shaw, Brendan Juba and Kuldeep S. Meel. All rights
 reserved.

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
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>
#include <string>

namespace SkolemFC {

struct Config
{
  double epsilon = 0.80;  // Tolerance.  AAAI-2024 paper default
  double delta = 0.8;     // Confidence. AAAI-2024 paper default
  unsigned verb = 0;
  unsigned verb_cls = 0;
  uint32_t seed = 1;
  int simplify = 1;
  std::vector<uint32_t> forall_vars;
  std::vector<uint32_t> exists_vars;
  std::string logfilename = "";
};

}  // namespace SkolemFC

// SKOLEMFC_CONFIG_H
#endif
