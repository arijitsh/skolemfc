/*
 SkolemFC

 Copyright (c) 2024, Arijit Shaw, Brendan Juba and Kuldeep S. Meel.

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
 */

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

#include <signal.h>

#include <atomic>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef USE_ZLIB
#include <zlib.h>
#endif

#include <dimacsparser.h>

#include "config.h"
#include "skolemfc.h"
#include "time_mem.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace CMSat;
using namespace SkolemFC;

po::options_description skolemfc_options =
    po::options_description("skolemfc options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;
double startTime;
SkolemFC::Config conf;
SkolemFC::SklFC* skolemfc = NULL;
string elimtofile;
string recover_file;

int recompute_sampling_set = 0;
uint32_t orig_sampling_set_size = 0;
uint32_t polar_mode = 0;
int sparsify = true;
int renumber = true;
bool gates = true;

// static void signal_handler(int) {
//     cout << endl << "c [skolemfc] INTERRUPTING ***" << endl << std::flush;
//     common.interrupt_asap = true;
// }

void add_skolemfc_options()
{
  conf.verb = 1;

  skolemfc_options.add_options()("help,h", "produce help message")(
      "absolute,a", "get the absolute count, no approximation")(
      "count,c", "count from already created sampling files")(
      "sample,s", "create sampling files only")("exact,x",
                                                "use in exact counting mode")(
      "delta,d",
      po::value<double>()->default_value(0.8),
      "value for delta [0 .. 1] (0.8)")("epsilon,e",
                                        po::value<double>()->default_value(0.8),
                                        "value for epsilon [0 .. 1] (0.8)")(
      "limit,l",
      po::value<int>()->default_value(0),
      "limit the number of sampling files [integer]")(
      "timeout,t",
      po::value<int>()->default_value(36000),
      "limit the timeout in seconds [integer]")(
      "filename", po::value<std::string>(), "input file");

  po::options_description debug_options("Debug options");
  debug_options.add_options()("exact,x", "use in exact counting mode");

  help_options.add(skolemfc_options);
  help_options.add(debug_options);
}

void add_supported_options(int argc, char** argv)
{
  add_skolemfc_options();
  p.add("input", -1);

  try
  {
    po::store(po::command_line_parser(argc, argv)
                  .options(help_options)
                  .positional(p)
                  .run(),
              vm);
    if (vm.count("help"))
    {
      cout << "Minimal projection set finder and simplifier." << endl
           << endl
           << "skolemfc [options] inputfile [outputfile]" << endl;

      cout << help_options << endl;
      std::exit(0);
    }

    if (vm.count("version"))
    {
      cout << "c [skolemfc] SHA revision: " << skolemfc->get_version_info()
           << endl;
      cout << "c [skolemfc] Compilation environment: "
           << skolemfc->get_compilation_env() << endl;
      std::exit(0);
    }

    po::notify(vm);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<po::unknown_option> >& c)
  {
    cerr << "ERROR: Some option you gave was wrong. Please give '--help' to "
            "get help"
         << endl
         << "       Unkown option: " << c.what() << endl;
    std::exit(-1);
  }
  catch (boost::bad_any_cast& e)
  {
    std::cerr << "ERROR! You probably gave a wrong argument type" << endl
              << "       Bad cast: " << e.what() << endl;

    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<
             po::invalid_option_value> >& what)
  {
    cerr << "ERROR: Invalid value '" << what.what() << "'" << endl
         << "       given to option '" << what.get_option_name() << "'" << endl;

    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<
             po::multiple_occurrences> >& what)
  {
    cerr << "ERROR: " << what.what() << " of option '" << what.get_option_name()
         << "'" << endl;

    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<po::required_option> >&
             what)
  {
    cerr << "ERROR: You forgot to give a required option '"
         << what.get_option_name() << "'" << endl;

    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<
             po::too_many_positional_options_error> >& what)
  {
    cerr << "ERROR: You gave too many positional arguments. Only the input CNF "
            "can be given as a positional option."
         << endl;
    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<po::ambiguous_option> >&
             what)
  {
    cerr << "ERROR: The option you gave was not fully written and matches"
         << endl
         << "       more than one option. Please give the full option name."
         << endl
         << "       The option you gave: '" << what.get_option_name() << "'"
         << endl
         << "       The alternatives are: ";
    for (size_t i = 0; i < what.alternatives().size(); i++)
    {
      cout << what.alternatives()[i];
      if (i + 1 < what.alternatives().size())
      {
        cout << ", ";
      }
    }
    cout << endl;

    std::exit(-1);
  }
  catch (boost::exception_detail::clone_impl<
         boost::exception_detail::error_info_injector<
             po::invalid_command_line_syntax> >& what)
  {
    cerr << "ERROR: The option you gave is missing the argument or the" << endl
         << "       argument is given with space between the equal sign."
         << endl
         << "       detailed error message: " << what.what() << endl;
    std::exit(-1);
  }
}

void readInAFile(const string& filename)
{
#ifndef USE_ZLIB
  FILE* in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, FN>, SklFC> parser(skolemfc, NULL, 0);
#else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, GZ>, SklFC> parser(skolemfc, NULL, 0);
#endif

  if (in == NULL)
  {
    std::cerr << "ERROR! Could not open file '" << filename
              << "' for reading: " << strerror(errno) << endl;

    std::exit(-1);
  }

  if (!parser.parse_DIMACS(in, true))
  {
    exit(-1);
  }

#ifndef USE_ZLIB
  fclose(in);
#else
  gzclose(in);
#endif
}

void set_config(SkolemFC::SklFC* skl)
{
  cout << "c [skolemfc] using seed: " << conf.seed << endl;
  skl->set_verbosity(conf.verb);
  skl->set_seed(conf.seed);
}

int main(int argc, char** argv)
{
// Die on division by zero etc.
#if defined(__GNUC__) && defined(__linux__)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  // Reconstruct the command line so we can emit it later if needed
  string command_line;
  for (int i = 0; i < argc; i++)
  {
    command_line += string(argv[i]);
    if (i + 1 < argc)
    {
      command_line += " ";
    }
  }
  add_supported_options(argc, argv);
  skolemfc =
      new SkolemFC::SklFC(conf.epsilon, conf.delta, conf.seed, conf.verb);

  cout << "c [sklfc] SklFC Version: " << skolemfc->get_version_info() << endl;
  if (conf.verb >= 2)
  {
    cout << "c [sklfc] compilation environment: "
         << skolemfc->get_compilation_env() << endl;

    cout << "c [sklfc] executed with command line: " << command_line << endl;
  }

  double starTime = cpuTime();
  cout << "c [sklfc] using epsilon: " << conf.epsilon
       << " delta: " << conf.delta << " seed: " << conf.seed << endl;

  // parsing the input
  if (vm.count("input") == 0)
  {
    cout << "ERROR: you must pass a file" << endl;
    exit(-1);
  }
  const string inp = vm["input"].as<string>();
  readInAFile(inp);

  cout << "c [sklfc] finished T: " << std::setprecision(2) << std::fixed
       << (cpuTime() - starTime) << endl;

  delete skolemfc;
  return 0;
}
