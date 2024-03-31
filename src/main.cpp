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

po::options_description main_options = po::options_description("Main options");
po::options_description help_options;
po::options_description oracle_options =
    po::options_description("Oracle options");
po::options_description hidden_options =
    po::options_description("Hidden options");

po::variables_map vm;
po::positional_options_description p;
double startTime;
uint32_t verbosity = 1;
bool count_unsat_inputs = false;
bool static_samp_est = true;
bool noguarantee = false;
uint32_t use_unisamp_sampling = 1;
uint32_t exactcount_f = 1;
uint32_t exactcount_g = 0;
uint32_t seed = 0;
uint32_t nthreads = 8;
double epsilon = 0.8;
double delta = 0.4;
double g_counter_epsilon = 0.08;
double g_counter_delta = 0.08;
double max_error_logcounter = 0.1;
double epsilon_weightage_fc = 0.6;
double delta_weightage_fc = 0.5;
string logfilename;
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
  std::ostringstream my_epsilon;

  std::ostringstream my_delta;
  std::ostringstream my_g_counter_epsilon;
  std::ostringstream my_max_error_logcounter;
  std::ostringstream my_epsilon_weightage_fc;
  std::ostringstream my_delta_weightage_fc;
  std::ostringstream my_g_counter_delta;

  my_epsilon << std::setprecision(8) << epsilon;
  my_delta << std::setprecision(8) << delta;
  my_g_counter_epsilon << std::setprecision(8) << g_counter_epsilon;
  my_epsilon_weightage_fc << std::setprecision(8) << epsilon_weightage_fc;
  my_delta_weightage_fc << std::setprecision(8) << delta_weightage_fc;
  my_g_counter_delta << std::setprecision(8) << g_counter_delta;
  my_max_error_logcounter << std::setprecision(8) << max_error_logcounter;

  main_options.add_options()("help,h", "Prints help")(
      "verb,v", po::value(&verbosity)->default_value(1), "verbosity")(
      "seed,s", po::value(&seed)->default_value(seed), "Seed")(
      "threads,j",
      po::value(&nthreads)->default_value(1),
      "Number of threads to use")("version", "Print version info")

      ("epsilon,e",
       po::value(&epsilon)->default_value(epsilon, my_epsilon.str()),
       "Tolerance parameter, i.e. how close is the count from the correct "
       "count? Count output is within bounds of (exact_count/(1+e)) < count < "
       "(exact_count*(1+e)). So e=0.8 means we'll output at most 180%% of "
       "exact count and at least 55%% of exact count. Lower value means more "
       "precise.")(
          "delta,d",
          po::value(&delta)->default_value(delta, my_delta.str()),
          "Confidence parameter, i.e. how sure are we of the result? (1-d) = "
          "probability the count is within range as per epsilon parameter. So "
          "d=0.2 means we are 80%% sure the count is within range as specified "
          "by epsilon. The lower, the higher confidence we have in the count.")(
          "log", po::value(&logfilename), "Logs of SkolemFC execution")(
          "count-unsat",
          po::bool_switch(&count_unsat_inputs)
              ->default_value(count_unsat_inputs),
          "Count for those input variables for which there is no output")(
          "no-static-samp",
          po::bool_switch(&static_samp_est)->default_value(static_samp_est),
          "Employ ApproxMC for sample number estimation beforehand");

  help_options.add(main_options);

  hidden_options.add_options()("input", po::value<string>(), "file to read");

  oracle_options.add_options()(
      "no-guarantee",
      po::bool_switch(&noguarantee)->default_value(noguarantee),
      "Run SkolemFC with extreme performance, but no theoretical guarantee")(
      "use-unisamp",
      po::value(&use_unisamp_sampling)->default_value(use_unisamp_sampling),
      "Use UniSamp for high precision sampling")(
      "exact-f",
      po::value(&exactcount_f)->default_value(exactcount_f),
      "Use Exact Counter to count size of set S0 / S1")(
      "exact-g",
      po::value(&exactcount_g)->default_value(exactcount_g),
      "Use Exact Counter to count size of set S2")(
      "epsilon-fc",
      po::value(&epsilon_weightage_fc)
          ->default_value(epsilon_weightage_fc, my_epsilon_weightage_fc.str()),
      "Weightage of error allowed for function counter")(
      "delta-fc",
      po::value(&delta_weightage_fc)
          ->default_value(delta_weightage_fc, my_delta_weightage_fc.str()),
      "Weightage of Tolerance allowed for function counter")(
      "delta-g",
      po::value(&g_counter_delta)
          ->default_value(g_counter_delta, my_g_counter_delta.str()),
      "Tolerance for approximate counter counting size of S2")(
      "epsilon-g",
      po::value(&g_counter_epsilon)
          ->default_value(g_counter_epsilon, my_g_counter_epsilon.str()),
      "Error for approximate counter counting size of S2")(
      "max-error-logcounter",
      po::value(&max_error_logcounter)
          ->default_value(max_error_logcounter, my_max_error_logcounter.str()),
      "Maximum error limit allowed to logcounter");

  help_options.add(oracle_options);
  help_options.add(hidden_options);
}

void add_supported_options(int argc, char** argv)
{
  add_skolemfc_options();
  p.add("input", 1);

  try
  {
    po::store(po::command_line_parser(argc, argv)
                  .options(help_options)
                  .positional(p)
                  .run(),
              vm);
    if (vm.count("help"))
    {
      cout << "Approximate Skolem Function Counter" << endl;

      cout << "skolemfc [options] inputfile" << endl << endl;

      cout << help_options << endl;
      std::exit(0);
    }

    if (vm.count("version"))
    {
      cout << skolemfc->get_version_info();
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
         << "       Unknown option: " << c.what() << endl;
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
  DimacsParser<StreamBuffer<FILE*, FN>, SklFC> parser(
      skolemfc, NULL, verbosity);
#else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, GZ>, SklFC> parser(
      skolemfc, NULL, verbosity);
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

  skolemfc = new SkolemFC::SklFC(epsilon, delta, seed, verbosity);

  cout << "c\nc ---- [ banner ] "
          "------------------------------------------------------------\nc\n";

  cout << "c [sklfc] SkolemFC Version: " << skolemfc->get_version_info()
       << endl;
  if (verbosity >= 3)
  {
    cout << "c [sklfc] compilation environment: "
         << skolemfc->get_compilation_env() << endl;
  }

  cout << "c [sklfc] executed with command line: " << command_line << endl;

  double starTime = cpuTime();
  cout << "c [sklfc] using epsilon: " << epsilon << " delta: " << delta
       << " seed: " << seed << endl;

  // parsing the input
  if (vm.count("input") == 0)
  {
    cout << "ERROR: you must pass a file" << endl;
    exit(-1);
  }
  const string inp = vm["input"].as<string>();
  readInAFile(inp);

  // The ordering of setting oracles are interdependent
  // Please do not change the order

  skolemfc->set_noguarntee_mode(noguarantee);

  skolemfc->set_oracles(use_unisamp_sampling, exactcount_f, exactcount_g);
  skolemfc->set_g_counter_parameters(g_counter_epsilon, g_counter_delta);

  skolemfc->check_ready();
  skolemfc->set_num_threads(nthreads);
  skolemfc->set_parameters();
  skolemfc->set_ignore_unsat(!count_unsat_inputs);
  skolemfc->set_static_samp(static_samp_est);
  skolemfc->set_dklr_parameters(
      epsilon_weightage_fc, delta_weightage_fc, max_error_logcounter);

  skolemfc->count();

  cout << "c\nc ---- [ profiling ] "
          "---------------------------------------------------------\nc\n";

  cout << "c [sklfc] finished T: " << std::setprecision(2) << std::fixed
       << (cpuTime() - starTime) << endl;
  cout << "c [sklfc] iterations: " << skolemfc->get_iteration() << endl;

  delete skolemfc;
  return 0;
}
