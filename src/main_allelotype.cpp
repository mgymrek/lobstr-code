/*
Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

This file is part of lobSTR.

lobSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lobSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <err.h>
#include <error.h>
#include <getopt.h>
#include <RInside.h>
#include <stdlib.h>

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>

#include "src/common.h"
#include "src/Genotyper.h"
#include "src/NoiseModel.h"
#include "src/ReadContainer.h"
#include "src/runtime_parameters.h"

using namespace std;

void show_help() {
  const char* help = "\n\nTo train the genotyping noise model " \
    "from a set of aligned reads:\n"                            \
    "allelotype --command train [OPTIONS] --bam <input.bam> "   \
    "--noise_model <noisemodel.txt> --sex M \n\n"                        \
    "To run str profiling on a set of aligned reads:\n"         \
    "allelotype --command classify [OPTIONS] --bam <input.bam> "  \
    "--noise_model <noisemodel.txt> [--no-rmdup] "                \
    "--out <output_prefix> --sex [M|F|U]\n\n"                     \
    "To run training and classification on a set of aligned reads:\n" \
    "allelotype --command both [OPTIONS] --bam <input.bam> "          \
    "--noise_model <noisemodel.txt> [--no-rmdup] "                    \
    "--out <output_prefix> --sex [M|F|U]\n\n"                         \
    "To allelotype without using a stutter noise model:\n"            \
    "allelotype simple [OPTIONS] --bam <input.bam> [--no-rmdup] "     \
    "--out <output_prefix> --sex [M|F|U]\n\n"                         \
    "Parameter description:\n"                                          \
    "--command [simple|train|classify|both]: specify which of the tasks\n" \
    "                                        described above to perform\n" \
    "--bam <file1.bam,[file2.bam, ...]>:     comma-separated list of bam files\n" \
    "                                        to analyze\n"              \
    "--noise_model <STRING>:                 file to write (for --command train or \n" \
    "                                        --command both) or read \n" \
    "                                        (--command classify) noise model parameters to.\n" \
    "--no-rmdup:                 don't remove pcr duplicates before allelotyping\n"     \
    "--min-het-freq <FLOAT>:     minimum frequency to make a heterozygous call\n" \
    "                            (default: 0.2)\n" \
    "--sex [U|M|F]:              Gender of sample, M=male, F=female, U=unknown\n"    \
    "                            If gender is irrelevant or doesn't make sense in \n " \
    "                            your usage case, specify --sex U.\n" \
    "-h: display this message\n"                                        \
    "-v: print out helpful progress messages\n\n" \
    "Options for filtering reads:\n" \
    "If not specified, no filters applied\n" \
    "--max-diff-ref <INT>:        filter reads differing from the\n"   \
    "                             reference allele by more than <INT> bp.\n" \
    "--unit:                      filter reads differing by a non-integer\n" \
    "                             number of repeat copies from reference\n" \
    "--mapq <INT>:                filter reads with mapq scores of more than\n" \
    "                             <INT>.\n"                             \
    "--max-matedist <INT>:        Filter reads with a mate distance larger than <INT> bp.\n\n";
  cout << help;
  exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc, char* argv[]) {
  enum LONG_OPTIONS {
    OPT_COMMAND,
    OPT_BAM,
    OPT_OUTPUT,
    OPT_NORMDUP,
    OPT_SEX,
    OPT_NOISEMODEL,
    OPT_UNIT,
    OPT_HELP,
    OPT_VERBOSE,
    OPT_DEBUG,
    OPT_NO_INCLUDE_FLANK,
    OPT_PROFILE,
    OPT_PRINT_READS,
    OPT_MIN_HET_FREQ,
    OPT_MAX_DIFF_REF,
    OPT_MAXMAPQ,
    OPT_MAXMATEDIST,
    OPT_PLOTINFO,
    OPT_POSTERIOR,
    OPT_MARGINAL,
  };

  int ch;
  int option_index = 0;

  static struct option long_options[] = {
    {"reads", 0, 0, OPT_PRINT_READS},
    {"bam", 1, 0, OPT_BAM},
    {"command", 1, 0, OPT_COMMAND},
    {"out", 1, 0, OPT_OUTPUT},
    {"no-rmdup", 0, 0, OPT_NORMDUP},
    {"sex", 1, 0, OPT_SEX},
    {"noise_model", 1, 0, OPT_NOISEMODEL},
    {"unit", 0, 0, OPT_UNIT},
    {"help", 1, 0, OPT_HELP},
    {"debug", 0, 0, OPT_DEBUG},
    {"no-include-flank", 0, 0, OPT_NO_INCLUDE_FLANK},
    {"profile", 0, 0, OPT_PROFILE},
    {"min-het-freq", 1, 0, OPT_MIN_HET_FREQ},
    {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
    {"mapq", 1, 0, OPT_MAXMAPQ},
    {"max-matedist", 1, 0, OPT_MAXMATEDIST},
    {"plot", 0, 0, OPT_PLOTINFO},
    {"min-genotype-score", 1, 0, OPT_POSTERIOR},
    {"min-allele-score", 1, 0, OPT_MARGINAL},
    {NULL, no_argument, NULL, 0},
  };

  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_MARGINAL:
      MIN_POSTERIOR = atof(optarg);
      AddOption("min-genotype-score",string(optarg), true,
                &user_defined_arguments_allelotyper);
      break;
    case OPT_POSTERIOR:
      MIN_POSTERIOR = atof(optarg);
      AddOption("min-allele-score",string(optarg), true,
                &user_defined_arguments_allelotyper);
      break;
    case OPT_PLOTINFO:
      plot_info++;
      break;
    case OPT_PRINT_READS:
      print_reads++;
      AddOption("print-reads", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_BAM:
      bam_files_string = string(optarg);
      AddOption("bam", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NO_INCLUDE_FLANK:
      include_flank = false;
      AddOption("no-include-flank", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_PROFILE:
      profile++;
      break;
    case OPT_MIN_HET_FREQ:
      min_het_freq = atof(optarg);
      AddOption("min-het-freq", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_COMMAND:
      command = string(optarg);
      AddOption("command", string(optarg), true, &user_defined_arguments_allelotyper);
      if ((command != "train") & (command != "classify") &
          (command != "both") & (command != "simple")) {
        cerr << "\n\nERROR: Command " << command
             << " is invalid. Command must be one of: train, " \
          "classify, both, simple";
        show_help();
        exit(1);
      }
      break;
    case OPT_OUTPUT:
      output_prefix = string(optarg);
      AddOption("out", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NORMDUP:
      rmdup = false;
      AddOption("normdup", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_SEX:
      if (string(optarg) != "F" &&
          string(optarg) != "M" &&
          string(optarg) != "U") {
        errx(1, "--sex must be F,M, or U");
      }
      if (string(optarg) == "F") male = false;
      if (string(optarg) == "U") sex_unknown = true;
      sex_set++;
      AddOption("sex", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NOISEMODEL:
      noise_model = string(optarg);
      if ((command == "train" || command == "both") &&
          fexists(noise_model.c_str())) {
        errx(1, "Cannot write to specified noise model file. " \
             "This file already exists.");
      }
      AddOption("noise-model", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_UNIT:
      unit = true;
      user_defined_arguments_allelotyper += "unit=True;";
      AddOption("unit", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_MAX_DIFF_REF:
      max_diff_ref = atoi(optarg);
      AddOption("max-diff-ref", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_MAXMAPQ:
      max_mapq = atoi(optarg);
      AddOption("mapq", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_MAXMATEDIST:
      max_matedist = atoi(optarg);
      AddOption("max-matedist", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case 'v':
    case OPT_VERBOSE:
      my_verbose++;
      break;
    case 'h':
    case OPT_HELP:
      show_help();
      exit(1);
      break;
    case OPT_DEBUG:
      debug = true;
      break;
    case '?':
      show_help();
      exit(1);
      break;
    default:
      show_help();
      exit(1);
    }
    ch = getopt_long(argc, argv, "hv?",
                     long_options, &option_index);
  }

  // any arguments left over are extra
  if (optind < argc) {
    cerr << "\n\nERROR: Unnecessary leftover arguments";
    show_help();
    exit(1);
  }
  // check that we have the mandatory parameters
  if (command.empty()) {
    cerr << "\n\nERROR: Must specify a command";
    show_help();
    exit(1);
  }
  if (command == "train") {
    if (bam_files_string.empty() || noise_model.empty()) {
      cerr << "\n\nERROR: Required arguments are missing. " \
        "Please specify a bam file and an output prefix";
      show_help();
      exit(1);
    }
    male = true;
  }
  if (command == "classify") {
    if (bam_files_string.empty() || noise_model.empty()
        || output_prefix.empty() || !sex_set) {
      cerr << "\n\nERROR: Required arguments are missing. " \
        "Please specify a bam file, " \
        "output prefix, noise model, and gender";
      show_help();
      exit(1);
    }
  }
  if (command == "both" || command == "simple") {
    if (bam_files_string.empty() || output_prefix.empty() || !sex_set)  {
      cerr << "\n\nERROR: Required arguments are missing. " \
        "Please specify a bam file, " \
        "output prefix, and gender";
      show_help();
      exit(1);
    }
  }
  // check that parameters make sense
  if ((command == "train" || command == "both") && !male) {
    cerr << "\n\nERROR: Cannot train on female sample";
    show_help();
    exit(1);
  }
}

/* copied from common.h */
/*bool fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
  }*/

int main(int argc, char* argv[]) {
  /* parse command line options */
  parse_commandline_options(argc, argv);

  if (my_verbose) cout << "Running allelotyping with command "
                       << command << endl;

  /* initialize noise model */
  RInside R(argc, argv);
  NoiseModel nm(&R);

  /* Add reads to read container */
  vector<string>bam_files;
  boost::split(bam_files, bam_files_string, boost::is_any_of(","));
  if (my_verbose) cout << "Adding reads to read container" << endl;
  ReadContainer read_container;
  read_container.AddReadsFromFile(bam_files);

  /* Perform pcr dup removal if specified */
  if (rmdup) {
    if (my_verbose) cout << "Performing pcr duplicate removal" << endl;
    read_container.RemovePCRDuplicates();
  }

  /* Train/classify */
  if (command == "train" or command == "both") {
    if (my_verbose) cout << "Training noise model..." << endl;
    nm.Train(&read_container);
    nm.WriteNoiseModelToFile(noise_model);
  } else if (command == "classify") {
    if (!nm.ReadNoiseModelFromFile(noise_model)) {
      errx(1, "Error reading noise file");
    }
  }
  if (command == "classify" or command == "both") {
    if (my_verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, male, false);
    genotyper.Genotype(read_container,
                       output_prefix + ".genotypes.tab");
  }
  if (command == "simple") {
    if (my_verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, male, true);
    genotyper.Genotype(read_container,
                       output_prefix + ".genotypes.tab");
  }
  return 0;
}
