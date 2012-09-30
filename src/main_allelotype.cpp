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
#include <stdlib.h>

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>

#include "src/common.h"
#include "src/Genotyper.h"
#include "src/linear.h"
#include "src/NoiseModel.h"
#include "src/ReadContainer.h"
#include "src/runtime_parameters.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

void show_help() {
  const char* help = "\n\nTo train the genotyping noise model " \
    "from a set of aligned reads:\n"                            \
    "allelotype --command train [OPTIONS] --bam <input.bam> "   \
    "--noise_model <noisemodelprefix> --strinfo <strinfo.tab> " \
    " --haploid chrY \n\n"                        \
    "To run str profiling on a set of aligned reads:\n"         \
    "allelotype --command classify [OPTIONS] --bam <input.bam> "  \
    "--noise_model <noisemodelprefix> [--no-rmdup] --strinfo <strinfo.tab> " \
    "--out <output_prefix> \n\n"                     \
    "To run training and classification on a set of aligned reads:\n" \
    "allelotype --command both [OPTIONS] --bam <input.bam> "          \
    "--noise_model <noisemodelprefix> [--no-rmdup] --strinfo <strininfo.tab> "  \
    "--out <output_prefix> --haploid chrY\n\n"                         \
    "To allelotype without using a stutter noise model:\n"            \
    "allelotype simple [OPTIONS] --bam <input.bam> [--no-rmdup] "     \
    "--out <output_prefix> \n\n"                         \
    "Parameter description:\n"                                          \
    "--command [simple|train|classify|both]: specify which of the tasks\n" \
    "                                        described above to perform\n" \
    "--bam <file1.bam,[file2.bam, ...]>:     comma-separated list of bam files\n" \
    "                                        to analyze\n"              \
    "--noise_model <STRING>:                 prefix of files to write (for --command train or \n" \
    "                                        --command both) or read \n" \
    "                                        (--command classify) noise model parameters to.\n" \
    "                                        An example is $PATH_TO_LOBSTR/models/illumina2\n" \
    "--no-rmdup:                 don't remove pcr duplicates before allelotyping\n"     \
    "--min-het-freq <FLOAT>:     minimum frequency to make a heterozygous call\n" \
    "                            (default: 0.25)\n" \
    "--haploid <chrX,[chrY,...]>:            comma-separated list of chromosomes\n" \
    "                                        that should be forced to have homozygous\n" \
    "                                        calls. Specify --haploid all if the organism\n" \
    "                                        is haploid\n" \
    "--sex [U|M|F]:              Gender of sample, M=male, F=female, U=unknown\n" \
    "                            If gender is irrelevant or doesn't make sense in \n" \
    "                            your usage case, specify --sex U.\n"   \
    "                            (Deprecated. this option will be removed.\n" \
    "                             Use --haploid instead).\n" \
    "--strinfo <strinfo.tab>     File containing statistics for each STR, available in the data/\n" \
    "                            directory of the lobSTR download.\n" \
    "--min-genotype-score <float> Minimum posterior probabililty for a genotype\n"\
    "                             to be reported (default 0)\n" \
    "--min-allele-score <float>   Minimum marginal probability for an allele to be\n"\
    "                             reported (default 0) \n" \
    "--include-flank              Include indels in flanking regions when determining\n" \
    "                             length of the STR allele.\n" \
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
    "--max-matedist <INT>:        Filter reads with a mate distance larger than <INT> bp.\n\n"
    "--exclude-partial            Do not report any information about partially\n" \
    "                             spanning reads.\n";
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
    OPT_HAPLOID,
    OPT_EXCLUDE_PARTIAL,
    OPT_STRINFO,
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
    {"include-flank", 0, 0, OPT_NO_INCLUDE_FLANK},
    {"profile", 0, 0, OPT_PROFILE},
    {"min-het-freq", 1, 0, OPT_MIN_HET_FREQ},
    {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
    {"mapq", 1, 0, OPT_MAXMAPQ},
    {"max-matedist", 1, 0, OPT_MAXMATEDIST},
    {"plot", 0, 0, OPT_PLOTINFO},
    {"min-genotype-score", 1, 0, OPT_POSTERIOR},
    {"min-allele-score", 1, 0, OPT_MARGINAL},
    {"haploid",1,0,OPT_HAPLOID},
    {"exclude-partial",0,0,OPT_EXCLUDE_PARTIAL},
    {"strinfo",1,0,OPT_STRINFO},
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
      AddOption("include-flank", "", false, &user_defined_arguments_allelotyper);
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
      if (string(optarg) == "M") {haploid_chroms_string = "chrY,chrX";}
      else if (string(optarg) == "F") {haploid_chroms_string = "chrX";}
      else {haploid_chroms_string = "";}
      sex_set++;
      AddOption("sex", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NOISEMODEL:
      noise_model = string(optarg);
      if ((command == "train" || command == "both") &&
          (fexists(noise_model.c_str()) || (fexists((noise_model+".stuttermodel").c_str())))) {
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
    case OPT_HAPLOID:
      haploid_chroms_string = string(optarg);
      AddOption("haploid",string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_EXCLUDE_PARTIAL:
      exclude_partial++;
      AddOption("exclude-partial", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_STRINFO:
      strinfofile = string(optarg);
      AddOption("strinfo",string(optarg), true, &user_defined_arguments_allelotyper);
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
        || output_prefix.empty()) {
      cerr << "\n\nERROR: Required arguments are missing. " \
        "Please specify a bam file, " \
        "output prefix, and noise model";
      show_help();
      exit(1);
    }
  }
  if (command == "both" || command == "simple") {
    if (bam_files_string.empty() || output_prefix.empty())  {
      cerr << "\n\nERROR: Required arguments are missing. " \
        "Please specify a bam file and " \
        "output prefix";
      show_help();
      exit(1);
    }
  }
  // check that parameters make sense
  if ((command == "train" || command == "both") &&
      haploid_chroms_string.empty()) {
    cerr << "\n\nERROR: Must specify which locus are " \
      "hemizygous using the --haploid option";
    show_help();
    exit(1);
  }

  if (strinfofile.empty()) {
    cerr << "\n\nERROR: Must specify --strinfo file \n";
    show_help();
    exit(1);
  }
}

void test_logr() {
  cerr << "test logr..." << endl;
  problem prob;
  prob.l = 5;
  prob.n = 2;
  prob.bias = -1;
  prob.y = Malloc(double, prob.l);
  prob.x = Malloc(struct feature_node*, prob.l);
  feature_node* x_space;
  x_space = Malloc(struct feature_node, prob.l*prob.n);
  int j = 0;
  for (int i = 0; i < prob.l; i++) {
    prob.y[i] = 1;
    prob.x[i] = &x_space[j];
    x_space[j].index = 1;
    x_space[j].value = 1;
    ++j;
    x_space[j].index = 2;
    x_space[j].value = 2;
    ++j;
  }
  x_space[j++].index = -1;

  parameter param;
  param.solver_type = L2R_LR;
  param.eps = 0.01;
  param.C = 1;
  param.nr_weight = 0;
  param.p = 0.1;
  param.weight_label = NULL;
	param.weight = NULL;
  train(&prob,&param);
}

int main(int argc, char* argv[]) {
  /* do test */
  // test_logr();

  /* parse command line options */
  parse_commandline_options(argc, argv);

  if (my_verbose) cout << "Running allelotyping with command "
                       << command << endl;

  /* Get haploid chromosomes */
  vector<string> haploid_chroms;
  boost::split(haploid_chroms,haploid_chroms_string,boost::is_any_of(","));

  /* initialize noise model */
  NoiseModel nm(strinfofile, haploid_chroms);

  /* Add reads to read container */
  // TODO add strand and mismatch infor for each read
  vector<string>bam_files;
  boost::split(bam_files, bam_files_string, boost::is_any_of(","));
  if (my_verbose) cout << "Adding reads to read container" << endl;
  ReadContainer read_container;
  read_container.AddReadsFromFile(bam_files, exclude_partial);

  /* Perform pcr dup removal if specified */
  if (rmdup) {
    if (my_verbose) cout << "Performing pcr duplicate removal" << endl;
    read_container.RemovePCRDuplicates();
  }

  /* Train/classify */
  if (command == "train" or command == "both") {
    if (my_verbose) cout << "Training noise model..." << endl;
    nm.Train(&read_container);
  } else if (command == "classify") {
    if (!nm.ReadNoiseModelFromFile(noise_model)) {
      errx(1, "Error reading noise file");
    }
  }
  if (command == "classify" or command == "both") {
    if (my_verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, haploid_chroms, false);
    genotyper.Genotype(read_container,
                       output_prefix + ".genotypes.tab");
  }
  if (command == "simple") {
    if (my_verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, haploid_chroms, true);
    genotyper.Genotype(read_container,
                       output_prefix + ".genotypes.tab");
  }
  return 0;
}
