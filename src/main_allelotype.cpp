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
#include "src/FastaFileReader.h"
#include "src/Genotyper.h"
#include "src/MSReadRecord.h"
#include "src/linear.h"
#include "src/NoiseModel.h"
#include "src/ReadContainer.h"
#include "src/runtime_parameters.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

// Keep track of reference nucleotide for each STR
map<pair<string,int>, char> ref_nucleotides;

void show_help() {
  const char* help = "\nTo train the genotyping noise model " \
    "from a set of aligned reads:\n"                            \
    "allelotype --command train [OPTIONS] --bam <input.bam> "   \
    "--noise_model <noisemodelprefix> --strinfo <strinfo.tab> " \
    " --haploid chrY --index-prefix $PATH_TO_INDEX/lobSTR_ \n\n" \
    "Training outputs model files: \n" \
    "   <noisemodelprefix>.stepmodel\n" \
    "   <noisemodelprefix>.stuttermodel\n" \
    "   <noisemodelprefix>.stutterproblem\n\n" \
    "To run str profiling on a set of aligned reads:\n" \
    "allelotype --command classify [OPTIONS] --bam <input.bam> "  \
    "--noise_model <noisemodelprefix> [--no-rmdup] --strinfo <strinfo.tab> " \
    "--out <output_prefix> --index-prefix $PATH_TO_INDEX/lobSTR_\n\n" \
    "Classifying outputs the files: \n" \
    "   <output_prefix>.vcf \n" \
    "   <output_prefix>.stats \n\n" \
    "Note: parameters are uploaded to Amazon S3 by default. This for\n" \
    "us see how people are using the tool and to help us continue to improve\n" \
    "lobSTR. To turn this function off, specify --noweb.\n\n" \
    "Parameter description:\n" \
    "--command [train|classify]:     (REQUIRED) specify which of the tasks\n" \
    "                                described above to perform\n" \
    "--bam <f1.bam,[f2.bam, ...]>:   (REQUIRED) comma-separated list\n" \
    "                                of bam files to analyze\n" \
    "--strinfo <strinfo.tab>:        (REQUIRED)\n" \
    "                                File containing statistics for each STR,\n" \
    "                                available in the data/ directory of the\n" \
    "                                lobSTR download.\n" \
    "--noise_model <STRING>:         (REQUIRED)\n" \
    "                                prefix of files to write (--command train)\n" \
    "                                or read (--command classify) noise model\n" \
    "                                parameters to.\n" \
    "                                An example is $PATH_TO_LOBSTR/models/illumina2\n" \
    "--index-prefix <STRING>         (REQUIRED) prefix for lobSTR's bwa reference\n" \
    "                                (must be same as for lobSTR alignment)\n" \
    "--no-rmdup:                     don't remove pcr duplicates before allelotyping\n" \
    "--min-het-freq <FLOAT>:         minimum frequency to make a heterozygous call\n" \
    "                                (default: NULL)\n" \
    "--haploid <chrX,[chrY,...]>:    comma-separated list of chromosomes\n" \
    "                                that should be forced to have homozygous\n" \
    "                                calls. Specify --haploid all if the organism\n" \
    "                                is haploid\n" \
    "--include-flank:                Include indels in flanking regions when\n" \
    "                                determining length of the STR allele.\n" \
    "-h,--help:                      display this message\n" \
    "-v,--verbose:                   print out helpful progress messages\n" \
    "--version:                      print out allelotype program version number\n\n" \
    "Options for calculating and reporting allelotypes:\n" \
    "--use-known-alleles <FILE>:     Use prior information about alleles in\n" \
    "                                the population when determining the allelotype.\n" \
    "                                NOTE: if you want to merge VCF files downstream,\n" \
    "                                they must have the same set of ALT alleles listed,\n" \
    "                                and so you must specify that here.\n" \
    "                                This is a tab delimited file with columns: \n" \
    "                                     1. chrom\n" \
    "                                     2. start\n" \
    "                                     3. end\n" \
    "                                     4. <allele1>:<count>;<allele2>:<count>...\n" \
    "                                where counts can be total counts or frequencies.\n" \
    "                                An example file based on 1000 genomes data is\n" \
    "                                available in the download at\n" \
    "                                <PATH-TO-LOBSTR>/data/lobstr.allelepriors.hg19.min50.tab.\n" \
    "                                If this option is specified, only loci and alleles\n" \
    "                                included in <FILE> will be considered. Allele counts\n" \
    "                                are not used unless --generate-posteriors is specified,\n" \
    "                                and can be set to dummy values if not used.\n" \
    "--generate-posteriors           Generate posterior probabilities for\n" \
    "                                possible allelotypes.\n" \
    "                                Requires allele frequencies to be set with \n" \
    "                                --use-known-alleles\n" \
    "Options for reporting allelotypes:\n" \
    "--tab                          Generate <output_prefix>.genotypes.tab file\n" \
    "--sample <STRING>:             Name of sample. Default to --out\n" \
    "--exclude-pos <FILE>:          File of \"chrom\\tpos\" of positions to exclude.\n" \
    "                               For downstream analysis, it is beneficial to exclude\n" \
    "                               any STRs with the same starting point but different motifs,\n" \
    "                               as this will cause errors when using vcftools.\n\n" \
    "Options for filtering reads:\n" \
    "If not specified, no filters applied\n" \
    "--chrom <STRING>:              only look at reads from this chromosome\n" \
    "--max-diff-ref <INT>:          filter reads differing from the\n" \
    "                               reference allele by more than <INT> bp.\n" \
    "--unit:                        filter reads differing by a non-integer\n" \
    "                               number of repeat copies from reference\n" \
    "--mapq <INT>:                  filter reads with mapq scores of more than\n" \
    "                               <INT>.\n" \
    "--max-matedist <INT>:          Filter reads with a mate distance larger than <INT> bp.\n\n"
    "--exclude-partial:             Do not report any information about partially\n" \
    "                               spanning reads.\n\n"
    "Additional options\n" \
    "--noweb                        Do not report any user information and paramters to Amazon S3.\n";
  cout << help;
  exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc, char* argv[]) {
  enum LONG_OPTIONS {
    OPT_BAM,
    OPT_CHROM,
    OPT_COMMAND,
    OPT_DEBUG,
    OPT_EXCLUDE_PARTIAL,
    OPT_EXCLUDE_POS,
    OPT_GENERATE_POSTERIORS,
    OPT_HAPLOID,
    OPT_HELP,
    OPT_INCLUDE_FLANK,
    OPT_INDEX,
    OPT_MAX_DIFF_REF,
    OPT_MAXMAPQ,
    OPT_MAXMATEDIST,
    OPT_MIN_HET_FREQ,
    OPT_NOISEMODEL,
    OPT_NORMDUP,
    OPT_NOWEB,
    OPT_OUTPUT,
    OPT_PRINT_READS,
    OPT_PROFILE,
    OPT_SAMPLE,
    OPT_STRINFO,
    OPT_TAB,
    OPT_UNIT,
    OPT_USE_KNOWN_ALLELES,
    OPT_VERBOSE,
    OPT_VERSION,
  };

  int ch;
  int option_index = 0;

  static struct option long_options[] = {
    {"bam", 1, 0, OPT_BAM},
    {"chrom", 1, 0, OPT_CHROM},
    {"command", 1, 0, OPT_COMMAND},
    {"debug", 0, 0, OPT_DEBUG},
    {"exclude-partial", 0, 0, OPT_EXCLUDE_PARTIAL},
    {"exclude-pos", 1, 0, OPT_EXCLUDE_POS},
    {"generate-posteriors", 0, 0, OPT_GENERATE_POSTERIORS},
    {"haploid", 1, 0, OPT_HAPLOID},
    {"help", 1, 0, OPT_HELP},
    {"include-flank", 0, 0, OPT_INCLUDE_FLANK},
    {"index-prefix", 1, 0, OPT_INDEX},
    {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
    {"mapq", 1, 0, OPT_MAXMAPQ},
    {"max-matedist", 1, 0, OPT_MAXMATEDIST},
    {"min-het-freq", 1, 0, OPT_MIN_HET_FREQ},
    {"noise_model", 1, 0, OPT_NOISEMODEL},
    {"no-rmdup", 0, 0, OPT_NORMDUP},
    {"no-web", 0, 0, OPT_NOWEB},
    {"out", 1, 0, OPT_OUTPUT},
    {"profile", 0, 0, OPT_PROFILE},
    {"reads", 0, 0, OPT_PRINT_READS},
    {"sample", 1, 0, OPT_SAMPLE},
    {"strinfo", 1, 0, OPT_STRINFO},
    {"tab", 0, 0, OPT_TAB},
    {"unit", 0, 0, OPT_UNIT},
    {"use-known-alleles", 1, 0, OPT_USE_KNOWN_ALLELES},
    {"verbose", 0, 0, OPT_VERBOSE},
    {"version", 0, 0, OPT_VERSION},
    {NULL, no_argument, NULL, 0},
  };
  program = ALLELOTYPE;
  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_BAM:
      bam_files_string = string(optarg);
      AddOption("bam", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_CHROM:
      use_chrom = string(optarg);
      AddOption("chrom", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_COMMAND:
      command = string(optarg);
      AddOption("command", string(optarg), true, &user_defined_arguments_allelotyper);
      if ((command != "train") & (command != "classify")) {
        cerr << "\n\n[allelotype-" << _GIT_VERSION << "]"
             <<  "ERROR: Command " << command
             << " is invalid. Command must be one of: train, classify";
        show_help();
        exit(1);
      }
      break;
    case OPT_DEBUG:
      debug = true;
      break;
    case OPT_EXCLUDE_PARTIAL:
      exclude_partial++;
      AddOption("exclude-partial", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_EXCLUDE_POS:
      exclude_positions_file = string(optarg);
      AddOption("exclude-pos", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_GENERATE_POSTERIORS:
      generate_posteriors++;
      AddOption("generate-posteriors", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_HAPLOID:
      haploid_chroms_string = string(optarg);
      AddOption("haploid",string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case 'h':
    case OPT_HELP:
      show_help();
      exit(1);
      break;
    case OPT_INCLUDE_FLANK:
      include_flank = false;
      AddOption("include-flank", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_INDEX:
      index_prefix = string(optarg);
      AddOption("index-prefix", string(optarg), true, &user_defined_arguments_allelotyper);
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
    case OPT_MIN_HET_FREQ:
      min_het_freq = atof(optarg);
      AddOption("min-het-freq", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NOISEMODEL:
      noise_model = string(optarg);
      if ((command == "train") &&
          (fexists(noise_model.c_str()) || (fexists((noise_model+".stuttermodel").c_str())))) {
        errx(1, "Cannot write to specified noise model file. " \
             "This file already exists.");
      }
      AddOption("noise-model", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_NORMDUP:
      rmdup = false;
      AddOption("normdup", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_NOWEB:
      noweb++;
      AddOption("noweb", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_OUTPUT:
      output_prefix = string(optarg);
      AddOption("out", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_PRINT_READS:
      print_reads++;
      AddOption("print-reads", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_PROFILE:
      profile++;
      break;
    case OPT_SAMPLE:
      sample = string(optarg);
      AddOption("sample", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_STRINFO:
      strinfofile = string(optarg);
      AddOption("strinfo", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case OPT_TAB:
      generate_tab++;
      AddOption("tab", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_UNIT:
      unit = true;
      AddOption("unit", "", false, &user_defined_arguments_allelotyper);
      break;
    case OPT_USE_KNOWN_ALLELES:
      known_alleles_file = string(optarg);
      AddOption("use-known-alleles", string(optarg), true, &user_defined_arguments_allelotyper);
      break;
    case 'v':
    case OPT_VERBOSE:
      my_verbose++;
      break;
    case OPT_VERSION:
      cerr << _GIT_VERSION << endl;
      exit(0);
      break;
    case '?':
      show_help();
      exit(1);
    default:
      show_help();
      exit(1);
    }
    ch = getopt_long(argc, argv, "hv?",
                     long_options, &option_index);
  }

  // any arguments left over are extra
  if (optind < argc) {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Unnecessary leftover arguments";
    show_help();
    exit(1);
  }
  if (command.empty()) {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Must specify a command";
    show_help();
    exit(1);
  }
  if (command == "train") {
    if (bam_files_string.empty() || noise_model.empty()) {
      cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Required arguments are missing. " \
        "Please specify a bam file and an output prefix";
      show_help();
      exit(1);
    }
  }
  if (command == "classify") {
    if (bam_files_string.empty() || noise_model.empty()
        || output_prefix.empty()) {
      cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Required arguments are missing. " \
        "Please specify a bam file, " \
        "output prefix, and noise model";
      show_help();
      exit(1);
    }
  }
  // check that parameters make sense
  if ((command == "train") &&
      (haploid_chroms_string.empty())) {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Must specify which locus are " \
      "hemizygous using the --haploid option";
    show_help();
    exit(1);
  }
  if (strinfofile.empty()) {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Must specify --strinfo to use noise model";
    show_help();
    exit(1);
  }
  if (index_prefix.empty() && command != "train") {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: Must specify --index-prefix";
    show_help();
    exit(1);
  }
  // set sample if not specified
  if (sample.empty()) sample = output_prefix;
  // can't generate posteriors without priors
  if (generate_posteriors && known_alleles_file.empty()) {
    cerr << "\n\n[allelotype-" << _GIT_VERSION << "] ERROR: cannot generate priors without "
      "known alleles file";
    show_help();
    exit(1);
  }
}

int main(int argc, char* argv[]) {
  /* parse command line options */
  parse_commandline_options(argc, argv);

  if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Running allelotyping with command "
                       << command << endl;

  /* Get haploid chromosomes */
  vector<string> haploid_chroms;
  boost::split(haploid_chroms,haploid_chroms_string,boost::is_any_of(","));

  /* initialize noise model */
  NoiseModel nm(strinfofile, haploid_chroms);

  /* Load ref character for each STR */
  if (command != "train") {
    FastaFileReader faReader(index_prefix+"ref.fa");
    MSReadRecord ref_record;
    while (faReader.GetNextRead(&ref_record)) {
      vector<string> items;
      string refstring = ref_record.ID;
      split(refstring, '$', items);
      if (items.size() == 7) {
        int start = atoi(items.at(2).c_str())+extend;
        string chrom = items.at(1);
        char refnuc = items.at(6).at(0);
        pair<string, int> locus = pair<string,int>(chrom, start);
        ref_nucleotides.insert(pair< pair<string, int>, char >(locus, refnuc));
      }
    }
  }

  /* Add reads to read container */
  vector<string>bam_files;
  boost::split(bam_files, bam_files_string, boost::is_any_of(","));
  if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Adding reads to read container" << endl;
  ReadContainer read_container;
  read_container.AddReadsFromFile(bam_files, exclude_partial);

  /* Perform pcr dup removal if specified */
  if (rmdup) {
    if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Performing pcr duplicate removal" << endl;
    read_container.RemovePCRDuplicates();
  }

  /* Train/classify */
  if (command == "train") {
    if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Training noise model..." << endl;
    nm.Train(&read_container);
  } else if (command == "classify") {
    if (!nm.ReadNoiseModelFromFile(noise_model)) {
      errx(1, "[allelotype-%s] ERROR: Problem reading noise file", _GIT_VERSION);
    }
  }
  if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Loading reference genome..." << endl;
  if (command == "classify") {
    Genotyper genotyper(&nm, haploid_chroms, &ref_nucleotides);
    // If available, load priors
    if (!known_alleles_file.empty()) {
      if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Loading known alleles..." << endl;
      genotyper.LoadPriors(known_alleles_file);
    }
    if (my_verbose) cerr << "[allelotype-" << _GIT_VERSION << "] Classifying allelotypes..." << endl;
    if (generate_tab) {
      genotyper.Genotype(read_container,
                         output_prefix + ".genotypes.tab",
                         output_prefix + ".vcf");
    } else {
      genotyper.Genotype(read_container,
                         "/dev/null",
                         output_prefix + ".vcf");
    }
  }
  return 0;
}
