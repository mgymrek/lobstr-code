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
#include <getopt.h>
#include <stdlib.h>

#include <boost/algorithm/string.hpp>
#include <ctime>
#include <iostream>
#include <string>

#include "src/common.h"
#include "src/FastaFileReader.h"
#include "src/Genotyper.h"
#include "src/MSReadRecord.h"
#include "src/linear.h"
#include "src/NoiseModel.h"
#include "src/ReadContainer.h"
#include "src/ReferenceSTR.h"
#include "src/runtime_parameters.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

// Keep track of reference nucleotide for each STR
map<pair<string,int>, string> ref_nucleotides;
// Keep track of reference repseq for each STR
map<pair<string,int>, string> ref_repseq;

void show_help() {
  const char* help = "\nTo train the genotyping noise model " \
    "from a set of aligned reads:\n"                            \
    "allelotype --command train [OPTIONS] --bam <input.bam> "   \
    "--noise_model <noisemodelprefix> --strinfo <strinfo.tab> " \
    " --haploid chrY\n\n" \
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
    "   <output_prefix>.allelotype.stats \n\n" \
    "Note: parameters are uploaded to Amazon S3 by default. This for\n" \
    "us see how people are using the tool and to help us continue to improve\n" \
    "lobSTR. To turn this function off, specify --noweb.\n\n" \
    "Parameter description:\n" \
    "--command [train|classify]:     (REQUIRED) specify which of the tasks\n" \
    "                                described above to perform\n" \
    "--bam <f1.bam,[f2.bam, ...]>:   (REQUIRED) comma-separated list\n" \
    "                                of bam files to analyze. Each sample should have\n \
    "                                "a unique read group.\n" \
    "--strinfo <strinfo.tab>:        (REQUIRED)\n" \
    "                                File containing statistics for each STR,\n" \
    "                                available in the data/ directory of the\n" \
    "                                lobSTR download.\n" \
    "--noise_model <STRING>:         (REQUIRED)\n" \
    "                                prefix of files to write (--command train)\n" \
    "                                or read (--command classify) noise model\n" \
    "                                parameters to.\n" \
    "                                An example is $PATH_TO_LOBSTR/models/illumina2\n" \
    "--index-prefix <STRING>         (REQUIRED for --command classify) prefix for lobSTR's bwa reference\n" \
    "                                (must be same as for lobSTR alignment)\n" \
    "--no-rmdup:                     don't remove pcr duplicates before allelotyping\n" \
    "--min-het-freq <FLOAT>:         minimum frequency to make a heterozygous call\n" \
    "                                (default: NULL)\n" \
    "--haploid <chrX,[chrY,...]>:    comma-separated list of chromosomes\n" \
    "                                that should be forced to have homozygous\n" \
    "                                calls. Specify --haploid all if the organism\n" \
    "                                is haploid. Will be applied to all samples.\n" \
    "--include-flank:                Include indels in flanking regions when\n" \
    "                                determining length of the STR allele.\n" \
    "-h,--help:                      display this message\n" \
    "-v,--verbose:                   print out helpful progress messages\n" \
    "--version:                      print out allelotype program version number\n\n" \
    "Options for calculating and reporting allelotypes:\n" \
    "--annotation <vcf file>         VCF file for STR set annotations (e.g. marshfield_hg19.vcf)\n" \
    "                                For more than one annotation, use comma-separated list of files\n" \
    "Options for reporting allelotypes:\n" \
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
    "--noweb                        Do not report any user information and parameters to Amazon S3.\n";
  cerr << help;
  exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc, char* argv[]) {
  enum LONG_OPTIONS {
    OPT_ANNOTATION,
    OPT_BAM,
    OPT_CHROM,
    OPT_COMMAND,
    OPT_DEBUG,
    OPT_EXCLUDE_PARTIAL,
    OPT_EXCLUDE_POS,
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
    OPT_UNIT,
    OPT_VERBOSE,
    OPT_VERSION,
  };

  int ch;
  int option_index = 0;

  static struct option long_options[] = {
    {"annotation", 1, 0, OPT_ANNOTATION},
    {"bam", 1, 0, OPT_BAM},
    {"chrom", 1, 0, OPT_CHROM},
    {"command", 1, 0, OPT_COMMAND},
    {"debug", 0, 0, OPT_DEBUG},
    {"exclude-partial", 0, 0, OPT_EXCLUDE_PARTIAL},
    {"exclude-pos", 1, 0, OPT_EXCLUDE_POS},
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
    {"unit", 0, 0, OPT_UNIT},
    {"verbose", 0, 0, OPT_VERBOSE},
    {"version", 0, 0, OPT_VERSION},
    {NULL, no_argument, NULL, 0},
  };
  program = ALLELOTYPE;
  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_ANNOTATION:
      annotation_files_string = string(optarg);
      AddOption("annotation", string(optarg), true, &user_defined_arguments_allelotyper);
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
        PrintMessageDieOnError("Command " + command + " is invalid. Must be train or classify", ERROR);
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
        PrintMessageDieOnError("Cannot write to specified noise model file. " \
                               "This file already exists.", ERROR);
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
    case OPT_UNIT:
      unit = true;
      AddOption("unit", "", false, &user_defined_arguments_allelotyper);
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
    PrintMessageDieOnError("Unnecessary leftover arguments", ERROR);
  }
  if (command.empty()) {
    PrintMessageDieOnError("Must specify a command", ERROR);
  }
  if (command == "train") {
    if (bam_files_string.empty() || noise_model.empty()) {
      PrintMessageDieOnError("Please specify a bam file and an output prefix", ERROR);
    }
  }
  if (command == "classify") {
    if (bam_files_string.empty() || noise_model.empty()
        || output_prefix.empty()) {
      PrintMessageDieOnError("Please specify a bam file, output prefix, and noise model", ERROR);
    }
  }
  // check that parameters make sense
  if ((command == "train") &&
      (haploid_chroms_string.empty())) {
    PrintMessageDieOnError("Must specify which loci are hemizygous using --haploid", ERROR);
  }
  if (strinfofile.empty()) {
    PrintMessageDieOnError("Must specify --strinfo", ERROR);
  }
  if (index_prefix.empty() && command != "train") {
    PrintMessageDieOnError("Must specify --index-prefix for classifying", ERROR);
  }
  // set sample if not specified
  if (sample.empty()) sample = output_prefix;
}

int main(int argc, char* argv[]) {
  time_t starttime, endtime;
  time(&starttime);
  PrintLobSTR();
  /* parse command line options */
  parse_commandline_options(argc, argv);
  PrintMessageDieOnError("Getting run info", PROGRESS);
  run_info.starttime = GetTime();
  if (_GIT_VERSION != NULL) {
    run_info.gitversion = _GIT_VERSION;
  } else {run_info.gitversion = "Not available";}
  if (_MACHTYPE != NULL) {
    run_info.machtype = _MACHTYPE;
  } else {run_info.machtype = "Not available";}
  run_info.params = user_defined_arguments_allelotyper;

  if (my_verbose) {
    PrintMessageDieOnError("Running allelotype with command " + command, PROGRESS);
  }
  /* Get haploid chromosomes */
  vector<string> haploid_chroms;
  boost::split(haploid_chroms,haploid_chroms_string,boost::is_any_of(","));

  /* initialize noise model */
  NoiseModel nm(strinfofile, haploid_chroms);

  /* Load ref character and ref object for each STR */
  if (my_verbose) {
    PrintMessageDieOnError("Loading reference STRs", PROGRESS);
  }
  vector<ReferenceSTR> reference_strs;
  if (command != "train") {
    FastaFileReader faReader(index_prefix+"ref.fa");
    MSReadRecord ref_record;
    while (faReader.GetNextRead(&ref_record)) {
      vector<string> items;
      string refstring = ref_record.ID;
      split(refstring, '$', items);
      if (items.size() == 7) {
	string chrom = items.at(1); 
        int start = atoi(items.at(2).c_str())+extend;
	int str_start = atoi(items.at(2).c_str());
	int str_end = atoi(items.at(3).c_str());
	ReferenceSTR ref_str;
	ref_str.start = str_start+extend;
	ref_str.stop = str_end-extend;
	ref_str.chrom = chrom;
	reference_strs.push_back(ref_str);
        string repseq_in_ref = items.at(6);
        string refnuc = ref_record.nucleotides.substr(extend, ref_record.nucleotides.length()-2*extend);
        pair<string, int> locus = pair<string,int>(chrom, start);
        ref_nucleotides.insert(pair< pair<string, int>, string>(locus, refnuc));
        ref_repseq.insert(pair< pair<string, int>, string>(locus, repseq_in_ref));
      }
    }
  }

  /* Load annotations */
  if (!annotation_files_string.empty()) {
    if (my_verbose) {
      PrintMessageDieOnError("Loading annotations", PROGRESS);
    }
    // TODO
  }

  vector<string>bam_files;
  boost::split(bam_files, bam_files_string, boost::is_any_of(","));

  /* Determine samples */
  if (my_verbose) {
    PrintMessageDieOnError("Determining samples to process", PROGRESS);
  }
  vector<string> samples_list;
  GetSamplesFromBamFiles(bam_files, &samples_list);
  if (samples_list.size() == 0) {
    PrintMessageDieOnError("Didn't find any read groups for samples in bam files", ERROR);
  }
  PrintMessageDieOnError("Debug",ERROR);

  /* Train/classify */
  if (command == "train") {
    ReadContainer read_container;
    ReferenceSTR dummy_ref_str;
    dummy_ref_str.chrom = "NA"; dummy_ref_str.start = -1; dummy_ref_str.stop = -1;
    read_container.AddReadsFromFile(bam_files, exclude_partial, dummy_ref_str);
    /* Perform pcr dup removal if specified */
    if (rmdup) {
      if (my_verbose) PrintMessageDieOnError("Performing PCR duplicate removal", PROGRESS);
      read_container.RemovePCRDuplicates();
    }
    if (my_verbose) PrintMessageDieOnError("Training noise model", PROGRESS);
    nm.Train(&read_container);
  } else if (command == "classify") {
    if (!nm.ReadNoiseModelFromFile(noise_model)) {
      PrintMessageDieOnError("Problem reading noise file", ERROR);
    }
  }
  if (command == "classify") {
    std::string tabfile = "";
    Genotyper genotyper(&nm, haploid_chroms, &ref_nucleotides, &ref_repseq,
			output_prefix + ".vcf");
    if (my_verbose) PrintMessageDieOnError("Classifying allelotypes", PROGRESS);
    for (size_t i = 0; i < reference_strs.size(); i++) {
      ReadContainer str_container;
      str_container.AddReadsFromFile(bam_files, exclude_partial, reference_strs.at(i));
      pair<string, int> coord(reference_strs.at(i).chrom, reference_strs.at(i).start);
      list<AlignedRead> aligned_reads;
      str_container.GetReadsAtCoord(coord, &aligned_reads);
      if (aligned_reads.size() > 0) {
	if (my_verbose) {
	  stringstream msg;
	  msg << "Processing locus " << reference_strs.at(i).chrom << ":"
	      << reference_strs.at(i).start << "-"
	      << reference_strs.at(i).stop;
	  PrintMessageDieOnError(msg.str(), PROGRESS);
	}
	if (rmdup) {
	  str_container.RemovePCRDuplicates();
	}
	genotyper.Genotype(aligned_reads);
      }
    }
  }
  run_info.endtime = GetTime();
  OutputRunStatistics();
  time(&endtime);
  stringstream msg;
  int seconds_elapsed = difftime(endtime, starttime);
  msg << "Done! " << seconds_elapsed/60/60/24 << ":"
      << (seconds_elapsed/60/60)%24 << ":" 
      << (seconds_elapsed/60)%60 << ":"
      << seconds_elapsed%60 << " elapsed";
  PrintMessageDieOnError(msg.str(), PROGRESS);
  return 0;
}
