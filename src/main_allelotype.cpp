//============================================================================
// Author      : Melissa Gymrek
//============================================================================
#include <err.h>
#include <error.h>
#include <getopt.h>
#include <iostream>
#include <RInside.h>
#include <stdlib.h>
#include <string>

#include "Genotyper.h"
#include "NoiseModel.h"
#include "runtime_parameters.h"
#include "ReadContainer.h"

using namespace std;

void show_help(){
  const char* help = "\n\nTo train the genotyping noise model from a set of aligned reads:\n" \
"allelotype --command train [OPTIONS] --bam <input.bam> --noise_model <noisemodel.txt>\n\n" \
"To run str profiling on a set of aligned reads:\n" \
"allelotype --command classify [OPTIONS] --bam <input.bam> --noise_model <noisemodel.txt> [--no-rmdup] --out <output_prefix> --sex [M|F]\n\n" \
"To run training and classification on a set of aligned reads:\n" \
"allelotype --command both [OPTIONS] --bam <input.bam> --noise_model <noisemodel.txt> [--no-rmdup] --out <output_prefix> --sex [M|F]\n\n" \
"To allelotype without using a stutter noise model:\n" \
"allelotype simple [OPTIONS] --bam <input.bam> [--no-rmdup] --out <output_prefix> --sex [M|F]\n\n" \
"Options:\n" \
"--no-rmdup: don't remove pcr duplicates before allelotyping\n" \
"-h: display this message\n" \
    "-v: print out helpful progress messages\n\n";
  cout << help;
  exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc,char* argv[]) {
	enum LONG_OPTIONS{
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
	};

	int ch;
	int option_index = 0;

	static struct option long_options[] = {
	  {"bam", 1, 0, OPT_BAM},
	  {"command", 1, 0, OPT_COMMAND},
	  {"out", 1, 0, OPT_OUTPUT},
	  {"no-rmdup", 0, 0, OPT_NORMDUP},
	  {"sex", 1, 0, OPT_SEX},
	  {"noise_model", 1, 0, OPT_NOISEMODEL},
	  {"unit", 0, 0, OPT_UNIT},
	  {"help", 1, 0, OPT_HELP},
	  {"debug", 0, 0, OPT_DEBUG},
	  {NULL, no_argument, NULL, 0},
	};

	ch = getopt_long(argc,argv,"hv?",
			 long_options,&option_index);
	while (ch != -1) { 
	  switch(ch){
	  case OPT_BAM:
	    bam_file = string(optarg);
	    break;
	  case OPT_COMMAND:
	    command = string(optarg);
	    if (command != "train" & command != "classify"
		& command != "both" & command != "simple") {
	      cerr << "\n\nERROR: Command " << command << " is invalid. Command must be one of: train, classify, both, simple";
	      show_help();
	      exit(1);
	    }
	    break;
	  case OPT_OUTPUT:
	    output_prefix = string(optarg);
	    break;
	  case OPT_NORMDUP:
	    rmdup = false;
	    break;
	  case OPT_SEX:
	    if (string(optarg) != "F" && string(optarg) != "M") {
	      errx(1, "--sex must be F or M");
	    }
	    if (string(optarg) == "F") male = false;
	    sex_set++;
	    break;
	  case OPT_NOISEMODEL:
	    noise_model = string(optarg);
	    break;
	  case OPT_UNIT:
	    unit = true;
	    break;
	  case 'v':
	  case OPT_VERBOSE:
	    verbose++;
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
	  ch = getopt_long(argc,argv,"hv?",
			   long_options,&option_index);

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
	  if (bam_file.empty() || noise_model.empty()) {
	    cerr << "\n\nERROR: Required arguments are missing. Please specify a bam file and an output prefix";
	    show_help();
	    exit(1);
	  }
	  male = true;
	}
	if (command == "classify") {
	  if (bam_file.empty() || noise_model.empty()
	      or output_prefix.empty() || !sex_set) {
	    cerr << "\n\nERROR: Required arguments are missing. Please specify a bam file, output prefix, noise model, and gender";
	    show_help();
	    exit(1);
	  }
	}
	if (command == "both" || command == "simple") {
	  if (bam_file.empty() || output_prefix.empty() || !sex_set)  {
	    cerr << "\n\nERROR: Required arguments are missing. Please specify a bam file, output prefix, and gender";
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


int main(int argc,char* argv[]) {
  /* parse command line options */
  parse_commandline_options(argc,argv);

  if (verbose) cout << "Running allelotyping with command "
		    << command << endl;

  /* initialize noise model */
  RInside R(argc, argv);
  NoiseModel nm(&R);

  /* Add reads to read container */
  if (verbose) cout << "Adding reads to read container" << endl;
  ReadContainer read_container;
  read_container.AddReadsFromFile(bam_file);

  /* Perform pcr dup removal if specified */
  if (rmdup) {
    if (verbose) cout << "Performing pcr duplicate removal" << endl;
    read_container.RemovePCRDuplicates();
  }

  /* Train/classify */
  if (command == "train" or command == "both") {
    if (verbose) cout << "Training noise model..." << endl;
    nm.Train(&read_container);
    nm.WriteNoiseModelToFile(noise_model);
  } else if (command == "classify") {
    if (!nm.ReadNoiseModelFromFile(noise_model)) {
      errx(1,"Error reading noise file");
    }
  }

  if (command == "classify" or command == "both") {
    if (verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, male, false);
    genotyper.Genotype(read_container,
		       output_prefix + ".genotypes.tab");
  }
  if (command == "simple") {
    if (verbose) cout << "Classifying allelotypes..." << endl;
    Genotyper genotyper(&nm, male, true);
    genotyper.Genotype(read_container,
		       output_prefix + ".genotypes.tab");
  }
  return 0;
}
