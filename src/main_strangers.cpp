//============================================================================
// Author      : Melissa Gymrek
//============================================================================

#include <boost/algorithm/string.hpp>
#include <err.h>
#include <error.h>
#include <getopt.h>
#include <iostream>
#include <limits.h>
#include <list>
#include <map>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "Anonymizer.h"
#include "AnonymizerMultiThreadData.h"
#include "runtime_parameters.h"

using namespace std;

// list of input files to process from
vector<string> input_files;
Anonymizer anonymizer;

void show_help(){
  const char* help = "strangers [OPTIONS] -f <file1[,file2,...]> -a <lobSTR alignment file> [-d <allele database filename> | --mask ]\n" \
"-o <output prefix>\n" \
"-f,--files     file or comma-separated list of files containing reads in fasta or fastq format\n" \
"-a, --str_alignment  file output from lobSTR containing STR alignments \n" \
"-d, --database database file of alleles to choose from during anonymization \n" \
"--mask,            mask reads with N's instead of choosing new allele. \n"\
"-o,--out       prefix for out put files. will output:\n" \
"                      <prefix>.<filename>.fast(a/q): file of raw reads with detected reads modified for each <filename>\n" \
"                      <prefix>.<filename>.ids: file with list of identifiers of reads that were changed\n"
"                      <prefix>.haplotype: allele database containing alleles of the randomized haplotype\n"
"-h,--help                  display this help screen\n" \
"-v,--verbose               print out useful progress messages\n" \
"-q,--fastq                 reads are in fastq format (default: fasta)\n" \
"-p,--threads <int>         number of threads (default: 1)\n\n" \
"--error-rate         error rate to simulate errors in anonymized reads (default 0.0001 substitutions/bp\n" \
"This program takes in raw reads and the alignment produced by lobSTR and randomizes the STR genotypes in detected reads.";
	cout << help;
	exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc,char* argv[]) {
	enum LONG_OPTIONS{
	  OPT_ANONYMIZER_DEBUG,
	  OPT_FILES,
	  OPT_ALIGNMENT,
	  OPT_DATABASE,
	  OPT_OUTPUT,
	  OPT_THREADS,
	  OPT_HELP,
	  OPT_VERBOSE,
	  OPT_DEBUG,
	  OPT_FASTQ,
	  OPT_MASKN,
	  OPT_ERROR_RATE,
	};

	int ch;
	int option_index = 0;

	static struct option long_options[] = {
	  {"anonymizer-debug", 0, 0, OPT_ANONYMIZER_DEBUG},
	  {"files", 1, 0, OPT_FILES},
	  {"alignment", 1, 0, OPT_ALIGNMENT},
	  {"database", 1, 0, OPT_DATABASE},
	  {"out", 1, 0, OPT_OUTPUT},
	  {"threads", 1, 0, OPT_THREADS},
	  {"help", 0, 0, OPT_HELP},
	  {"verbose", 0, 0, OPT_VERBOSE},
	  {"debug", 0, 0, OPT_DEBUG},
	  {"fastq", 0, 0, OPT_FASTQ},
	  {"mask", 0, 0, OPT_MASKN},
	  {"error-rate", 1, 0, OPT_ERROR_RATE},
	};
	while ((ch = getopt_long(argc,argv,"hvqp:f:o:s:d:a:",
				 long_options,&option_index))!= -1) { 
	  switch(ch){
	  case OPT_ERROR_RATE:
	    error_rate = atof(optarg);
	    if (error_rate < 0 || error_rate > 1) {
	      errx(1,"Error rate must be between 0 and 1");
	    }
	    break;
	  case OPT_ANONYMIZER_DEBUG:
	    anonymizer_debug++;
	    break;
	  case 'a':
	  case OPT_ALIGNMENT:
	    str_alignment_file = optarg;
	    break;
	  case 'd':
	  case OPT_DATABASE:
	    database_file = optarg;
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
	  case 'q':
	  case OPT_FASTQ:
	    input_type = INPUT_FASTQ;
	    break;
	  case OPT_MASKN:
	    mask_n++;
	    break;
	  case 'p':
	  case OPT_THREADS:
	    threads = atoi(optarg);
	    if (threads <=0)
	      errx(1,"Error: invalid number of threads");
	    break;
	  case 'f':
	  case OPT_FILES:
	    input_files_string = optarg;
	    break;
	  case 'o':
	  case OPT_OUTPUT:
	    output_prefix = string(optarg);
	    break;
	  default:
	    show_help();
	    exit(1);
	  }
	}
	// any arguments left over are extra
	if (optind < argc) {
	  cout << "Unnecessary leftover arguments...\n";
	  show_help();
	  exit(1);
	}
	// check that we have the mandatory parameters
	if (input_files_string.empty() || output_prefix.empty()) {
	  errx(1, "Required arguments are mising");
	}
}

/*
 * process read in single thread
 */
void single_thread_process_loop(const vector<string>& files) {
  for (vector<string>::iterator it = input_files.begin();
       it != input_files.end() ; it++) {
    // call anonymizer
    vector<string> items;
    boost::split(items, *it, boost::is_any_of("/"));
    string filename = items.back();
    cout << "processing file " << *it << " ..." << endl;
    anonymizer.AnonymizeReads(*it, output_prefix + "."+ filename + ".anonymized",
			      output_prefix + "."+ filename + ".modified_read_ids");
  }  
}

void* satellite_process_consumer_thread(void *arg)
{
  AnonymizerMultithreadData *pMT_DATA = (AnonymizerMultithreadData*)(arg);
  while (1) {
    string whole_filename = pMT_DATA->get_new_input();
    // call anonymizer on file
    vector<string> items;
    boost::split(items, whole_filename, boost::is_any_of("/"));
    string filename = items.back();

    anonymizer.AnonymizeReads(whole_filename, output_prefix + filename + ".anonymized",
			      output_prefix + filename + ".modified_read_ids");
    pMT_DATA->increment_output_counter();
  }
}

void multi_thread_process_loop(vector<string> files)
{
  AnonymizerMultithreadData mtdata(threads);
  list<pthread_t> satellite_threads;
  pthread_t writer_thread;

  for (size_t i=0; i<threads; ++i) {
    pthread_t id;
    if (pthread_create(&id, NULL, satellite_process_consumer_thread, (void*)&mtdata))
      err(1,"Failed to create thread");
    satellite_threads.push_back(id);
  }
  
  size_t counter = 1 ;

  for (vector<string>::iterator it = input_files.begin();
       it != input_files.end() ; it++) {
    counter++;
    mtdata.increment_input_counter();
    mtdata.post_new_input(*it);
  }

  while ( 1 ) {
    sleep(1); //OMG, the horror...
    if ( mtdata.input_output_counters_equal())
      break;
  }
}

int main(int argc,char* argv[]) {
  parse_commandline_options(argc,argv);
  // get the input files
  boost::split(input_files, input_files_string, boost::is_any_of(","));

  // get the map of identifier -> read record
  anonymizer.CreateReadIDToRecordMap(str_alignment_file);

  // create the YDatabase
  if (!mask_n)
    anonymizer.CreateDatabase(database_file);

  if (threads == 1) {
    single_thread_process_loop(input_files);
  } else {
    multi_thread_process_loop(input_files);
  }
  return 0;
}
