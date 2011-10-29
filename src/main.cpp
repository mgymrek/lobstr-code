//============================================================================
// Author      : Melissa Gymrek
/*
TODO:
use cigar score
 */


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

#include "AlignmentFileReader.h"
#include "BWAReadAligner.h"
#include "bwtaln.h"
#include "bwase.h"
#include "common.h"
#include "FastaFileReader.h"
#include "FastqFileReader.h"
#include "FFT_nuc_vectors.h"
#include "FFT_four_nuc_vectors.h"
#include "Genotyper.h"
#include "HammingWindowGenerator.h"
#include "IFileReader.h"
#include "MSReadRecord.h"
#include "MultithreadData.h"
#include "SamFileWriter.h"
#include "STRDetector.h"
#include "runtime_parameters.h"
#include "TabFileWriter.h"
#include "TukeyWindowGenerator.h"

using namespace std;

// list of input files to process from
vector<string> input_files;

// Perform STRDetection
STRDetector str_detector;

// keep track of genotypes
// Genotyper genotyper;

// keep track of # bases so we can calculate coverage
long bases = 0;

// keep track of reference genome to use for sam format
map<string, int> chrom_sizes;

// alignment references, keep global
void LoadReference(const std::string& repseq);
void DestroyReferences();
map<std::string, BWT> bwt_references;
map<std::string, BNT> bnt_annotations;
gap_opt_t *opts;

static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;

void show_help(){
  const char* help = "lobSTR [OPTIONS] -f <file1[,file2,...]> --index-prefix <index prefix> -o <output prefix>\n" \
"-f,--files     file or comma-separated list of files containing reads in fasta or fastq format\n" \
"-o,--out       prefix for out put files. will output:\n" \
"--index-prefix     prefix for bwa reference (must run lobstr_index.py to create index)\n" \
"                      <prefix>.<filename>.fast(a/q): file of raw reads with detected reads modified for each <filename>\n" \
"                      <prefix>.msalign: file of alignments\n" \
"\n\nOptions:\n" \
"-h,--help                  display this help screen\n" \
"-v,--verbose               print out useful progress messages\n" \
"-q,--fastq                 reads are in fastq format (default: fasta)\n" \
"-p,--threads <int>         number of threads (default: 1)\n" \
"\n\nAdvanced options - general:\n" \
"--min-read-length    minimum number of nucleotides for a read to be processed. This should be at least two times fft-window-size (default: 48)\n" \
"--max-read-length    maximum number of nucleotides for a read to be processed. (default: 1000)\n" \
"\n\nAdvanced options - detection:\n" \
"--fft-window-size    size of fft window (default: 24)\n" \
"--fft-window-step    step size of sliding window (default: 12)\n" \
"--entropy-threshold     threshold score to call a window periodic (defualt: 0.3)\n" \
"--minperiod          minimum period to attempt to detect (default: 1)\n" \
"--maxperiod          maximum period to attempt to detect (default: 8)\n" \
"--minflank           minimum length of flanking region to try to align (default: 10)\n" \
"--maxflank           length to trim flanking regions to if they exceed that length (default: 1000)\n" \
"\n\nAdvanced options - alignment:\n" \
"--max-diff-ref       maximum difference in length from the reference sequence to allow for alignment (default 50) (will take the absolute value)\n" \
"-u                   require length difference to be a multiple of the repeat unit\n" \
    "-m,--mismatch <int>  edit distance allowed during alignment of each flanking region (default: 0)\n" \
    "-g <int>             maximum number of gap opens allowed in each flanking region (default: 0)\n" \
    "-e <int>             maximum number of gap extensions allowed in each flanking region (default: 0)\n" \
    "-r <float>           edit distance allowed during alignment of each flanking region (ignored if -m is set) (default 0)\n" \
"This program takes in raw reads, detects and aligns reads containing microsatellites, and genotypes STR locations.\n\n";
	cout << help;
	exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc,char* argv[]) {
	enum LONG_OPTIONS{
	  OPT_FILES,
	  OPT_TABLE,
	  OPT_GENOME,
	  OPT_OUTPUT,
	  OPT_HELP,
	  OPT_VERBOSE,
	  OPT_DEBUG,
	  OPT_ALIGN_DEBUG,
	  OPT_FASTQ,
	  OPT_THREADS,
	  OPT_MISMATCH,
	  OPT_SAM,
	  OPT_RMDUP,
	  OPT_FFT_WINDOW_SIZE,
	  OPT_FFT_WINDOW_STEP,
	  OPT_LOBE_THRESHOLD,
	  OPT_EXTEND,
	  OPT_MAX_PERIOD,
	  OPT_MIN_PERIOD,
	  OPT_MIN_FLANK_LEN,
	  OPT_MAX_FLANK_LEN,
	  OPT_MAX_DIFF_REF,
	  OPT_FFTW_DEBUG,
	  OPT_LOBE_DEBUG,
	  OPT_GENOTYPER_DEBUG,
	  OPT_PI,
	  OPT_MU,
	  OPT_MIN_READ_LENGTH,
	  OPT_MAX_READ_LENGTH,
	  OPT_ERROR_RATE,
	  OPT_FEMALE,
	  OPT_MIN_COVERAGE,
	  OPT_GENOTYPE_ONLY,
	  OPT_ALIGNED_FILE,
	  OPT_WHY_NOT,
	  OPT_USE_ENTROPY,
	  OPT_ENTROPY_THRESHOLD,
	  OPT_DEBUG_ENTROPY,
	  OPT_PROFILE,
	  OPT_INDEX,
	  OPT_GAP_OPEN,
	  OPT_GAP_EXTEND,
	  OPT_FPR,
	};

	int ch;
	int option_index = 0;

	static struct option long_options[] = {
	  {"files", 1, 0, OPT_FILES},
	  {"table", 1, 0, OPT_TABLE},
	  {"genome", 1, 0, OPT_GENOME},
	  {"out", 1, 0, OPT_OUTPUT},
	  {"threads", 1, 0, OPT_THREADS},
	  {"mismatch", 1, 0, OPT_MISMATCH},
	  {"bam", 0, 0, OPT_SAM},
	  {"fft-window-size", 1, 0, OPT_FFT_WINDOW_SIZE},
	  {"fft-window-step", 1, 0, OPT_FFT_WINDOW_STEP},
	  {"lobe-threshold", 1, 0, OPT_LOBE_THRESHOLD},
	  {"extend", 1, 0, OPT_EXTEND},
	  {"maxperiod", 1, 0, OPT_MAX_PERIOD},
	  {"minperiod", 1, 0, OPT_MIN_PERIOD},
	  {"minflank", 1, 0, OPT_MIN_FLANK_LEN},
	  {"maxflank", 1, 0, OPT_MAX_FLANK_LEN},
	  {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
	  {"pi", 1, 0, OPT_PI},
	  {"mu", 1, 0, OPT_MU},
	  {"help", 0, 0, OPT_HELP},
	  {"verbose", 0, 0, OPT_VERBOSE},
	  {"debug", 0, 0, OPT_DEBUG},
	  {"fastq", 0, 0, OPT_FASTQ},
	  {"rmdup", 0, 0, OPT_RMDUP},
	  {"fftw-debug", 0, 0, OPT_FFTW_DEBUG},
	  {"lobe-debug", 0, 0, OPT_LOBE_DEBUG},
	  {"align-debug", 0, 0, OPT_ALIGN_DEBUG},
	  {"genotyper-debug", 0, 0, OPT_GENOTYPER_DEBUG},
	  {"min-read-length", 1, 0, OPT_MIN_READ_LENGTH},
	  {"max-read-length", 1, 0, OPT_MAX_READ_LENGTH},
	  {"min-coverage", 1, 0, OPT_MIN_COVERAGE},
	  {"female", 0, 0, OPT_FEMALE},
	  {"genotype-only", 0, 0, OPT_GENOTYPE_ONLY},
	  {"aligned-file", 1, 0, OPT_ALIGNED_FILE},
	  {"why-not", 0, 0, OPT_WHY_NOT},
	  {"use-entropy", 0, 0, OPT_USE_ENTROPY},
	  {"entropy-threshold", 1, 0, OPT_ENTROPY_THRESHOLD},
	  {"entropy-debug", 0, 0, OPT_DEBUG_ENTROPY},
	  {"profile", 0, 0, OPT_PROFILE},
	  {"index-prefix", 1, 0, OPT_INDEX},
	};
	while ((ch = getopt_long(argc,argv,"hvqp:f:t:g:o:m:s:d:e:g:r:u",
				 long_options,&option_index))!= -1) { 
	  switch(ch){
	  case 'u':
	    unit++;
	  case OPT_PROFILE:
	    profile++;
	    break;
	  case OPT_ALIGN_DEBUG:
	    align_debug ++;
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
	  fastq++;
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
	  case 'm':
	  case OPT_MISMATCH:
	    allowed_mismatches = atoi(optarg);
	    if (allowed_mismatches < 0)
	      errx(1,"Error: invalid number of mismatches");
	    break;
	  case 'b':
	  case OPT_SAM:
	    sam++;
	    break;
	  case OPT_FFT_WINDOW_SIZE:
	    fft_window_size = atoi(optarg);
	    break;
	  case OPT_FFT_WINDOW_STEP:
	    fft_window_step = atoi(optarg);
	    break;
	  case OPT_LOBE_THRESHOLD:
	    fft_lobe_threshold = atof(optarg);
	    if (fft_lobe_threshold <=0)
	      errx(1,"Error: invalid lobe threshold");
	    break;
	  case OPT_EXTEND:
	    extend = atoi(optarg);
	    if (extend <=0)
	      errx(1, "Error: invalid extension length");
	    break;
	  case OPT_MIN_PERIOD:
	    min_period = atoi(optarg);
	    if (min_period <= 0)
	      errx(1, "Error: invalid min period");
	    break;
	  case OPT_MAX_PERIOD:
	    max_period = atoi(optarg);
	    if (max_period <= 0)
	      errx(1, "Error: invalid max period");
	    break;
	  case OPT_MAX_FLANK_LEN:
	    max_flank_len = atoi(optarg);
	    if (max_flank_len <= 0)
	      errx(1, "Error: invalid max flank length");
	    break;
	  case OPT_MIN_FLANK_LEN:
	    min_flank_len = atoi(optarg);
	    if (min_flank_len <= 0)
	      errx(1, "Error: invalid min flank length");
	    break;
	  case OPT_MAX_DIFF_REF:
	    max_diff_ref = atoi(optarg);
	    if (max_diff_ref <=0 )
	      errx(1, "Error: invalid max diff ref");
	    break;
	  case OPT_RMDUP:
	    rmdup++;
	    break;
	  case OPT_FFTW_DEBUG:
	    fftw_debug = true;
	    break;
	  case OPT_LOBE_DEBUG:
	    lobe_debug = true;
	    break;
	  case OPT_GENOTYPER_DEBUG:
	    genotyper_debug = true;
	    break;
	  case OPT_PI:
	    pi_string = string(optarg);
	    break;
	  case OPT_MU:
	    mu_string = string(optarg);
	    break;
	  case OPT_MIN_READ_LENGTH:
	    min_read_length = atoi(optarg);
	    break;
	  case OPT_MAX_READ_LENGTH:
	    max_read_length = atoi(optarg);
	    break;
	  case OPT_MIN_COVERAGE:
	    min_coverage = atoi(optarg);
	    break;
	  case OPT_FEMALE:
	    male = false;
	    break;
	  case OPT_GENOTYPE_ONLY:
	    genotype_only++;
	    break;
	  case OPT_ALIGNED_FILE:
	    aligned_file = string(optarg);
	    break;
	  case OPT_WHY_NOT:
	    why_not_debug++;
	    break;
	  case OPT_USE_ENTROPY:
	    use_entropy++;
	    break;
	  case OPT_ENTROPY_THRESHOLD:
	    entropy_threshold = atof(optarg);
	    break;
	  case OPT_DEBUG_ENTROPY:
	    entropy_debug++;
	    break;
	  case OPT_INDEX:
	    index_prefix = string(optarg);
	    break;
	  case OPT_GAP_OPEN:
	  case 'g':
	    gap_open = atoi(optarg);
	    break;
	  case OPT_GAP_EXTEND:
	  case 'e':
	    gap_extend = atoi(optarg);
	    break;
	  case OPT_FPR:
	  case 'r':
	    fpr = atof(optarg);
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
	// make sure arguments make sense
	if (fft_window_step > fft_window_size) {
	  errx(1, "fft_window_step must be <=fft_window_size");
	}
	if (min_period > max_period) {
	  errx(1, "min_period must be <= max_period");
	}
	if (min_flank_len > max_flank_len) {
	  errx(1, "min_flank_len must be <=max_flank_len");
	}
	// check that we have the mandatory parameters
	if ((input_files_string.empty()||  output_prefix.empty() ||
	     index_prefix.empty()) &&
	    !genotype_only) {
	  errx(1, "Required arguments are mising");
	} else if (genotype_only && (aligned_file.empty()
				     || output_prefix.empty())) {
	  errx(1, "Required arguments are missing");
	}
}

void LoadReference(const std::string& repseq) {
  // Load BWT index
  char* prefix = (char*)calloc(strlen(index_prefix.c_str()) +
			       strlen(repseq.c_str()) + 20, 1);
  strcpy(prefix, index_prefix.c_str()); strcat(prefix, repseq.c_str());
  strcat(prefix, ".fa");
  cout << "loading ref " << repseq << endl;
  bwt_t *bwt[2];
  char *str = (char*)calloc(strlen(prefix) + 10, 1);
  strcpy(str, prefix); strcat(str, ".bwt");
  bwt[0] = bwt_restore_bwt(str);
  strcpy(str, prefix); strcat(str, ".rbwt");
  bwt[1] = bwt_restore_bwt(str);
  strcpy(str, prefix); strcat(str, ".sa");
  bwt_restore_sa(str, bwt[0]);
  strcpy(str, prefix); strcat(str, ".rsa");
  bwt_restore_sa(str, bwt[1]);

  BWT bwt_ref;
  bwt_ref.bwt[0] = bwt[0];
  bwt_ref.bwt[1] = bwt[1];
  bwt_references.insert(pair<string, BWT>(repseq,bwt_ref));

  // Load BNT annotations
  BNT bnt;
  bntseq_t *bns, *ntbns = 0;
  bns = bns_restore(prefix);
  bnt.bns = bns;
  bnt_annotations.insert(pair<string, BNT>(repseq,bnt));
  free(str);
  free(prefix);
}

void DestroyReferences() {
  for (map<string, BWT>::iterator it = bwt_references.begin();
       it != bwt_references.end(); ++it) {
    bwt_destroy(it->second.bwt[0]);
    bwt_destroy(it->second.bwt[1]);
  }
  for (map<string, BNT>::iterator it = bnt_annotations.begin();
       it != bnt_annotations.end(); ++it) {
    bns_destroy(it->second.bns);
  }
}


/*
 * process read in single thread
 */
void single_thread_process_loop(const vector<string>& files) {
  MSReadRecord msread;
  TabFileWriter pWriter(output_prefix + ".aligned.tab");
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references, &bnt_annotations, opts);
  for (vector<string>::const_iterator it = files.begin(); it != files.end(); it++) {
    if (fexists((*it).c_str())) {
      cout << "processing file " <<  *it << " ...\n";
      IFileReader* pReader = create_file_reader(*it);
      int aligned = false;
      while (pReader->GetNextRecord(&msread)) {
	aligned = false;
	msread.msRepeat = "";
	if (!(msread.nucleotides.length() >= min_read_length &&
	      msread.nucleotides.length() <= max_read_length)) {
	  continue;
	} else {
	  bases += msread.nucleotides.length();
	}
	if (!pDetector->ProcessRead(&msread)) {
	  //	  	   cout << "not aligned (not detected) " << msread.nucleotides << endl;
	  continue;
	}
	string repseqfw;
	string repseqrev;
	string repseq;
	getMSSeq(msread.detected_ms_region_nuc, msread.ms_repeat_best_period, &repseqfw);
	getCanonicalMS(reverseComplement(repseqfw), &repseqrev);
	repseq = getFirstString(repseqfw, repseqrev);
	msread.repseq = repseq;
	msread.repseq_reverse = (repseq==repseqrev);
	// alignment
	if (pAligner->ProcessRead(&msread)) {
	  aligned = true;
	} else if (msread.ms_repeat_next_best_period != 0) {
	  msread.ms_repeat_best_period = msread.ms_repeat_next_best_period;
	  getMSSeq(msread.detected_ms_region_nuc, msread.ms_repeat_best_period, &repseqfw);
	  getCanonicalMS(reverseComplement(repseqfw), &repseqrev);
	  repseq = getFirstString(repseqfw, repseqrev);
	  msread.repseq = repseq;
	  msread.repseq_reverse = (repseq==repseqrev);
	  if (pAligner->ProcessRead(&msread)) {
	    aligned = true;
	  }
	} 
	if (aligned) {
	  // genotyper.AddRead(&msread);
	  pWriter.WriteRecord(msread);
	  if (sam) samWriter.WriteRecord(msread);
	}
      }
      delete pReader;
    } else {
      cerr << "Warning: file " << *it << " does not exist" << endl;
    }
  }
  delete pDetector;
  delete pAligner;
}

void* satellite_process_consumer_thread(void *arg) {
  MultithreadData *pMT_DATA = (MultithreadData*)(arg);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references, &bnt_annotations,opts);
  int aligned = false;
  while (1) {
    aligned = false;
    MSReadRecord *pReadRecord = pMT_DATA->get_new_input();
    if (pReadRecord->nucleotides.length() >= min_read_length &&
	pReadRecord->nucleotides.length() <= max_read_length) {
      bases += pReadRecord->nucleotides.length();
      if (pDetector->ProcessRead(pReadRecord)) {
	string repseqfw;
	string repseqrev;
	string repseq;
	getMSSeq(pReadRecord->detected_ms_region_nuc, pReadRecord->ms_repeat_best_period, &repseqfw);
	getCanonicalMS(reverseComplement(repseqfw), &repseqrev);
	repseq = getFirstString(repseqfw, repseqrev);
	pReadRecord->repseq = repseq;
	pReadRecord->repseq_reverse = (repseq == repseqrev);
	// alignment
	if (pAligner->ProcessRead(pReadRecord)) {
	  aligned = true;
	} else if (pReadRecord->ms_repeat_next_best_period != 0) {
	  pReadRecord->ms_repeat_best_period = pReadRecord->ms_repeat_next_best_period;
	  getMSSeq(pReadRecord->detected_ms_region_nuc, pReadRecord->ms_repeat_best_period, &repseqfw);
	  getCanonicalMS(reverseComplement(repseqfw), &repseqrev);
	  repseq = getFirstString(repseqfw, repseqrev);
	  pReadRecord->repseq = repseq;	pReadRecord->repseq_reverse = (repseq == repseqrev);
	  if (pAligner->ProcessRead(pReadRecord)) {
	    aligned = true;
	  }
	}
      }
    }	  
    if (aligned) {
      //genotyper.AddRead(pReadRecord);
      pMT_DATA->post_new_output_read(pReadRecord);
    } else {
      //Don't pass this read on to the output thread,
      //discard it here
      delete pReadRecord;
      pMT_DATA->increment_output_counter();
    }
  }
  pthread_exit((void*) arg);
}

void* output_writer_thread(void *arg) {
  MultithreadData *pMT_DATA = (MultithreadData*)(arg);
  
  TabFileWriter *pWriter = new TabFileWriter(output_prefix + ".aligned.tab");
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  while (1) {
    MSReadRecord *pReadRecord = pMT_DATA->get_new_output();
    pWriter->WriteRecord(*pReadRecord);
    samWriter.WriteRecord(*pReadRecord);
    delete pReadRecord;
    pMT_DATA->increment_output_counter();
  }
  delete pWriter;
  delete pMT_DATA;
  pthread_exit((void*) arg);
}

void multi_thread_process_loop(vector<string> files) {
  MultithreadData mtdata(threads);
  list<pthread_t> satellite_threads;
  pthread_t writer_thread;
  if (files.size() == 0) return;
  
  for (size_t i=0; i<threads; ++i) {
    pthread_t id;
    if (pthread_create(&id, NULL, satellite_process_consumer_thread, (void*)&mtdata))
      err(1,"Failed to create thread");
    satellite_threads.push_back(id);
  }
  
  if (pthread_create(&writer_thread, NULL, output_writer_thread, (void*)&mtdata))
    err(1,"failed to create output writer thread");
  
  size_t counter = 1 ;
  
  for (vector<string>::const_iterator it = files.begin();
       it != files.end(); ++it) {
    if (fexists((*it).c_str())) {
      cout << "processing file " << *it << endl;
      IFileReader *pReader = create_file_reader(*it); 
      
      do {
	MSReadRecord *pRecord = new MSReadRecord;
	if (!pReader->GetNextRecord(pRecord))
	  break; //no more reads
	
	pRecord->read_counter = counter;
	counter++;
	mtdata.increment_input_counter();
	
	mtdata.post_new_input_read(pRecord);
	pRecord = NULL ; //the consumers will take it from here, and free it
      } while(1);
      delete pReader;
    } else {
      cerr << "File " << *it << " does not exist" << endl;
    }
  }
  while ( 1 ) {
    sleep(1); //OMG, the horror...
    if ( mtdata.input_output_counters_equal())
      break;
    sleep(10);
    break;
  }
}


int main(int argc,char* argv[]) {
  parse_commandline_options(argc,argv);
  
  if (verbose) cout << "Initializing..." << endl;
  // open file with all names
  TextFileReader tReader(index_prefix+"strdict.txt");
  string line;
  while(tReader.GetNextLine(&line)) {
    // set chrom sizes
    if (line.substr(0,3) == "REF") {
      vector<string> items;
      split(line,'\t', items);
      chrom_sizes[items[1]] = atoi(items[2].c_str());
    } else {
      LoadReference(line);
    }
  }
  // set up options
  opts = gap_init_opt();
  // edit distance per flank
  //opts->fnr = 0.1;
  opts->max_diff = allowed_mismatches;
  if (allowed_mismatches == 0) {
    opts->fnr = fpr;
  } else {
    opts->fnr = -1.0;
  }
  // only 1 indel per flank
  opts->max_gapo=gap_open;
  opts->max_gape=gap_extend;
  // take the first INT subsequence as seed
  opts->seed_len = 5;
  // all hits with no more than maxdiff found
  // maybe remove, supposed to be slow
  //opts->mode |= BWA_MODE_NONSTOP; opts->max_top2 = 0x7fffffff;



  if (!genotype_only) {
    // get the input files
    boost::split(input_files, input_files_string, boost::is_any_of(","));
    // Initialize fft  
    // Create the singleton HammingWindow (before starting any threads)
    HammingWindowGenerator* hamgen =
      HammingWindowGenerator::GetHammingWindowSingleton();
    TukeyWindowGenerator* tukgen =
      TukeyWindowGenerator::GetTukeyWindowSingleton();
    
    // Initialize global FFTW plans
    FFT_NUC_VECTOR::initialize_fftw_plans();
    
    // run detection/alignment
   if (verbose) {cerr << "Running detection/alignment..." << endl;}
    if (threads == 1) {
      single_thread_process_loop(input_files);
    } else {
      multi_thread_process_loop(input_files);
    }
    delete hamgen;
    delete tukgen;
  } /*else {
    // read alignment file and add reads to genotyper
    AlignmentFileReader alignment_reader(aligned_file);
    MSReadRecord aligned_read_record;
    // get rid of header
    while(alignment_reader.GetNextRecord(&aligned_read_record)) {
      genotyper.AddRead(&aligned_read_record);
    }
  }
  // init priors
  if (!mu_string.empty()) {
    vector<string> mu_values;
    boost::split(mu_values, mu_string, boost::is_any_of(","));
    mu_0 = atof(mu_values.at(0).c_str());
    mu_1 = atof(mu_values.at(1).c_str());
    mu_2 = atof(mu_values.at(2).c_str());
    genotyper.ResetMu(mu_0, mu_1, mu_2);
  }
  if (!pi_string.empty()) {
    vector<string> pi_values;
    boost::split(pi_values, pi_string, boost::is_any_of(","));
    pi_0 = atof(pi_values.at(0).c_str());
    pi_1 = atof(pi_values.at(1).c_str());
    pi_2 = atof(pi_values.at(2).c_str());
    genotyper.ResetPi(pi_1, pi_1, pi_2);
    }

  // run genotyping and write output
  genotyper.WriteOutput(output_prefix + ".genotypes.tab");*/

  return 0;
}
