/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 

#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "src/BamFileReader.h"
#include "src/BamPairedFileReader.h"
#include "src/common.h"
#include "src/FastaFileReader.h"
#include "src/FastaPairedFileReader.h"
#include "src/FastqFileReader.h"
#include "src/FastqPairedFileReader.h"
#include "src/TextFileWriter.h"
#include "src/ZippedFastaFileReader.h"
#include "src/ZippedFastqFileReader.h"
#include "src/runtime_parameters.h"

using namespace std;

const char* NUCLEOTIDES[4] = {"A","C","G","T"};

void PrintLobSTR() {
  stringstream msg;
  msg << endl << endl;
  msg << "                       _______" << endl;
  msg << "             \\\\ //  /     -^--\\ |" << endl;
  msg << "             ||||  / /\\_____/ /" << endl;
  msg << " {\\         ______{ }        /         lobSTR: profiling short tandem repeats" << endl;
  msg << " {_}{\\{\\{\\{|         \\=@____/          from high-throughput sequencing data" << endl;
  msg << "<{_{-{-{-{-| ====---- >>>" << endl;
  msg << " { }{/{/{/{|______  _/=@_____" << endl;
  msg << " {/               { }        \\         Copyright (C) 2011-2014 Melissa Gymrek" << endl;
  msg << "            ||||  \\ \\______  \\         <mgymrek@mit.edu>" << endl;
  msg << "             // \\\\  \\    _^_\\  |" << endl;
  msg << "                     \\______/" << endl << endl;
  cout << msg.str();
}

void AddOption(const string& optname, const string& optval,
               bool hasvalue, string* paramstring) {
  *paramstring += optname;
  if (hasvalue) {
    *paramstring += "=";
    *paramstring += optval;
  }
  *paramstring += ";";
  return;
}

void OutputRunStatistics() {
  PrintMessageDieOnError("Outputting run statistics", PROGRESS);
  TextFileWriter sWriter(output_prefix + (program == LOBSTR? ".aligned.stats":".allelotype.stats"));
  // Output run statistics to stats file
  string stats_string = run_info.PrintToString((program == LOBSTR ? 0: 1), filter_counter);
  sWriter.Write(stats_string);
}

void PrintMessageDieOnError(const string& msg, MSGTYPE msgtype) {
  string typestring = "";
  switch (msgtype) {
  case ERROR:
    typestring = "ERROR: ";
    break;
  case WARNING:
    typestring = "WARNING: ";
    break;
  case PROGRESS:
    typestring = "ProgressMeter: ";
    break;
  case DEBUG:
    typestring = "DEBUG: ";
    break;
  default:
    errx(1,"Invalid message type. This should never happen");
  }
  stringstream ss;
  ss  << "[" << (program == LOBSTR ? "lobSTR":"allelotype")
      << "-" << _GIT_VERSION << "] " << currentDateTime() << " " << typestring << msg << endl;
  if (!quiet) {
    cerr << ss.str();
  }
  if (msgtype == ERROR) {
    run_info.error = ss.str();
    run_info.endtime = GetTime();
    if (!output_prefix.empty()) {
      OutputRunStatistics();
    }
    exit(1);
  }
}

void CheckIndexVersion() {
  if (fexists((index_prefix+"strdict.txt").c_str())) {
    PrintMessageDieOnError("It appears you are using an outdated lobSTR index. Please " \
			   "use the index for version 3.0.0 or above. You can find the hg19 " \
			   "index at lobstr.teamerlich.org/download.html", ERROR);
  }
}

std::string GetReadDebug(const ReadPair& read_pair,
                         const std::string& detector_err,
                         const std::string& detector_msg,
                         const std::string& aln_err,
                         const std::string& aln_msg) {
  stringstream msg;
  msg << "[STRDetector]: processing "
      << read_pair.reads[0].ID
      << " motif 1 ";
  if (read_pair.reads[0].repseq.empty()) {
    msg << "NA";
  } else {
    msg << read_pair.reads[0].repseq;
  }
  msg << " " << detector_err << " " << detector_msg << " " << aln_err << " " << aln_msg;
  return msg.str();
}

string GetReadGroup() {
  stringstream read_group;
  read_group << "lobSTR";
  if (!read_group_sample.empty()) {
    read_group << ";" << read_group_sample;
  }
  if (!read_group_library.empty()) {
    read_group << ";" << read_group_library;
  }
  return read_group.str();
}

void TrimRead(const string& input_nucs,
              const string& input_quals,
              string* trimmed_nucs,
              string* trimmed_quals,
              int cutoff) {
  // if last bp is fine, return as is
  size_t l = input_nucs.length();
  if (static_cast<int>(input_quals[l - 1] - QUALITY_CONSTANT)
      >= cutoff) {
    *trimmed_nucs = input_nucs;
    *trimmed_quals = input_quals;
    return;
  }

  // else find the best place to chop
  // done according to bwa manual -q option
  // don't let read length go below minimum
  size_t max_x = min_read_length;
  int max_score = 0;
  for (size_t x = min_read_length; x <= l; x++) {
    int score = 0;
    for (size_t i = x+1; i < l; i++) {
      score += (cutoff-(input_quals[i]-QUALITY_CONSTANT));
    }
    if (score >= max_score) {
      max_score = score;
      max_x = x;
    }
  }
  *trimmed_nucs = input_nucs.substr(0, max_x + 1);
  *trimmed_quals = input_quals.substr(0, max_x + 1);
}

bool fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
}

bool valid_nucleotides_string(const string &str) {
  if (str.empty()) {
    return false;
  }
  return str.find_first_not_of("ACGTNacgtn");
}

double calculate_N_percentage(const std::string& nuc) {
  size_t n_count = 0;
  for (size_t i = 0; i < nuc.length(); i++)
    if (nuc[i] == 'N' || nuc[i] == 'n')
      n_count++;
  return (static_cast<double>(n_count))/
    (static_cast<double>(nuc.length()));
}

string reverseComplement(const string& nucs) {
  string rev;
  size_t size = nucs.size();
  rev.resize(size);
  for (size_t i = 0; i < size; i++) {
    rev.replace(size-i-1, 1, 1, complement(nucs[i]));
  }
  return rev;
}

char complement(const char nucleotide) {
  switch (nucleotide) {
  case 'A':
  case 'a':
    return 'T';
  case 'T':
  case 't':
    return 'A';
  case 'G':
  case 'g':
    return 'C';
  case 'C':
  case 'c':
    return 'G';
  default:
    return 'N';
  }
  return 'N';
}

int nucToNumber(const char& nuc) {
  switch (nuc) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  default:  return 4;
  }
}

bool CheckRepeatCount(const std::string& nucs, const size_t& k, const size_t& minlen, std::string* bestkmer) {
  if (k < 1 || k > 6) {
    PrintMessageDieOnError("Invalid kmer size for CheckRepeatCount", ERROR);
  }
  if (nucs.size() < minlen) {
    PrintMessageDieOnError("Nucs smaller than kmer*count size", ERROR);
  }
  // Get most common kmer
  map<string, int> countKMers;
  std::string subseq;
  subseq.resize(k);
  std::string kmer = "";
  int maxkmer = 0;
  for (size_t i = 0; i < nucs.size() - k; i=i+k) {
    subseq = nucs.substr(i, k);
    countKMers[subseq]++;
    if (countKMers[subseq] > maxkmer) {
      kmer = subseq;
      maxkmer = countKMers[subseq];
    }
  }
  *bestkmer = kmer;
  // Check that we have enough consecutive occurrences of the kmer
  stringstream kmerrep;
  size_t s = 0;
  while (s < minlen) {
    kmerrep << kmer;
    s += kmer.size();
  }
  return nucs.find(kmerrep.str().substr(0, minlen)) != std::string::npos;
}

std::string reverse(const std::string& s) {
  string rev;
  size_t size = s.size();
  rev.resize(size);
  for (size_t i = 0; i < size; i++) {
    rev.replace(size-i-1, 1, s.substr(i, 1));
  }
  return rev;
}

IFileReader* create_file_reader(const string& filename1,
                                const string& filename2) {
  switch (input_type) {
    case INPUT_FASTA:
      if (paired) {
        return new FastaPairedFileReader(filename1, filename2);
      } else {
        if (gzip) {
          return new ZippedFastaFileReader(filename1);
        } else {
          return new FastaFileReader(filename1);
        }
      }
  case INPUT_FASTQ:
    if (paired) {
      return new FastqPairedFileReader(filename1, filename2);
    } else {
      if (gzip) {
        return new ZippedFastqFileReader(filename1);
      } else {
        return new FastqFileReader(filename1);
      }
    }
  case INPUT_BAM:
    if (paired) {
      return new BamPairedFileReader(filename1);
    } else {
      return new BamFileReader(filename1);
    }
  default:
    // This should really never happen
    errx(1, "Internal error, unknown 'input_type' (%d)",
         static_cast<int>(input_type));
  }
}

void GenerateCorrectCigar(CIGAR_LIST* cigar_list,
                          const std::string& nucs,
                          bool* added_s,
                          bool* cigar_had_s) {
  // how many nucleotides the current cigar counts for
  size_t cigar_length = 0;
  *cigar_had_s = false;
  for (size_t i = 0; i < cigar_list->cigars.size(); i++) {
    CIGAR cig = cigar_list->cigars[i];
    if (cig.cigar_type == 'M' || cig.cigar_type == 'I' || cig.cigar_type == 'S') {
      cigar_length += cig.num;
    }
    if (cig.cigar_type == 'S') *cigar_had_s = true;
  }
  // bp not covered by the cigar score
  int diff = nucs.length() - cigar_length;
  if (diff > 0) {
    *added_s = true;
    if (cigar_list->cigars[cigar_list->cigars.size()-1].cigar_type == 'M') {
      cigar_list->cigars[cigar_list->cigars.size()-1].num += diff;
    } else {
      CIGAR cig;
      cig.num = diff;
      cig.cigar_type = 'M';
      cigar_list->cigars.push_back(cig);
      cigar_list->ResetString();
    }
  } else {
    *added_s = false;
  }
}

bool ExtractCigar(const CIGAR_LIST& cigar_list, const int& cigar_start,
		  const int& region_start, const int& region_end,
		  CIGAR_LIST* cigar_subset) {
  int pos = cigar_start; // beginning position of current cigar item
  size_t cigar_start_index = 0; // position in cigar string where region starts
  size_t last_match_index = 0; // position of last 'M' cigar
  int bp = 0;
  char cigar_type; // type of the cigar item

  // *** Determine cigar end *** //
  int cigar_region_length = 0;
  for (size_t i = 0; i < cigar_list.cigars.size(); i++) {
    bp = cigar_list.cigars.at(i).num;
    cigar_type = cigar_list.cigars.at(i).cigar_type;
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D') {
      cigar_region_length += bp;
    }
  }

  // *** Checks for boundaries *** //
  if (cigar_start < 0) { return false; }
  if (region_start < cigar_start || region_start >= cigar_start+cigar_region_length) { return false; }
  if (region_end < cigar_start || region_end >= cigar_start+cigar_region_length) { return false; }
  if (region_end < region_start) { return false; }

  // *** Set start index *** //
  while (pos < region_start && cigar_start_index < cigar_list.cigars.size()) {
    bp = cigar_list.cigars.at(cigar_start_index).num;
    cigar_type = cigar_list.cigars.at(cigar_start_index).cigar_type;
    // If match or del, increment position. All other CIGAR strings (I, P, H, S, N) don't increment
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D') {
      pos += bp;
    }
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X') {
      last_match_index = cigar_start_index;
    }
    cigar_start_index += 1;
  }
  cigar_start_index = last_match_index;
  if (cigar_start_index == 0 && !(cigar_list.cigars.at(cigar_start_index).cigar_type == 'M' ||
				  cigar_list.cigars.at(cigar_start_index).cigar_type == '=' ||
				  cigar_list.cigars.at(cigar_start_index).cigar_type == 'X')) {
    // No match at the left enpd
    return false;
  }
  // *** Set end index *** //
  size_t cigar_end_index = cigar_list.cigars.size() - 1;
  last_match_index = cigar_list.cigars.size() - 1;
  pos = cigar_start + cigar_region_length;
  while (pos > region_end) {
    bp = cigar_list.cigars.at(cigar_end_index).num;
    cigar_type = cigar_list.cigars.at(cigar_end_index).cigar_type;
    // If match or del, decrement position. All other CIGAR strings (I, P, H, S, N) don't decrement
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D') {
      pos -= bp;
    }
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X') {
      last_match_index = cigar_end_index;
    }
    if (cigar_end_index == 0) { break;}
    cigar_end_index -= 1;
  }
  cigar_end_index = last_match_index;
  if (cigar_end_index == cigar_list.cigars.size()-1 && !(cigar_list.cigars.at(cigar_end_index).cigar_type == 'M' ||
							 cigar_list.cigars.at(cigar_end_index).cigar_type == '=' ||
							 cigar_list.cigars.at(cigar_end_index).cigar_type == 'X')) {
    // No match at the right end
    return false;
  }
  cigar_subset->cigars.clear();
  for (size_t i = cigar_start_index; i <= cigar_end_index; i++) {
    cigar_subset->cigars.push_back(cigar_list.cigars.at(i));
  }
  cigar_subset->ResetString();
  return true;
}

int CountOccurrences(const std::string& str, const std::string& substr) {
  int n = 0;
  size_t pos = str.find(substr, 0);
  while (pos != std::string::npos) {
    n = n + 1;
    pos += 1;
    pos = str.find(substr, pos);
  }
  return n;
}

std::string currentDateTime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  return buf;
}

std::string string_replace(std::string src,
                           const std::string& target,
                           const std::string& replace) {
  if (target.length() == 0) {
    // searching for a match to the empty string will result in 
    //  an infinite loop
    //  it might make sense to throw an exception for this case
    return src;
  }
  if (src.length() == 0) {
    return src;  // nothing to match against
  }
  size_t idx = 0;
  for (;;) {
    idx = src.find( target, idx);
    if (idx == string::npos)  break;
    src.replace(idx, target.length(), replace);
    idx += replace.length();
  }
  return src;
}

std::vector<std::string> &split(const std::string &s,
                                char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
}

std::string GetTime() {
  stringstream t;
  time_t et; time(&et);
  t << ctime(&et);
  string tstring = t.str();
  tstring.erase(tstring.find_last_not_of(" \n\r\t")+1);
  return tstring;
}

std::string GetDurationString(const size_t duration)
{
  const size_t days = duration/60/60/24;
  const size_t hours = (duration/60/60)%24;
  const size_t minutes = (duration/60)%60;
  const size_t seconds = duration%60;

  stringstream ss;
  if (days>0)
    ss << days << " days and " ;
  ss << setw(2) << setfill('0') << hours << ':'
     << setw(2) << setfill('0') << minutes << ':'
     << setw(2) << setfill('0') << seconds ;

  return ss.str();
}

// Prints Running time information
void OutputRunningTimeInformation(const size_t start_time,
                                  const size_t processing_start_time,
                                  const size_t end_time,
                                  const size_t num_threads,
                                  const size_t units_processed)
{
  stringstream msg;
  if (num_threads<=0 || (start_time>end_time) || (processing_start_time>end_time)) {
    //Should never happen, but an error is better than invalid output (or division by zero)
    msg << "Internal Error: invalid values for OutputRunningTimeInformation ("
        << "start_time=" << start_time << " processing_start_time="<<processing_start_time
        << "end_time=" << end_time << " num_threads=" << num_threads << ")";
    PrintMessageDieOnError(msg.str(), ERROR);
    return;
  }
  size_t total_seconds = difftime(end_time, start_time);
  size_t processing_seconds = difftime(end_time, processing_start_time);

  msg << "Total Running Time " << GetDurationString(total_seconds);
  PrintMessageDieOnError(msg.str(), PROGRESS);

  msg.str("");
  msg.clear();
  msg << "Processing time: " << GetDurationString(processing_seconds)
      << " (" << processing_seconds << " seconds)";
  PrintMessageDieOnError(msg.str(), PROGRESS);

  msg.str("");
  msg.clear();
  msg << "Processing speed (avg.): ";
  if (processing_seconds>0) {
    double units_seconds = (double)(units_processed) / processing_seconds / num_threads;
    msg << units_seconds ;
  } else {
    msg << "<1";
  }
  msg << " units/seconds/thread";
  PrintMessageDieOnError(msg.str(), PROGRESS);
}
