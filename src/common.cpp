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
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 

#include <map>
#include <iostream>
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
  string stats_string = run_info.PrintToString();
  sWriter.Write(stats_string);
  // Upload to AWS S3
  if (!noweb) {
    // Set upt POST data
    int size = stats_string.size();
    stringstream pd;
    pd << "POST /test.py HTTP/1.1\r\n";
    pd << "Host: mgymrek.scripts.mit.edu\r\n";
    pd << "Connection: Keep-Alive\r\n";
    pd << "Content-Type: text/plain, name=\"lobSTRstats\"\r\n";
    pd << "Content-Length: " << size;
    pd << "\r\n\r\n";
    string post_header = pd.str();
    // PrintMessageDieOnError(pd.str(), DEBUG);
    struct hostent *server;
    server = gethostbyname("mgymrek.scripts.mit.edu");
    if (server == NULL) {
      PrintMessageDieOnError("Server not found", ERROR);
    }
    struct sockaddr_in serverAddr;
    memset(&serverAddr, 0, sizeof(serverAddr));
    serverAddr.sin_family = AF_INET;
    memcpy(&serverAddr.sin_addr.s_addr,
           server->h_addr, 
           server->h_length);
    serverAddr.sin_port = htons(80);
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
      PrintMessageDieOnError("Couldn't initialize socket", ERROR);
    }
    if (connect(sock, (struct sockaddr *) &serverAddr, sizeof(serverAddr)) < 0) {
      PrintMessageDieOnError("Failed internet connection", ERROR);
    }
    if (send(sock, post_header.c_str(), post_header.size(),0) < 0) {
      PrintMessageDieOnError("Couldn't send header", ERROR);
    }
    if (send(sock, stats_string.c_str(), size, 0) < 0) {
      PrintMessageDieOnError("Couldn't send data", ERROR);
    }
    /*char buf[1024];
    if (recv(sock, buf, 1024, 0) < 0) {
      PrintMessageDieOnError("Error receiving", DEBUG);
      }
      cerr << buf << endl;*/
    close(sock);
  }
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
  cerr << "[" << (program == LOBSTR ? "lobSTR":"allelotype")
       << "-" << _GIT_VERSION << "] " << typestring << msg << endl;
  if (msgtype == ERROR) {
    run_info.error = msg;
    run_info.endtime = GetTime();
    OutputRunStatistics();
    exit(1);
  }
}

void TrimRead(const string& input_nucs,
              const string& input_quals,
              string* trimmed_nucs,
              string* trimmed_quals,
              int cutoff) {
  // if last bp is fine, return as is
  size_t l = input_nucs.length();
  if (static_cast<int>(input_quals.at(l - 1) - QUALITY_CONSTANT)
      >= cutoff) {
    *trimmed_nucs = input_nucs;
    *trimmed_quals = input_quals;
    return;
  }

  // else find the best place to chop
  // done according to bwa manual -q option
  // don't let read length go below minimum
  size_t max_x;
  int max_score = 0;
  for (size_t x = min_read_length; x <= l; x++) {
    int score = 0;
    for (size_t i = x+1; i < l; i++) {
      score += (cutoff-(input_quals.at(i)-QUALITY_CONSTANT));
    }
    if (score >= max_score) {
      max_score = score;
      max_x = x;
    }
  }
  *trimmed_nucs = input_nucs.substr(0, max_x + 1);
  *trimmed_quals = input_quals.substr(0, max_x + 1);
}

size_t count(const string& s, const char& c) {
  size_t num = 0;
  for (size_t i = 0; i < s.length(); i++) {
    if (s.at(i) == c) num++;
  }
  return num;
}


bool getMSSeq(const string& nucs, int k, string* repeat) {
  if (k < 1 || k > 6) return false;
  if (static_cast<int>(nucs.size()) < k) return false;
  map<string, int> countKMers;
  size_t i;
  string subseq;
  string kmer;
  int maxkmer = 0;
  subseq.resize(k);
  for (i = 0; i < nucs.size() - k; i++) {
    subseq = nucs.substr(i, k);
    countKMers[subseq]++;
    if (countKMers.at(subseq) > maxkmer) {
      kmer = subseq;
      maxkmer = countKMers.at(subseq);
    }
  }
  string repseqfw;
  string repseqrev;
  string repseq;
  if (kmer.size() < 1 || kmer.size() > 6) return false;
  getCanonicalMS(kmer, &repseqfw);
  getCanonicalMS(reverseComplement(repseqfw), &repseqrev);
  repseq = getFirstString(repseqfw, repseqrev);
  *repeat = repseq;
  return true;
}

string getFirstString(const std::string& seq1, const std::string& seq2) {
  for (size_t i = 0; i < seq1.size(); i++) {
    if (nucToNumber(seq1[i]) < nucToNumber(seq2[i])) return seq1;
    if (nucToNumber(seq1[i]) > nucToNumber(seq2[i])) return seq2;
  }
  return seq1;
}

bool IsPerfectRepeat(const std::string& sequence,
                     const std::string& repeat) {
  // find first occurrence of repeat
  size_t found;
  found = sequence.find(repeat);

  if (found == string::npos || found > repeat.length() - 1) return false;
  // check the part before found
  if (sequence.substr(0, found) != repeat.substr(repeat.length() - found,
                                                 found)) return false;
  for (size_t i = found; i < sequence.length()
         - repeat.length() + 1; i += repeat.length()) {
    string test_seq = sequence.substr(i, repeat.length());
    if (test_seq != repeat) return false;
  }
  // check the part after
  return true;
}

float GetAverageQualityScore(const vector<string>& qualities) {
  if (qualities.size() == 0) { return 0;}
  float average_quality = 0;
  for (vector<string>::const_iterator it = qualities.begin();
       it != qualities.end(); ++it) {
    average_quality += GetQualityScore(*it);
  }
  return average_quality/qualities.size();
}

float GetQualityScore(const std::string& quality_score) {
  if (quality_score.length() == 0) return 0;
  float total_quality = 0;
  for (size_t i = 0; i < quality_score.length(); ++i) {
    int qs = quality_score.at(i);
    total_quality += qs - 33;
  }
  return total_quality/quality_score.size();
}

int GetChromNumber(string chromosome) {
  // get whatever is after "chr"
  string chrom_string = chromosome.substr(3);
  // convert this to a number
  if (chrom_string == "X") {
    return 23;
  } else if (chrom_string == "Y") {
    return 24;
  } else {
    return atoi(chrom_string.c_str());
  }
}

bool fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
}

bool valid_nucleotides_string(const string &str) {
  if (str.empty())
    return false;
  for (size_t i = 0 ; i < str.length(); ++i) {
    const char ch = str[i];
    if ( (ch != 'A') && (ch != 'C') && (ch != 'G') && (ch != 'T') &&
         (ch != 'N') &&
         (ch != 'a') && (ch != 'c') && (ch != 'g') && (ch != 't') &&
         (ch != 'n') )
      return false;
  }
  return true;
}

char OneAbundantNucleotide(const std::string& nuc, float perc_threshold) {
  size_t countA = 0, countC = 0, countG = 0, countT = 0;
  for (size_t i = 0; i < nuc.length(); i++) {
    switch (nuc.at(i)) {
      case 'A':
      case 'a':
        countA++;
      break;
      case 'C':
      case 'c':
        countC++;
      break;
      case 'G':
      case 'g':
        countG++;
      break;
      case 'T':
      case 't':
        countT++;
      break;
      case 'N':
      case 'n':
        break;
      default:
        errx(1, "Internal error: OneAbundantNucleotide " \
             "called with invalid nucleotide string '%s'" \
             ", characte '%c'", nuc.c_str(), nuc.at(i));
      }
  }
  size_t threshold = nuc.length()*perc_threshold;
  if (countA >= threshold)
    return 'A';
  if (countC >= threshold)
    return 'C';
  if (countG >= threshold)
    return 'G';
  if (countT >= threshold)
    return 'T';
  return 0;
}

double calculate_N_percentage(const std::string& nuc) {
  size_t n_count = 0;
  for (size_t i = 0; i < nuc.length(); i++)
    if (nuc.at(i) == 'N' || nuc.at(i) == 'n')
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

std::string reverse(const std::string& s) {
  string rev;
  size_t size = s.size();
  rev.resize(size);
  for (size_t i = 0; i < size; i++) {
    rev.replace(size-i-1, 1, s.substr(i, 1));
  }
  return rev;
}

void getCanonicalMS(const string& msnucs, string* canonical) {
  // common ones
  // first check to see if it is hashed already
  if (canonicalMSTable.find(msnucs) != canonicalMSTable.end()) {
    *canonical = canonicalMSTable.at(msnucs);
    return;
  }
  string newseq;
  size_t size = msnucs.size();
  size_t i;
  size_t j;
  *canonical = msnucs;
  newseq.resize(size);
  for (i = 1; i < size ; i++) {
    newseq = msnucs.substr(size-i, size) + msnucs.substr(0, size-i);
    // if newsq > canon, make it canon
    for (j = 0; j < size; j++) {
      if (nucToNumber(newseq[j]) < nucToNumber((*canonical)[j])) {
        *canonical = newseq;
        break;
      } else if (nucToNumber(newseq[j]) > nucToNumber((*canonical)[j])) {
        break;
      }
    }
  }
  canonicalMSTable.insert(pair<string, string>
                          (msnucs, *canonical));
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

string GenerateS3Command(const string& bucket,
                         const string& filename,
                         const string& configfile) {
  stringstream s;
  s << "s3cmd -c " << configfile << " get --skip-existing "
    << bucket << "/" << filename << " "
    << "/mnt/lobstr/" << filename << endl;
  return s.str();
}

void GenerateCorrectCigar(CIGAR_LIST* cigar_list,
                          const std::string& nucs,
                          bool* added_s,
                          bool* cigar_had_s) {
  // how many nucleotides the current cigar counts for
  size_t cigar_length = 0;
  *cigar_had_s = false;
  for (int i = 0; i < cigar_list->cigars.size(); i++) {
    CIGAR cig = cigar_list->cigars.at(i);
    if (cig.cigar_type == 'M' || cig.cigar_type == 'I' || cig.cigar_type == 'S') {
      cigar_length += cig.num;
    }
    if (cig.cigar_type == 'S') *cigar_had_s = true;
  }
  // bp not covered by the cigar score
  int diff = nucs.length() - cigar_length;
  if (diff > 0) {
    *added_s = true;
    if (cigar_list->cigars.at(cigar_list->cigars.size()-1).cigar_type == 'M') {
      cigar_list->cigars.at(cigar_list->cigars.size()-1).num += diff;
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

std::string currentDateTime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
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

std::string fftw_complex_to_string(fftw_complex v) {
  stringstream s;
  s.setf(ios::fixed, ios::floatfield);
  s.width(7);
  s.precision(5);
  s << v[0] << " ";
  s << ((v[1] >= 0) ? "+ " : "- ");
  s << (std::abs(v[1])) << "i";
  return s.str();
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
