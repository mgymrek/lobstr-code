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

#ifndef SRC_RUNINFO_H_
#define SRC_RUNINFO_H_

#include <sstream>
#include <string>
#include <vector>

/*
  RunInfo stores information about the run that will be uploaded to Amazon for further analysis
 */

class RunInfo {
 public:
  // time at which run started
  std::string starttime;
  // time at which run ended
  std::string endtime;
  // git version
  std::string gitversion;
  // machine type
  std::string machtype;
  // run time params
  std::string params;
  // error msg
  std::string error;

  // Alignment stats
  int num_aligned_reads;
  int num_stitched;
  int num_single;
  int num_mates;
  int num_reverse;
  int total_insert;
  int num_nonunit;

  // Allelotype stats
  std::vector<std::string> samples;
  std::vector<int> num_calls; // vectors are one per sample
  std::vector<int> num_calls5x;
  std::vector<int> total_coverage;
  std::vector<int> total_agree;
  std::vector<std::vector<int> > calltype_by_period; // all samples

  // Initialize values
  void Reset() {
    starttime = "";
    endtime = "";
    gitversion = "";
    machtype = "";
    params = "";
    error = "";
    num_aligned_reads = 0;
    num_stitched = 0;
    num_single = 0;
    num_mates = 0;
    num_reverse = 0;
    total_insert = 0;
    num_nonunit = 0;
    samples.clear();
    num_calls.clear();
    num_calls5x.clear();
    total_coverage.clear();
    total_agree.clear();
    calltype_by_period.clear();
  }
  
  // Print to string
  std::string PrintToString(int program) {
    std::stringstream ss;
    ss << "Starttime: " << starttime << std::endl;
    ss << "Endtime: " << endtime << std::endl;
    ss << "Git version: " << gitversion << std::endl;
    ss << "Machtype: " << machtype << std::endl;
    ss << "Params: " << params << std::endl;
    if (!error.empty()) {
      ss << "Error: " << error << std::endl;
    }
    // Print quality metrics
    if (program == 0) {
      if (num_aligned_reads > 0) {
	ss << "Alignment stats" << std::endl;
	ss << "Reads aligned\t" << num_aligned_reads << std::endl;
	ss << "Perc. stitched\t" << static_cast<float>(num_stitched)/static_cast<float>(num_aligned_reads) << std::endl;
	ss << "Perc. single\t" << static_cast<float>(num_single)/static_cast<float>(num_aligned_reads) << std::endl;
	ss << "Perc. mate\t" << static_cast<float>(num_mates)/static_cast<float>(num_aligned_reads) << std::endl;
	ss << "Perc. reverse\t" << static_cast<float>(num_reverse)/static_cast<float>(num_single+num_stitched) << std::endl;
	ss << "Perc. nonunit\t" << static_cast<float>(num_nonunit)/static_cast<float>(num_single+num_stitched) << std::endl;
	if (num_mates > 0) {
	  ss << "Mean insert size: " << static_cast<float>(total_insert)/static_cast<float>(num_mates) << std::endl;
	}
      }
    } else {
      ss << "Allelotype stats" << std::endl;
      for (size_t i = 0; i < samples.size(); i++) {
	ss << "Sample\t" << samples.at(i) << std::endl;
	ss << "Num calls\t" << num_calls.at(i) << std::endl;
	ss << "Num calls >=5x\t" << num_calls5x.at(i) << std::endl;
	ss << "Mean coverage\t" << static_cast<float>(total_coverage.at(i))/static_cast<float>(num_calls.at(i)) << std::endl;
	ss << "Mean perc. agree\t" << static_cast<float>(total_agree.at(i))/static_cast<float>(num_calls.at(i)) << std::endl;
	ss << std::endl;
      }
      ss << "Call type by period (0/0, 0/1, 1/1, 1/2)" << std::endl;
      for (size_t i = 0; i < calltype_by_period.size(); i++) {
	int period = static_cast<int>(i)+1;
	ss << "Period " << period;
	for (size_t j = 0; j < calltype_by_period.at(i).size(); j++) {
	  ss << "\t" << calltype_by_period.at(i).at(j);
	}
	ss << std::endl;
      }
    }

    return ss.str();
  }
};

#endif  // SRC_RUNINFO_H_
