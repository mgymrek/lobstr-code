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

#include <algorithm>

#include "src/common.h"
#include "src/EntropyDetection.h"
#include "src/runtime_parameters.h"
#include "src/STRDetector.h"

using namespace std;

const size_t EXTEND_FLANK = 6;
const size_t MIN_REP_LENGTH[] = {8, 8, 8, 10, 10, 10};

STRDetector::STRDetector() {}

bool STRDetector::ProcessReadPair(ReadPair* read_pair, string* err, string* messages) {
  *err = "Detection-errors-here:";
  *messages = "Detection-notes-here:";
  read_pair->read1_passed_detection = false;
  read_pair->read2_passed_detection = false;
  if (ProcessRead(&read_pair->reads.at(0), err)) {
    read_pair->read1_passed_detection = true;
  }
  // Returns true if at least one read in the pair is detected
  if (read_pair->reads.at(0).paired) {
    if (ProcessRead(&read_pair->reads.at(1), err)) {
      read_pair->read2_passed_detection = true;
    }
    return (read_pair->read1_passed_detection ||
            read_pair->read2_passed_detection);
  } else {
    return read_pair->read1_passed_detection;
  }
}

bool STRDetector::ProcessRead(MSReadRecord* read, string* err) {
  // Get the size of the read so we don't keep computing
  size_t read_length = read->nucleotides.size();
  // Preprocessing checks
  if (read_length < (fft_window_size-1) ||
      calculate_N_percentage(read->nucleotides) > percent_N_discard) {
    if (debug) {
      *err += "failed-read-length-check";
    }
    return false;
  }

  //  Step 1 - sliding window
  size_t nuc_start;
  size_t nuc_end;
  string detected_microsatellite_nucleotides = "";
  EntropyDetection ed_filter(read->nucleotides,
                             fft_window_size, fft_window_step);
  if (!ed_filter.EntropyIsAboveThreshold()) {
    if (debug) {
      stringstream msg;
      msg << "failed-entropy-threshold-" << ed_filter.GetMaxEntropy();
      *err += msg.str();
    }
    return false;
  }
  size_t start, end;
  bool rep_end = false;
  ed_filter.FindStartEnd(&start, &end, &rep_end);

  nuc_start = start * fft_window_step;
  nuc_end = (end + 2) * fft_window_step;

  if (nuc_start >= read_length ||
      nuc_start <= 0 ||
      nuc_end-nuc_start + 1 <= 0 ||
      nuc_start >= nuc_end ||
      nuc_end-nuc_start+1 >= read_length ||
      nuc_end >= read_length) {
    if (debug) {
      *err += "failed-STR-region-location-sanity-check";
    }
    return false;
  }

  // Step 2 - Check if evidence of repeats
  bool passed_repeat_check = false;
  for (size_t i = 1; i <= 6; i++) {
    std::string bestkmer;
    if (CheckRepeatCount(read->nucleotides, i, MIN_REP_LENGTH[i-MIN_PERIOD], &bestkmer)) {
      passed_repeat_check = true;
      break;
    }
  }
  if (!passed_repeat_check) {
    if (debug) {
      *err += "failed-repeat-check";
    }
    return false;
  }

  // Step 3 - Store values in the Read Record
  read->ms_start = nuc_start;
  read->ms_end   = nuc_end;

  // set indices of left, STR, and right regions
  if ((read->ms_start >=
       static_cast<int>(read_length)) ||
      (read->ms_end+1 >=
       static_cast<int>(read_length))) {
    if (debug) {
      *err += "failed-second-STR-region-location-sanity-check;";
    }
    return false;
  }

  // allow more nucleotides in detection step
  if ((EXTEND_FLANK <= nuc_start) &&
      (nuc_start-EXTEND_FLANK < read_length) &&
      (nuc_end + EXTEND_FLANK + 1 < read_length)) {
    detected_microsatellite_nucleotides =
      read->nucleotides.substr(nuc_start-EXTEND_FLANK,
                               nuc_end - nuc_start+1+2*EXTEND_FLANK);
  } else {
    detected_microsatellite_nucleotides =
      read->nucleotides.substr(nuc_start, nuc_end - nuc_start+1);
  }

  read->left_flank_nuc = (nuc_start>0) ?
    read->nucleotides.substr(0, read->ms_start) : "-";
  read->right_flank_nuc = read->nucleotides.substr(read->ms_end+1);
  read->detected_ms_region_nuc = detected_microsatellite_nucleotides;

  // adjust for max flank region lengths, if repetitive end, don't trim
  read->left_flank_index_from_start = 0;
  read->right_flank_index_from_end = 0;
  if ((read->left_flank_nuc.size() > max_flank_len)) {
    string left_flank = read->left_flank_nuc;
    read->left_flank_nuc = left_flank.substr(left_flank.length()-
                                             max_flank_len, max_flank_len);
    read->left_flank_index_from_start = left_flank.length() -
      max_flank_len;
  }
  if ((read->right_flank_nuc.size() > max_flank_len)) {
    string right_flank = read->right_flank_nuc;
    read->right_flank_nuc = right_flank.substr(0, max_flank_len);
    read->right_flank_index_from_end = right_flank.length() - max_flank_len;
  }

  size_t nuc_len = (read->orig_nucleotides.length() -
                    read->right_flank_index_from_end -
                    read->left_flank_index_from_start);

  if (((read->left_flank_index_from_start+nuc_len) <= read_length)) {
    read->nucleotides = read->
      orig_nucleotides.substr(read->left_flank_index_from_start, nuc_len);
    read_length = read->nucleotides.size();
    read->quality_scores = read->
      quality_scores.substr(read->left_flank_index_from_start,
                            read_length);
  } else {
    if (debug) {
      *err += "failed-after-trimming-flanking-regions;";
    }
    return false;
  }
  if ( read->left_flank_nuc.length() < min_flank_len ||
       read->right_flank_nuc.length() < min_flank_len) {
    if (debug) {
      *err += "failed-min-flank-len;";
    }
    return false;
  }
  return true;
}

