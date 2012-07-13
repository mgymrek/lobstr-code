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

#include <math.h>

#include <algorithm>
#include <string>
#include <vector>

#include "src/common.h"
#include "src/EntropyDetection.h"
#include "src/runtime_parameters.h"

using namespace std;


static double MinusPlogP(double value) {
  return -1 * value * log(value)/log(2);
}


EntropyDetection::EntropyDetection(const string& nucleotides,
                                   int size, int step) {
  _nucs = nucleotides;
  _window_size = size;
  _window_step = step;
}

EntropyDetection::~EntropyDetection() {}

double EntropyDetection::EntropyOneWindowDinuc(const std::string& window_nucs) {
  size_t window_length = window_nucs.length();
  vector<int> kmer_counts(44, 0);
  float subseqs = 0;
  for (size_t i = 0; i < window_length - 2; ++i) {
    char nuc1 = window_nucs.at(i);
    char nuc2 = window_nucs.at(i+1);
    if (nuc1 != 'N' && nuc2 != 'N') {
      kmer_counts[nucToNumber(nuc1) + nucToNumber(nuc2)*10] += 1;
      subseqs+= 1;
    }
  }
  float entropy = 0;
  for (size_t i = 0; i < kmer_counts.size(); ++i) {
    int count = kmer_counts.at(i);
    if (count !=0) {
      float p = static_cast<float>(count)/subseqs;
      entropy += MinusPlogP(p);
    }
  }
  return (4-entropy)/4;
}

void EntropyDetection::CalculateEntropyWindow() {
  _entropy_window.clear();
  string window_nucleotides = "";
  for (size_t i = 0; i < _nucs.length() - _window_size; i += _window_step) {
    window_nucleotides = _nucs.substr(i, _window_size);
    double entropy = EntropyOneWindowDinuc(window_nucleotides);
    if (entropy > 0.8)
      entropy = 0;
    _entropy_window.push_back(entropy);
  }
}

bool EntropyDetection::EntropyIsAboveThreshold() {
  CalculateEntropyWindow();
  return *max_element(_entropy_window.begin(), _entropy_window.end()) >
    entropy_threshold;
}

void EntropyDetection::FindStartEnd(size_t* start, size_t* end,
                                    bool* repetitive_end) {
  const vector<double> entropy_window = _entropy_window;
  vector<double>::const_iterator it = max_element(entropy_window.begin(),
                                                  entropy_window.end());
  size_t index_of_max = distance(entropy_window.begin(), it);
  float max_entropy = *it;
  // Go backwards (to the left)
  int starti;
  for (starti = index_of_max ; starti >= 0; starti--)
    if (entropy_window[starti] < 0.8*max_entropy)
      break;
  *start = starti+1;
  if (*start == 0) {
    *start += 1;
    *repetitive_end = true;
  }
  // Go forwards (to the right)
  int endi;
  for (endi = index_of_max ; endi < entropy_window.size(); endi++)  {
    if (entropy_window[endi] < 0.8*max_entropy)
      break;
  }
  (endi)--;
  if (endi == entropy_window.size() -1) {
    *end = endi-1;
    *repetitive_end = true;
  } else {
    *end = endi;
    *repetitive_end = false;
  }
}

