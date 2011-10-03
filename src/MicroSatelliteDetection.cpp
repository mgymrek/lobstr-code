/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <algorithm>
#include <iostream>

#include "common.h"
#include "MicroSatelliteDetection.h"
#include "runtime_parameters.h"

using namespace std;

MicroSatelliteDetection::MicroSatelliteDetection(const std::vector<double>& _v) :
  v(_v) {}

double MicroSatelliteDetection::find_max_lobe_detection_value() {
  std::vector<double>::const_iterator it = max_element(v.begin(), v.end());
  return *it;
}

bool MicroSatelliteDetection::is_above_threshold() {
  double max = find_max_lobe_detection_value();
  return max>fft_lobe_threshold;
}

void MicroSatelliteDetection::find_start_end(size_t &start, size_t &end) {
  std::vector<double>::const_iterator it = max_element(v.begin(), v.end());
  size_t index_of_max = distance(v.begin(), it);
  
  //Go backwards (to the left)
  int starti;
  for (starti = index_of_max ; starti>=0; starti--)
    if (v[starti] < fft_lobe_threshold)
      break;
  start = starti+1;
  
  //Go forwards (to the right)
  for (end = index_of_max ; end < v.size(); end++)  {
    if ( v[end] < fft_lobe_threshold )
      break;
  }
  end--;
}

bool MicroSatelliteDetection::is_whole_read_microsatellite(size_t start, size_t end) {
  return (start == 0 && end >= (v.size()-1));
}
