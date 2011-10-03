/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __MICROSATELLITE_DETECTION_H__
#define __MICROSATELLITE_DETECTION_H__

#include <vector>

class MicroSatelliteDetection {
 public:
  MicroSatelliteDetection(const std::vector<double>& _v);
  
  double find_max_lobe_detection_value();
  bool is_above_threshold();
  
  void find_start_end(size_t /*output*/ &start, size_t /*output*/ &end);
  bool is_whole_read_microsatellite(size_t start, size_t end);

 private:
  const std::vector<double> &v;
};

#endif /* __MICROSATELLITE_DETECTION_H__ */
