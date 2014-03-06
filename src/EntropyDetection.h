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

#ifndef SRC_ENTROPYDETECTION_H__
#define SRC_ENTROPYDETECTION_H__

#include <string>
#include <vector>

class EntropyDetection {
 private:
  std::string _nucs;
  std::vector<double> _entropy_window;
  int _window_size;
  int _window_step;

 public:
  EntropyDetection(const std::string& nucleotides, int size, int step);
  ~EntropyDetection();
  // dinuc entropy, not generalized but optimized
  double EntropyOneWindowDinuc(const std::string& window_nucs);

  // vector of entropy values for each window
  void CalculateEntropyWindow();

  // determine if there is a window above the entropy threshold
  bool EntropyIsAboveThreshold();

  // determine start and end region
  void FindStartEnd(size_t* start, size_t* end, bool* repetitive_end);
  
  // return the max entropy
  float GetMaxEntropy();
};

#endif  // SRC_ENTROPYDETECTION_H__
