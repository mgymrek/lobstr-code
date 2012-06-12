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

#include <cmath>

#include "src/CachedWindowGenerator.h"
#include "src/common.h"
#include "src/HammingWindowGenerator.h"
#include "src/runtime_parameters.h"

HammingWindowGenerator* HammingWindowGenerator::singleton;

// TODO(mgymrek): This one is not thread safe (when it's first called).
// add a static lock.
// Current the first creation is done before any threads are created,
// so it kinda works.
HammingWindowGenerator* HammingWindowGenerator::GetHammingWindowSingleton() {
  if (singleton == NULL)
    singleton = new HammingWindowGenerator();

  return singleton;
}

void HammingWindowGenerator::CalculateWindow(size_t window_size,
                                             WindowVector &v) {
  double pi = 3.1416;
  int N = window_size-1;

  // Calculate the actual values, store them in "in_v"
  v.resize(window_size);
  for (size_t i = 0 ; i < window_size; ++i) {
    v[i] = 0.54-0.46*cos(2*pi*i/N);
  }
}
