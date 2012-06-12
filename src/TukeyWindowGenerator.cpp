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
#include <sstream>

#include "src/CachedWindowGenerator.h"
#include "src/common.h"
#include "src/runtime_parameters.h"
#include "src/TukeyWindowGenerator.h"

using namespace std;

TukeyWindowGenerator* TukeyWindowGenerator::singleton;

// TODO(mgymrek): This one is not thread safe (when it's first called).
// add a static lock.
// Current the first creation is done before any
// threads are created, so it kinda works.
TukeyWindowGenerator* TukeyWindowGenerator::GetTukeyWindowSingleton() {
  if ( singleton == NULL )
    singleton = new TukeyWindowGenerator();
  return singleton;
}

void TukeyWindowGenerator::CalculateWindow(size_t window_size,
                                           WindowVector &v) {
  double a = tukey_alpha;
  double pi = 3.1416;

  // Calculate the actual values, store them in "in_v"
  v.resize(window_size);
  int N = window_size -1;
  for (size_t i = 0 ; i < window_size; ++i) {
    if (i <= a*N/2) {
      v[i] = 0.5*(1+cos(pi*(2*i/(a*N)-1)));
    } else if (i <= N*(1-a/2)) {
      v[i] = 1;
    } else {
      v[i] = 0.5*(1+cos(pi*(2*i/(a*N)+1-2/a)));
    }
  }
}
