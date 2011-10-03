/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <cmath>
#include <iostream>
#include <sstream>

#include "CachedWindowGenerator.h"
#include "common.h"
#include "HammingWindowGenerator.h"
#include "runtime_parameters.h"

using namespace std;

HammingWindowGenerator* HammingWindowGenerator::singleton;

//TODO: This one is not thread safe (when it's first called).
//      add a static lock.
//      Current the first creation is done before any threads are created, so it kinda works.
HammingWindowGenerator* HammingWindowGenerator::GetHammingWindowSingleton() {
  if ( singleton == NULL )
    singleton = new HammingWindowGenerator();
  
  return singleton;
}

void HammingWindowGenerator::CalculateWindow(size_t window_size, WindowVector &v) {
  double pi = 3.1416;
  int N = window_size-1;
  
  //Calculate the actual values, store them in "in_v"
  v.resize(window_size);
  for (size_t i = 0 ; i < window_size; ++i ) {
    v[i] = 0.54-0.46*cos(2*pi*i/N);
  }
}
