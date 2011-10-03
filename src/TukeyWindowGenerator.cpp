/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <cmath>
#include <iostream>
#include <sstream>

#include "CachedWindowGenerator.h"
#include "common.h"
#include "runtime_parameters.h"
#include "TukeyWindowGenerator.h"

using namespace std;

TukeyWindowGenerator* TukeyWindowGenerator::singleton;

//TODO: This one is not thread safe (when it's first called).
//      add a static lock.
//      Current the first creation is done before any threads are created, so it kinda works.
TukeyWindowGenerator* TukeyWindowGenerator::GetTukeyWindowSingleton() {
  if ( singleton == NULL )
    singleton = new TukeyWindowGenerator();
  
  return singleton;
}

void TukeyWindowGenerator::CalculateWindow(size_t window_size, WindowVector &v) {
  double a = tukey_alpha;
  double pi = 3.1416;
  
  //Calculate the actual values, store them in "in_v"
  v.resize(window_size);
  int N = window_size -1;
  for (size_t i = 0 ; i < window_size; ++i ) {
    if (i <= a*N/2){
      v[i] = 0.5*(1+cos(pi*(2*i/(a*N)-1)));
    }
    else if( i <= N*(1-a/2)){
      v[i] = 1;  
    }
    else{
      v[i] = 0.5*(1+cos(pi*(2*i/(a*N)+1-2/a)));
    }
  }
}
