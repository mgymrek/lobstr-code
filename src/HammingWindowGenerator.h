/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __HAMMING_WINDOW_GENERATOR_H__
#define __HAMMING_WINDOW_GENERATOR_H__

#include "CachedWindowGenerator.h"

class HammingWindowGenerator : public CachedWindowGenerator {
 public:
  static HammingWindowGenerator* GetHammingWindowSingleton();
  
 protected:
  static HammingWindowGenerator* singleton;
  
 HammingWindowGenerator() :
  CachedWindowGenerator("hamming") {}
  void CalculateWindow(size_t window_size, WindowVector /*output*/ &v) ;
};

#endif /* __HAMMING_WINDOW_GENERATOR_H__ */
