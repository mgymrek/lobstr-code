/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_HAMMINGWINDOWGENERATOR_H__
#define SRC_HAMMINGWINDOWGENERATOR_H__

#include "src/CachedWindowGenerator.h"

class HammingWindowGenerator : public CachedWindowGenerator {
 public:
  static HammingWindowGenerator* GetHammingWindowSingleton();

 protected:
  static HammingWindowGenerator* singleton;

  HammingWindowGenerator() :
  CachedWindowGenerator("hamming") {}
  void CalculateWindow(size_t window_size, WindowVector& v);
};

#endif  // SRC_HAMMINGWINDOWGENERATOR_H__
