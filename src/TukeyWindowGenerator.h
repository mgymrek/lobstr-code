/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

*/

#ifndef SRC_TUKEYWINDOWGENERATOR_H__
#define SRC_TUKEYWINDOWGENERATOR_H__

#include "src/CachedWindowGenerator.h"

class TukeyWindowGenerator : public CachedWindowGenerator {
 public:
  static TukeyWindowGenerator* GetTukeyWindowSingleton();

 protected:
  static TukeyWindowGenerator* singleton;

  TukeyWindowGenerator() :
  CachedWindowGenerator("tukey") {}
  void CalculateWindow(size_t window_size, WindowVector &v);
};

#endif  // SRC_TUKEYWINDOWGENERATOR_H__
