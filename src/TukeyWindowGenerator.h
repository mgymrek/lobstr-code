/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

*/

#ifndef __TUKEY_WINDOW_GENERATOR_H__
#define __TUKEY_WINDOW_GENERATOR_H__

#include "CachedWindowGenerator.h"

class TukeyWindowGenerator : public CachedWindowGenerator {
 public:
  static TukeyWindowGenerator* GetTukeyWindowSingleton();

 protected:
  static TukeyWindowGenerator* singleton;
  
 TukeyWindowGenerator() :
  CachedWindowGenerator("tukey")
    {
    }
  
  void CalculateWindow(size_t window_size, WindowVector /*output*/ &v) ;
};


#endif
