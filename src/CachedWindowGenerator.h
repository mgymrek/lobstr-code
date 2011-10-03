/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __CACHED_WINDOW_GENERATOR_H__
#define __CACHED_WINDOW_GENERATOR_H__

#include <pthread.h>
#include <tr1/unordered_map>
#include <vector>

using namespace std;

typedef vector<double> WindowVector;

class CachedWindowGenerator {
 public:
  virtual ~CachedWindowGenerator();
  const WindowVector* GetWindow(size_t window_size);
  
 protected:
  pthread_mutex_t access_lock;
  typedef std::tr1::unordered_map<size_t, WindowVector> WINDOWS_HASH;
  WINDOWS_HASH cached_windows;
  const string window_function_name;

  CachedWindowGenerator(const string& name);
  
  const WindowVector* CreateWindow(size_t window_size);
  
  void lock();
  void unlock();

  virtual void CalculateWindow(size_t window_size, WindowVector /*output*/ &v) = 0 ;
};


#endif /* __CACHED_WINDOW_GENERATOR_H__ */
