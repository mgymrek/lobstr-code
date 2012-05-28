/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_CACHEDWINDOWGENERATOR_H__
#define SRC_CACHEDWINDOWGENERATOR_H__

#include <pthread.h>
#include <tr1/unordered_map>

#include <string>
#include <vector>

typedef std::vector<double> WindowVector;

class CachedWindowGenerator {
 public:
  virtual ~CachedWindowGenerator();
  const WindowVector* GetWindow(size_t window_size);

 protected:
  pthread_mutex_t access_lock;
  typedef std::tr1::unordered_map<size_t, WindowVector> WINDOWS_HASH;
  WINDOWS_HASH cached_windows;
  const std::string window_function_name;
  explicit CachedWindowGenerator(const std::string& name);
  const WindowVector* CreateWindow(size_t window_size);

  void lock();
  void unlock();

  virtual void CalculateWindow(size_t window_size,
                               WindowVector &v) = 0;
};


#endif  // SRC_CACHEDWINDOWGENERATOR_H__
