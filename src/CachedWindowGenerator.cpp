/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <cmath>
#include <iostream>
#include <sstream>

#include "CachedWindowGenerator.h"
#include "common.h"
#include "runtime_parameters.h"

using namespace std;

CachedWindowGenerator::CachedWindowGenerator(const std::string &name) :
  window_function_name(name) {
  pthread_mutex_init(&access_lock, NULL);
}

CachedWindowGenerator::~CachedWindowGenerator() {
  pthread_mutex_destroy(&access_lock);
}

void CachedWindowGenerator::lock() {
  pthread_mutex_lock(&access_lock);
}

void CachedWindowGenerator::unlock() {
  pthread_mutex_unlock(&access_lock);
}

const WindowVector* CachedWindowGenerator::CreateWindow(size_t window_size) {
  WindowVector in_v;

  CalculateWindow(window_size, in_v);

  std::pair<WINDOWS_HASH::iterator, bool> result =
    cached_windows.insert(std::pair<size_t, WindowVector>(window_size, in_v));

  const WindowVector &out_v = result.first->second; //gotta love STL :)
  const WindowVector *p = &out_v;
  return p;
}

const WindowVector* CachedWindowGenerator::GetWindow(size_t window_size) {
  lock();

  // If we already calculated this window, just return it
  WINDOWS_HASH::const_iterator it = cached_windows.find(window_size);
  if (it != cached_windows.end()) {
    const WindowVector& v = it->second;
    const WindowVector *p = &v;
    unlock();
    return p;
  }
  
  // Not found, we need to create it (this will also store it in the cache)
  
  const WindowVector *p = CreateWindow(window_size);
  
  unlock();
  return p;
}
