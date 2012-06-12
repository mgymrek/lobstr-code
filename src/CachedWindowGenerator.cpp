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
#include <string>
#include <utility>

#include "src/CachedWindowGenerator.h"
#include "src/common.h"
#include "src/runtime_parameters.h"

using namespace std;

CachedWindowGenerator::CachedWindowGenerator(const std::string &name)
  : window_function_name(name) {
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

  const WindowVector &out_v = result.first->second;
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
