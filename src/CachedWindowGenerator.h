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

#ifndef SRC_CACHEDWINDOWGENERATOR_H__
#define SRC_CACHEDWINDOWGENERATOR_H__

#include "config.h"
#include <pthread.h>

#include <string>
#include <vector>

typedef std::vector<double> WindowVector;

#if defined HAVE_CXX11
  #include <unordered_map>
  typedef std::unordered_map<size_t, WindowVector> WINDOWS_HASH;
#elif defined STDCXX_TR1_HEADERS
  #include <tr1/unordered_map>
  typedef std::tr1::unordered_map<size_t, WindowVector> WINDOWS_HASH;
#else
  /* Fall back to regular map */
  #include <map>
  typedef std::map<size_t, WindowVector> WINDOWS_HASH_TYPE;
#endif


class CachedWindowGenerator {
 public:
  virtual ~CachedWindowGenerator();
  const WindowVector* GetWindow(size_t window_size);

 protected:
  pthread_mutex_t access_lock;
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
