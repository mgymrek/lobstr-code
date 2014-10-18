/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#ifndef SRC_STRDETECTOR_H__
#define SRC_STRDETECTOR_H__

#include "src/MSReadRecord.h"

class STRDetector {
 public:
  STRDetector();
  bool ProcessReadPair(ReadPair* read_pair, std::string* err, std::string* messages);
 private:
  bool ProcessRead(MSReadRecord* read, std::string* err);
};

#endif  // SRC_STRDETECTOR_H__

