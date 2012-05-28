/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_STRDETECTOR_H__
#define SRC_STRDETECTOR_H__

#include "src/ISatellite.h"
#include "src/MSReadRecord.h"

class STRDetector : public ISatellite {
 public:
  STRDetector();
  bool ProcessReadPair(ReadPair* read_pair);
 private:
  bool ProcessRead(MSReadRecord* read);
};

#endif  // SRC_STRDETECTOR_H__

