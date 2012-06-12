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

#include <err.h>

#include <string>

#include "src/common.h"
#include "src/BamPairedFileReader.h"

BamPairedFileReader::BamPairedFileReader(const std::string& _filename) {
  if (!reader.Open(_filename)) {
    errx(1, "Could not open bam file");
  }
}

bool BamPairedFileReader::GetNextRecord(ReadPair* read_pair) {
  // TODO(mgymrek)
  errx(1, "Paired end bam file input not yet implemented.");
  return false;
}

// Does not apply to paired read readers
bool BamPairedFileReader::GetNextRead(MSReadRecord* read) {
  return false;
}
