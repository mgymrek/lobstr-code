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

#ifndef SRC_REFERENCEGENOME_H__
#define SRC_REFERENCEGENOME_H__

#include <string>

class ReferenceGenome {
 public:
  ReferenceGenome(const std::string& ref_fasta);
  ~ReferenceGenome();

  // pull out positions from genome
  const std::string GetPositions(const std::string& chrom,
                                 const int& start,
                                 const int& end);
};

#endif  // SRC_REFERENCEGENOME_H__
