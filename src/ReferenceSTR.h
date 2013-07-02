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

#ifndef SRC_REFERENCESTR_H_
#define SRC_REFERENCESTR_H_

using namespace std;

/*
  Struct to keep track of info for a single STR locus
 */
struct ReferenceSTR {
  // Locus properties
  std::string chrom;
  int start;
  int stop;
};

#endif  // SRC_REFERENCESTR_H_
