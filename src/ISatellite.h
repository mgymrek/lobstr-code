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

#ifndef SRC_ISATELLITE_H__
#define SRC_ISATELLITE_H__

#include "src/MSReadRecord.h"

/*
  Note: this class is provided for the case
  where we may want to include other detection
  methods. Right now this class is only
  implemented by STRDetector.
 */

class ISatellite {
 public:
  virtual ~ISatellite() { }
	/* returns TRUE if the read should be written to output.
	   FALSE if the read can be discarded */
  virtual bool ProcessRead(MSReadRecord* read) = 0;
};

#endif  // SRC_ISATELLITE_H__
