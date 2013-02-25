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

#ifndef SRC_RUNINFO_H_
#define SRC_RUNINFO_H_

#include <sstream>
#include <string>
#include <vector>

/*
  RunInfo stores information about the run that will be uploaded to Amazon for further analysis
 */

class RunInfo {
 public:
  // time at which run started
  std::string starttime;
  // time at which run ended
  std::string endtime;
  // git version
  std::string gitversion;
  // machine type
  std::string machtype;
  // run time params
  std::string params;
  // error msg
  std::string error;
  // print to string
  std::string PrintToString() {
    std::stringstream ss;
    ss << "Starttime: " << starttime << std::endl;
    ss << "Endtime: " << endtime << std::endl;
    ss << "Git version: " << gitversion << std::endl;
    ss << "Machtype: " << machtype << std::endl;
    ss << "Params: " << params << std::endl;
    if (!error.empty()) {
      ss << "Error: " << error << std::endl;
    }
    return ss.str();
  }
};

#endif  // SRC_RUNINFO_H_
