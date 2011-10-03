/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

 This file is part of MicroSatelliteDetector.

 MicroSatelliteDetector is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MicroSatelliteDetector is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MicroSatelliteDetector.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __DEBUG_SATELLITE_H__
#define __DEBUG_SATELLITE_H__

#include "ReadRecord.h"

class DebugSatellite : public ISatellite
{
public:
	DebugSatellite();

	bool ProcessRead ( ReadRecord& read ) ;
};

#endif
