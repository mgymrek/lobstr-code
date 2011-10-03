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
#ifndef __I_SATELLITE_H__
#define __I_SATELLITE_H__

#include "MSReadRecord.h"

class ISatellite
{
public:
	virtual ~ISatellite() { } ;

	/* returns TRUE if the read should be written to output.
	   FALSE if the read can be discarded */
	virtual bool ProcessRead ( MSReadRecord* read ) = 0 ;
};

#endif
