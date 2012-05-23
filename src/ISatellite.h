/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
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
