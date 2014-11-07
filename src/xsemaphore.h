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
#ifndef __XSEMAPHORE_H__
#define __XSEMAPHORE_H__

/* A simple semaphore wrapper class,
   hiding the differences between Mac-OS-X and GNU/Linux semaphores.

   Mac-OS-X does not implement that 'sem_init()' and related functions
   (functions return errno=ENOSYS).

   Instead, use Mac-OS's dispatch semaphores.

   Detection of Mac-OS is performed in './configure' (see 'configure.ac').
   Value of 'USE_MACOS_DISPTACH' is stored in 'config.h' by './configure'.
*/

#include "config.h"

#ifdef USE_MACOS_DISPTACH
#include <dispatch/dispatch.h>

class XSemaphore
{
	dispatch_semaphore_t s;
public:
	XSemaphore(int count) :
		s(NULL)
	{
		s = dispatch_semaphore_create(count);
		if (s==NULL)
			err(1,"%s:%d: dispatch_semaphore_create() failed",__FILE__,__LINE__); \
	}

	~XSemaphore()
	{
		if (s==NULL)
			errx(1,"%s:%d: s==NULL (before disp-release)",__FILE__,__LINE__); \
		dispatch_release(s);
		s = NULL;
	}

	void post()
	{
		if (s==NULL)
			errx(1,"%s:%d: s==NULL (before disp-sem-sig)",__FILE__,__LINE__); \
		dispatch_semaphore_signal(s);
	}

	void wait()
	{
		if (s==NULL)
			errx(1,"%s:%d: s==NULL (before disp-sem-wait)",__FILE__,__LINE__); \
		dispatch_semaphore_wait(s,DISPATCH_TIME_FOREVER);
	}
};

#else

#include <semaphore.h>

class XSemaphore
{
	sem_t s;

public:
	XSemaphore(int count)
	{
		if (sem_init(&s, 0, count)!=0)
			err(1,"%s:%d: sem_init failed",__FILE__,__LINE__); \
	}


	~XSemaphore()
	{
		if (sem_destroy(&s) !=0)
			err(1,"%s:%d: sem_destroy failed",__FILE__,__LINE__); \
	}

	void post()
	{
		if (sem_post(&s) != 0)
			err(1,"%s:%d: sem_post failed",__FILE__,__LINE__); \
	}

	void wait()
	{
		if (sem_wait(&s) != 0)
			err(1,"%s:%d: sem_wait failed",__FILE__,__LINE__); \
	}
};

#endif

#endif
