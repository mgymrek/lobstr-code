#include <stdio.h>
#include <string.h>
#include "main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0.2"
#endif

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: lobSTRIndex (build Burrows-Wheeler transformation using BWA)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Melissa Gymrek <mgymrek@mit.edu>\n\n");
	fprintf(stderr, "Usage:   lobSTRIndex index [OPTIONS]\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	else { return bwa_index(argc-1, argv+1);}
	return 0;
}
