/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "GST.h"
#include "MSRecord.h"

#ifndef GST_CONTAINER_H_
#define GST_CONTAINER_H_

class GSTContainer {
 public:
  GSTContainer();
  ~GSTContainer();
  // clear the GSTS
  void clear();

  // get the right or left GST for a given STR
  bool getGST(const std::string& str, GST* gst, bool left);
  
  // get the MS locus record for an STR locus
  bool getMSLociRecord(const int& key, MSRecord* ms_loci_record);

  // Build GSTS from ms table file
  void buildGSTs(const std::string& ms_table_filename,
		 const std::string& genome_filename);

  // map of STR to left GST
  map<std::string, GST*> gsts_left_flanks;

  // map of STR to right GST
  map<std::string, GST*> gsts_right_flanks;

  // map of STR id # to its record
  map<int, MSRecord> ms_loci_dict;
};

#endif /* GST_CONTAINER_H_ */
