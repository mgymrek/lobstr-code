/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "GSTContainer.h"
#include "MSTableFileReader.h"
#include "runtime_parameters.h"

using namespace std;

GSTContainer::GSTContainer() {};
GSTContainer::~GSTContainer() {};

void GSTContainer::clear() {
  for (map<string, GST*>::iterator it = gsts_left_flanks.begin();
       it != gsts_left_flanks.end(); ++it) {
    if (it->second != NULL) delete it->second;
  }
  for (map<string, GST*>::iterator it = gsts_right_flanks.begin();
       it != gsts_right_flanks.end(); ++it) {
    if (it->second != NULL) delete it->second;
  }
  gsts_left_flanks.clear();
  gsts_right_flanks.clear();
  ms_loci_dict.clear();
};

bool GSTContainer::getGST(const string& str, GST* gst, bool left) {
  if (left) {
    if (gsts_left_flanks.find(str) != gsts_left_flanks.end()) {
      gst = gsts_left_flanks.at(str);
      return true;
    }
  } else {
    if (gsts_right_flanks.find(str) != gsts_right_flanks.end()) {
      gst = gsts_right_flanks.at(str);
      return true;
    }
  }
  return false;
}

bool GSTContainer::getMSLociRecord(const int& key,
				   MSRecord* ms_loci_record) {
  if (ms_loci_dict.find(key) != ms_loci_dict.end()) {
    *ms_loci_record = ms_loci_dict.at(key);
    return true;
  }
  return false;
}

void GSTContainer::buildGSTs(const string& ms_table_filename,
			     const string& genome_filename) {
  MSTableFileReader ms_table_reader = MSTableFileReader(ms_table_filename,
							genome_filename);
  map<string, list<MSRecord> > seq_to_records_dict;
  int ms_counter = 0;
  MSRecord ms_record;
  // make map of repeat -> msrecord for that repeat
  while(ms_table_reader.GetNextRecord(&ms_record)) {
    ms_record.seqid = ms_counter;
    // process it and populate dictionary
    if (ms_record.repeat.length() <= max_period &&
	ms_record.repeat.length() >= min_period) {
      if (seq_to_records_dict.count(ms_record.repeat) != 0) {
	seq_to_records_dict.at(ms_record.repeat).push_back(ms_record);
      } else {
	list<MSRecord> ms_record_list;
	ms_record_list.push_back(ms_record);
	seq_to_records_dict.insert(pair<string, list<MSRecord> >
				   (ms_record.repeat, ms_record_list));
      }
      ms_loci_dict.insert(pair<int, MSRecord>(ms_counter, ms_record));
    }
    ms_counter++;
  }
  // build GSTs
  for (map<string, list<MSRecord> >::const_iterator
	 it = seq_to_records_dict.begin();
       it != seq_to_records_dict.end(); it++) {
    gsts_left_flanks.insert(pair<string, GST*>(it->first, new GST(it->second, true)));
    gsts_right_flanks.insert(pair<string, GST*>(it->first, new GST(it->second, false)));
  }
}
