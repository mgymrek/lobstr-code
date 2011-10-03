/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>

#include "Anonymizer.h"
#include "common.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "runtime_parameters.h"

using namespace std;

Anonymizer::Anonymizer() {
  // seeding the generators
  srand((unsigned int)time(NULL));
  eng.seed((unsigned int)time(NULL));
}

Anonymizer::~Anonymizer() {}

void Anonymizer::CreateDatabase(const string& database_file) {
  TextFileReader database_reader(database_file);
  string line;
  while (database_reader.GetNextLine(&line)) {
    string STR_name;
    int allele_num;
    int start_coord;
    int copy_number;
    string allele_nucleotides;
    if (ParseAlleleDatabaseLine(line, &STR_name, &allele_num,
				&start_coord, &copy_number,
				&allele_nucleotides)) {
      if (str_id_to_allele_nucs_.find(STR_name) !=
	  str_id_to_allele_nucs_.end()) {
	str_id_to_allele_nucs_.at(STR_name).alleles.
	  push_back(allele_nucleotides);
	str_id_to_allele_nucs_.at(STR_name).copy_numbers.
	  push_back(copy_number);
      } else{
	// populate AlleleRecord
	AlleleSet allele_set;
	allele_set.start_coord = start_coord;
	allele_set.alleles.push_back(allele_nucleotides);
	allele_set.copy_numbers.push_back(copy_number);
	str_id_to_allele_nucs_.insert(pair<string, AlleleSet>
				     (STR_name, allele_set));
      }
    }
  }

  // set which allele will be used for each STR and write haplotype egg
  TextFileWriter hap_egg_writer(output_prefix + ".haplotype.tab");
  for (map<string, AlleleSet>::const_iterator it = str_id_to_allele_nucs_.begin();
       it != str_id_to_allele_nucs_.end(); ++it) {
    const AlleleSet& allele_set = it->second;
    str_id_to_allele_nucs_.at(it->first).allele_used = rand()%allele_set.alleles.size();
    // write haplotype egg
    stringstream ss;
    ss << it->first
       << "\t0\t" 
       << allele_set.start_coord << "\t"
       << allele_set.copy_numbers[allele_set.allele_used] << "\t"
       << allele_set.alleles[allele_set.allele_used];
    string result = ss.str();
    hap_egg_writer.Write(result);
  }
}

void Anonymizer::AnonymizeReads(const string& input_file,
				const string& output_file,
				const string& ids_file) {
  // get file reader depending on input type
  IFileReader* pReader = create_file_reader(input_file);
  TextFileWriter idWriter(ids_file);
  TextFileWriter pWriter(output_file);
  MSReadRecord ms_record;
  while(pReader->GetNextRecord(&ms_record)) {
    // for each read, if in list of IDS
    if (read_id_to_alignment_.find(ms_record.ID) !=
	read_id_to_alignment_.end()) {
      // call anonymize single read
      AnonymizeSingleRead(&ms_record);
      // and output ID/anonymized read
      idWriter.Write(ms_record.ID);
    }
    // output read
    if (input_type == INPUT_FASTQ) {
      pWriter.Write("@" + ms_record.ID);
    } else {
      pWriter.Write(">" + ms_record.ID);
    }
    pWriter.Write(ms_record.nucleotides);
    if (input_type == INPUT_FASTQ) { 
      pWriter.Write("+" + ms_record.ID);
      pWriter.Write(ms_record.quality_scores);
    }
  }
  delete pReader;
}

void Anonymizer::CreateReadIDToRecordMap(const std::string& str_alignment_file) {
  // open file for reading
  TextFileReader alignment_file_reader(str_alignment_file);
  string line;
  while(alignment_file_reader.GetNextLine(&line)) {
    // parse each line to get ID, start location, STR
    string STR_name;
    int start_coord;
    string read_id;
    bool reverse_complement;
    if (ParseAlignmentLine(line, &STR_name, &start_coord,
			   &read_id, &reverse_complement)) {
      AlignedRead aligned_read;
      aligned_read.read_identifier = read_id;
      aligned_read.STR_name = STR_name;
      aligned_read.start_coord = start_coord;
      aligned_read.is_reverse_complement = reverse_complement;
      if (anonymizer_debug) {
	cout << "adding " << read_id << " to dictionary" << endl;
      }
      read_id_to_alignment_.insert(pair<string, AlignedRead>
				   (read_id, aligned_read));
    }
  }
}

void Anonymizer::AnonymizeSingleRead(MSReadRecord* record_to_modify) {
  size_t length = record_to_modify->nucleotides.length();
  if (mask_n) {
    record_to_modify->nucleotides = string(length, 'N');
    return;
  } else {
    try {
      AlignedRead aligned_read = read_id_to_alignment_.at(record_to_modify->ID);
      AlleleSet allele_set = str_id_to_allele_nucs_.at(aligned_read.STR_name);
      string new_allele_nucs = allele_set.alleles.at(allele_set.allele_used).
	substr(aligned_read.start_coord-allele_set.start_coord, length);
      if (aligned_read.is_reverse_complement) {
	record_to_modify->nucleotides = reverseComplement(new_allele_nucs);
      } else{
	record_to_modify->nucleotides = new_allele_nucs;
      }
      // introduce errors
      if (anonymizer_debug) {
	cerr << "create errors in this read: " << record_to_modify->nucleotides << endl;
      }
      string nucleotides_with_errors =
	ReadWithRandomErrors(*record_to_modify);
      record_to_modify->nucleotides = nucleotides_with_errors;
      if (anonymizer_debug) {
	cerr << "after adding error: " << nucleotides_with_errors << endl << endl;
      }
    } catch (exception& e){
      cerr << "Read not in dictionary...\n";
    }
  }
}

bool Anonymizer::ParseAlleleDatabaseLine(const string& line, string* STR_name,
					 int* allele_num, int* start_coord,
					 int* copy_number,
					 string* allele_nucleotides) {
  // split string at tab
  vector<string> items;
  boost::split(items, line, boost::is_any_of("\t"));
  if (items.size() != 5) return false;
  *STR_name = items[0];
  *allele_num = atoi(items[1].c_str());
  *start_coord = atoi(items[2].c_str());
  *copy_number = atoi(items[3].c_str());
  *allele_nucleotides = items[4];
  return true;
}

bool Anonymizer::ParseAlignmentLine(const std::string& line, std::string* STR_name,
				    int* start_coord, std::string* read_id,
				    bool* reverse_complement) {
  vector<string> items;
  boost::split(items, line, boost::is_any_of("\t"));
  if (items[0] == "ID") return false;
  if (items.size() != 22) return false;
  *read_id = items[0];
  *reverse_complement = false;
  if (atoi(items[19].c_str()) == 1) {
    *reverse_complement = true;
  }
  *start_coord = atoi(items[14].c_str());
  *STR_name = items[21];
  return true;
}

string Anonymizer::ReadWithRandomErrors(const MSReadRecord& read_record) {
  std::tr1::binomial_distribution<int, double>
    error_generator(read_record.nucleotides.length(),
		    error_rate);
  // get the number of errors to introduce in the read
  int num_errors = error_generator(eng);
  
  // read nucleotides
  string read = read_record.nucleotides;
  if (num_errors == 0) {
    return read;
  }
  // generate errors
  for (int i = 0; i < num_errors; ++i) {
    // randomly choose index to change
    // TODO might change the same index twice, it's ok
    // this will happen VERY rarely, and we plan to
    // change how this is done in a future version anyway
    int index_to_corrupt = rand()%read_record.nucleotides.length();
    // randomly pick which nucleotide to switch to
    // don't switch to what's already there
    int nuc_there_already = NucToNumber(read.at(index_to_corrupt));
    int new_nuc_int = nuc_there_already;
    while (new_nuc_int == nuc_there_already) {
      // pick a new one randomly
      new_nuc_int = rand()%4;
    }
    // convert new nuc to a character
    char new_nuc_character = NumberToNuc(new_nuc_int);
    if (anonymizer_debug) {
      cerr << "changing " << read.at(index_to_corrupt)
	   << " to " << new_nuc_character
	   << " at index " << index_to_corrupt << endl;
    }
    // set that character
    read.replace(index_to_corrupt, 1, 1, new_nuc_character);
  }
  return read;
}

int Anonymizer::NucToNumber(char nuc) {
  switch(nuc) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g': 
    return 2;
  case 'T':
  case 't':
    return 3;
  case 'N':
  case 'n':
    return 4;
  default:
    // NOTE THIS BETTER NOT HAPPEN
    return -1;
  }
}

char Anonymizer::NumberToNuc(int nuc) {
  switch(nuc) {
  case 0:
    return 'A';
  case 1:
    return 'C';
  case 2:
    return 'G';
  case 3:
    return 'T';
  default:
    return 'N';
  }
}
