/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef ANONYMIZER_H_
#define ANONYMIZER_H_

#include "MSReadRecord.h"

#include <boost/tr1/random.hpp>
#include <map>
#include <string>
#include <vector>

// information about a set of STR alleles

struct AlleleSet {
  std::vector<std::string> alleles;
  std::vector<int> copy_numbers;
  int start_coord;
  int allele_used;
};

struct AlignedRead {
  std::string read_identifier;
  std::string STR_name;
  int start_coord;
  bool is_reverse_complement;
};

class Anonymizer {
 public:
  Anonymizer();
  ~Anonymizer();

  // create database of marker -> alleles
  void CreateDatabase(const std::string& database_file);

  // anonymize reads from one file and write them to another
  // writes IDs of modified reads to a file
  void AnonymizeReads(const std::string& input_file,
		      const std::string& output_file,
		      const std::string& ids_file);

  // create map of identifier -> read record
  void CreateReadIDToRecordMap(const std::string& str_alignment_file);

 protected:
  // Randomize a single read
  void AnonymizeSingleRead(MSReadRecord* record_to_modify, std::string* name);
  
  // parse a line of the database file
  static bool ParseAlleleDatabaseLine(const std::string& line, std::string* STR_name,
				      int* allele_num, int* start_coord,
				      int* copy_number,
				      std::string* allele_nucleotides);

  // parse a line of the alignment file
  static bool ParseAlignmentLine(const std::string& line, std::string* STR_name,
				 int* start_coord, std::string* read_id,
				 bool* reverse_complement);

  // map of STR identifier to alleles and start coordinate
  std::map<std::string, AlleleSet> str_id_to_allele_nucs_;

  // map of ID to alignment
  std::map<std::string, AlignedRead> read_id_to_alignment_;

  // introduce random errors
  std::string ReadWithRandomErrors(const MSReadRecord& read_record);

  // helper to convert nucleotide to number
  static int NucToNumber(char nuc);

  // helper to convert number to nucleotide
  static char NumberToNuc(int nuc);

  // random number generator
  std::tr1::mt19937 eng;

  friend class AnonymizerTest;
};

#endif /* GENOTYPER_H_ */
