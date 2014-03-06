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

#include <err.h>
#include <errno.h>
#include <math.h>

#include <algorithm>
#include <string>

#include "src/common.h"
#include "src/runtime_parameters.h"
#include "src/TextFileReader.h"
#include "src/VCFWriter.h"

using namespace std;

const float NEGATIVE_INFINITY = -1000;
const float POSITIVE_INFINITY = 1000;
const float NOQUAL = -1;

static const std::string GetReference() {
  vector<string> items;
  split(user_defined_arguments, ';', items);
  for (size_t i = 0; i < items.size(); i++) {
    vector<string> paramitems;
    split(items[i], '=', paramitems);
    if (paramitems[0] == "index-prefix") {
      return paramitems[1]+"ref.fa";
    }
  }
  return "";
}

VCFWriter::VCFWriter(const string& filename, const vector<string>& samples)
  : TextFileWriter(filename) {
  // update run info
  run_info.samples = samples;
  run_info.num_calls.resize(samples.size());
  run_info.num_calls5x.resize(samples.size());
  run_info.total_coverage.resize(samples.size());
  run_info.total_agree.resize(samples.size());
  run_info.calltype_by_period.resize(6);
  for (size_t i = 0; i < 6; i++) {
    run_info.calltype_by_period.at(i).resize(4);
  }
  // VCF header
  string lobstr_header = user_defined_arguments;
  lobstr_header = string_replace(lobstr_header, "@CO\t# version=", "");
  lobstr_header = string_replace(lobstr_header, "\n", "");
  string allelotype_header = user_defined_arguments_allelotyper;
  allelotype_header = string_replace(allelotype_header, "# version=", "");
  output_stream << "##fileformat=VCFv4.1" << endl;
  output_stream << "##fileDate=" << currentDateTime() << endl;
  output_stream << "##source=" << lobstr_header << endl;
  output_stream << "##source=" << allelotype_header << endl;
  if (!GetReference().empty()) {
    output_stream << "##reference=file://" << GetReference() << endl;
  }
  // INFO fields
  output_stream << "##INFO=<ID=RPA,Number=A,Type=Float,Description=\"Repeats per allele\">" << endl;
  output_stream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">" << endl;
  output_stream << "##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Canonical repeat motif\">" << endl;
  output_stream << "##INFO=<ID=REF,Number=1,Type=Float,Description=\"Reference copy number\">" << endl;
  output_stream << "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference STR track length in bp\">" << endl;
  output_stream << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat motif\">" << endl;
  output_stream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type\">" << endl;
  // FORMAT fields
  output_stream << "##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description=\"All reads aligned to locus\">" << endl;
  output_stream << "##FORMAT=<ID=AML,Number=1,Type=String,Description=\"Allele marginal likelihood ratio scores\">" << endl;
  output_stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
  output_stream << "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">" << endl;
  output_stream << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << endl;
  output_stream << "##FORMAT=<ID=Q,Number=1,Type=Float,Description=\"Likelihood ratio score of allelotype call\">" << endl;
  output_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  output_stream << "##FORMAT=<ID=STITCH,Number=1,Type=Integer,Description=\"Number of stitched reads\">"<< endl;

  if (include_gl)
    output_stream << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods for genotypes as defined in the VCF specification\">" << endl;

  // header columns
  output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (vector<string>::const_iterator it = samples.begin();
       it != samples.end(); it++) {
    output_stream << "\t" << (*it);
  }
  output_stream << endl;
}

std::string VCFWriter::GetSTRVar(const string& refseq, const string& ref_repseq, int allele) {
  stringstream strvar;
  if (allele == 0) {
    return refseq;
  } else if (allele < 0) {
    return refseq.substr(0,refseq.size()+allele);
  } else {
    strvar << refseq;
    int offset = refseq.size() % ref_repseq.size();
    stringstream toadd;
    for (int i = 0; i < allele/static_cast<int>(ref_repseq.length()) + 2; i++) {
      toadd << ref_repseq;
    }
    strvar << toadd.str().substr(offset, allele);
    return strvar.str();
  }
}

VCFWriter::~VCFWriter() {}

void VCFWriter::WriteRecord(const STRRecord& str_record) {
  // CHROM
  output_stream << str_record.chrom << "\t";
  // POS
  output_stream << str_record.start << "\t";
  // ID
  if (str_record.name.empty()) {
    output_stream << ".\t";
  } else {
    output_stream << str_record.name << "\t";
  }
  // REF
  output_stream << str_record.ref_allele << "\t";
  // ALT
  stringstream repeats_per_allele;
  map<int,int> allele_to_index;
  allele_to_index[0] = 0;
  int ind = 1;
  if (str_record.alleles_to_include.size() == 0) return;
  if (str_record.alleles_to_include.size() == 1) { // if only 1 allele, will always be reference
    output_stream << ".\t";
    repeats_per_allele << ".";
  } else {
    for (size_t i = 1; i < str_record.alleles_to_include.size(); i++) {
      const int allele = str_record.alleles_to_include.at(i);
      allele_to_index[allele] = ind;
      ind++;
      output_stream << GetSTRVar(str_record.ref_allele, str_record.repseq_in_ref, allele);
      repeats_per_allele << (static_cast<float>(allele)/static_cast<float>(str_record.repseq.size()) + str_record.refcopy);
      if (i != str_record.alleles_to_include.size() - 1) {
	output_stream << ",";
	repeats_per_allele << ",";
      } 
    }
    output_stream << "\t";
  }

  // QUAL
  float qual = 0;
  for (size_t i = 0; i < str_record.prob_ref.size(); i++) {
    qual += -10*(str_record.prob_ref.at(i));
  }
  if (qual < 0) qual = 0;
  output_stream << qual << "\t";
  // FILTER
  output_stream << ".\t";
  // INFO
  output_stream << "END=" << str_record.stop << ";"
		<< "MOTIF=" << str_record.repseq << ";"
		<< "REF=" << str_record.refcopy << ";"
		<< "RL=" << str_record.stop - str_record.start<< ";"
		<< "RPA=" << repeats_per_allele.str() << ";"
		<< "RU=" << str_record.repseq_in_ref << ";"
		<< "VT=STR" << "\t";
  // FORMAT
  if (include_gl)
    output_stream << "GT:ALLREADS:AML:DP:GB:GL:PL:Q:STITCH";
  else
    output_stream << "GT:ALLREADS:AML:DP:GB:PL:Q:STITCH";

  // Sample info  
  for (size_t i = 0; i < str_record.samples.size(); i++) {
    output_stream << "\t";
    WriteSample(str_record, i, allele_to_index);
  }
  output_stream << endl;
}

void VCFWriter::WriteSample(const STRRecord& str_record, size_t sample_index,
			    map<int,int> allele_to_index) {
  stringstream genotype_string;
  if (str_record.coverage.at(sample_index) == 0) {
    output_stream << "./.";
    return;
  } else {
    run_info.num_calls.at(sample_index)++;
    if (str_record.coverage.at(sample_index) >= 5) {
      run_info.num_calls5x.at(sample_index)++;
    }
    run_info.total_coverage.at(sample_index) += str_record.coverage.at(sample_index);
    run_info.total_agree.at(sample_index) += static_cast<float>(str_record.agreeing.at(sample_index))/
      static_cast<float>(str_record.coverage.at(sample_index));
    size_t calltype_index = -1;
    const int& allele1 = allele_to_index[str_record.allele1.at(sample_index)];
    const int& allele2 = allele_to_index[str_record.allele2.at(sample_index)];
    if (allele1 == allele2) {
      if (allele1 == 0) {
	calltype_index = 0;
      } else {
	calltype_index = 2;
      }
    } else {
      if (allele1 == 0 || allele2 == 0) {
	calltype_index = 1;
      } else {
	calltype_index = 3;
      }
    }
    run_info.calltype_by_period.at(str_record.period-1).at(calltype_index)++;
    genotype_string << allele1 << "/" << allele2;
  }

  size_t num_alleles = str_record.alleles_to_include.size();
  size_t num_allele_pairs = (num_alleles*(num_alleles+1))/2;
  vector<float> genotype_likelihoods; // log10(lik)
  genotype_likelihoods.resize((size_t)(num_allele_pairs));

  // go over possible genotypes. assume 0 at front of alleles to include
  for (size_t i = 0; i < num_alleles; i++) {
    for (size_t j = 0; j <= i; j++) {
      size_t index = i*(i+1)/2+j;
      int a1 = str_record.alleles_to_include.at(i);
      int a2 = str_record.alleles_to_include.at(j);
      pair<int, int> atype;
      if (a1 < a2) atype = pair<int,int>(a1,a2);
      else atype = pair<int,int>(a2,a1);
      if (str_record.likelihood_grid.at(sample_index).find(atype)
	  != str_record.likelihood_grid.at(sample_index).end()) {
        float lik = (str_record.likelihood_grid.at(sample_index).at(atype));
        if (lik < NEGATIVE_INFINITY) {lik = NEGATIVE_INFINITY;}
        genotype_likelihoods[index] = lik;
      }
    }
  }

  stringstream genotype_likelihoods_string;
  if (include_gl){
    genotype_likelihoods_string.precision(2);
    if (str_record.coverage.at(sample_index) == 0)
      genotype_likelihoods_string << ".";
    else {
      for (size_t i = 0; i < genotype_likelihoods.size(); i++) {
	genotype_likelihoods_string << fixed << genotype_likelihoods.at(i);
	if (i != genotype_likelihoods.size()-1)
	  genotype_likelihoods_string << ",";
      }  
    }
  }

  stringstream genotype_scaled_likelihoods_string;
  if (str_record.coverage.at(sample_index) == 0) {
    genotype_scaled_likelihoods_string << ".";
  } else {
    for (size_t i = 0; i < genotype_likelihoods.size(); i++) {
      genotype_scaled_likelihoods_string << static_cast<int>
        (-10*(genotype_likelihoods.at(i)-str_record.max_log_lik.at(sample_index)));
      if (i != genotype_likelihoods.size()-1) {
        genotype_scaled_likelihoods_string << ",";
      }
    }  
  }

  stringstream gbstring;
  stringstream marginal_lik_score_string;
  if (str_record.coverage.at(sample_index) == 0) {
    gbstring << "./.";
    marginal_lik_score_string << "./.";
  } else {
    gbstring << str_record.allele1.at(sample_index) << "/"
             << str_record.allele2.at(sample_index);
    marginal_lik_score_string << str_record.allele1_marginal_lik_score.at(sample_index) << "/"
                              << str_record.allele2_marginal_lik_score.at(sample_index);
  }

  // Output format fields
  output_stream << genotype_string.str() << ":"
                << str_record.readstring.at(sample_index) << ":"
                << marginal_lik_score_string.str() << ":"
                << str_record.coverage.at(sample_index) << ":"
                << gbstring.str() << ":";

  if (include_gl)
    output_stream << genotype_likelihoods_string.str() << ":";
  
  output_stream << genotype_scaled_likelihoods_string.str() << ":"
                << str_record.max_lik_score.at(sample_index) << ":"
                << str_record.num_stitched.at(sample_index);
}
