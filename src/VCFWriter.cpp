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

static const std::string GetSTRVar(int allele) {
  stringstream strvar;
  strvar << "<STRVAR:" << allele << ">";
  return strvar.str();
}

VCFWriter::VCFWriter(const string& filename)
  : TextFileWriter(filename) {
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
  output_stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count in genotypes, for each ALT allele, in the same order as listed\">" << endl;
  output_stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl;
  output_stream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">" << endl;
  output_stream << "##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Repeat motif\">" << endl;
  output_stream << "##INFO=<ID=REF,Number=1,Type=Float,Description=\"Reference copy number\">" << endl;
  output_stream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type\">" << endl;
  // ALT fields
  output_stream << "##ALT=<ID=STRVAR,Description=\"Short tandem variation\">" << endl;
  // FORMAT fields
  output_stream << "##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description=\"All reads aligned to locus\">" << endl;
  output_stream << "##FORMAT=<ID=ALLPARTIALREADS,Number=1,Type=String,Description=\"All partially spanning reads aligned to locus\">" << endl;  
  output_stream << "##FORMAT=<ID=CONFLICT,Number=1,Type=Integer,Description=\"Number of reads conflicting genotype call\">" << endl;
  output_stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
  output_stream << "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">" << endl;
  output_stream << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods (log10 scaled)\">" << endl;
  output_stream << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << endl;
  output_stream << "##FORMAT=<ID=GPP,Number=G,Type=Float,Description=\"Genotype Posterior probabilities (phred scaled, -10log10)\">" << endl;
  output_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  output_stream << "##FORMAT=<ID=PP,Number=1,Type=Float,Description=\"Posterior probability of call\">" << endl;
  output_stream << "##FORMAT=<ID=MP,Number=1,Type=Float,Description=\"Upper bound on maximum partially spanning allele\">" << endl;
  output_stream << "##FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Coverage by partially spanning reads\">" << endl;
  output_stream << "##FORMAT=<ID=S1,Number=1,Type=Float,Description=\"Allele 1 marginal posterior\">" << endl;
  output_stream << "##FORMAT=<ID=S2,Number=1,Type=Float,Description=\"Allele 2 marginal posterior\">" << endl;
  output_stream << "##FORMAT=<ID=STITCH,Number=1,Type=Integer,Description=\"Number of stitched reads\">"<< endl;
  output_stream << "##FORMAT=<ID=SUPP,Number=1,Type=Integer,Description=\"Number of reads supporting genotype call\">" << endl;
  // header columns
  output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample << endl;
  LoadPositionsToExclude();
}

void VCFWriter::LoadPositionsToExclude() {
  TextFileReader expos(exclude_positions_file);
  string line;
  while (expos.GetNextLine(&line)) {
    vector<string> items;
    split(line, '\t', items);
    if (items.size() == 0) break;
    if (items.size() != 2) {
      errx(1, "Error, exclude-pos file has invalid format.");
    }
    string chrom;
    int pos;
    chrom = items[0];
    pos = atoi(items[1].c_str());
    if (pos_to_exclude.find(chrom) == pos_to_exclude.end()) {
      set<int> pos_list;
      pos_to_exclude[chrom] = pos_list;
    }
    pos_to_exclude[chrom].insert(pos);
  }
}

VCFWriter::~VCFWriter() {}

void VCFWriter::WriteRecord(const STRRecord& str_record) {
  if (pos_to_exclude.find(str_record.chrom) !=
      pos_to_exclude.end()) {
    if (pos_to_exclude[str_record.chrom].count(str_record.start) != 0) {
      return;
    }  
  }
  if ((str_record.allele1 == -10000 || str_record.allele2 == -10000) &&
      str_record.partial_coverage == 0) {
    return;
  }
  if (!(str_record.stop > 0 && str_record.start > 8 &&
        str_record.stop > str_record.start)) {
    return;
  }
  // CHROM
  output_stream << str_record.chrom << "\t";
  // POS
  output_stream << str_record.start+1 << "\t";
  // ID
  output_stream << ".\t";
  // REF
  output_stream << str_record.ref_allele << "\t";
  // ALT
  stringstream ac;
  stringstream genotype_string;
  int alt_allele_count = 1;
  int allele1_num = 0;
  int allele2_num = 0;
  if (str_record.alleles_to_include.size() == 0) {
    output_stream << ".\t";
    ac << ".";
  } else {
    for (vector<int>::const_iterator it = str_record.alleles_to_include.begin();
         it != str_record.alleles_to_include.end(); it++) {
      if (*it != 0) {
        output_stream << GetSTRVar(*it);
        int count = 0;
        if (str_record.allele1 == *it) {
          count++;
          allele1_num = alt_allele_count;
        }
        if (str_record.allele2 == *it) {
          count++;
          allele2_num = alt_allele_count;
        }
        ac << count;
        if (*it != str_record.alleles_to_include.at(str_record.alleles_to_include.size()-1)) {
          ac << ",";
          output_stream << ",";
        } 
        alt_allele_count++;
      }
    }
    output_stream << "\t";
  }
  if (str_record.coverage == 0) {
    genotype_string << "./.";
  } else {
    genotype_string << allele1_num << "/" << allele2_num;
  }
  // QUAL
  float qual = -10*log10(1-str_record.score); // score is posterior prob of call
  if (qual > POSITIVE_INFINITY) {qual = POSITIVE_INFINITY;}
  if (str_record.coverage == 0) {qual = NOQUAL;}
  output_stream << qual << "\t";
  // FILTER
  output_stream << ".\t";
  // INFO
  if (ac.str() != ".") output_stream << "AC=" << ac.str() << ";";
  output_stream  << "AN=" << (str_record.coverage==0?0:2) << ";"
                 << "END=" << str_record.stop+1 << ";"
                 << "MOTIF=" << str_record.repseq << ";"
                 << "REF=" << str_record.refcopy << ";"
                 << "VT=STR" << "\t";
  // FORMAT
  output_stream << "GT:ALLREADS:ALLPARTIALREADS:CONFLICT:DP:GB:GL:PL:GPP:MP:PC:PP:S1:S2:STITCH:SUPP" << "\t";
  // Sample info
  size_t num_alleles = str_record.alleles_to_include.size()+1;
  size_t num_allele_pairs = (num_alleles*(num_alleles+1))/2;
  vector<float> genotype_likelihoods; // log10(lik)
  genotype_likelihoods.resize((size_t)(num_allele_pairs));
  vector<float> genotype_posteriors; // -10*log10(posteriors)
  genotype_posteriors.resize((size_t)(num_allele_pairs));
  for (size_t i = 0; i < num_alleles; i++) {
    for (size_t j = 0; j <= i; j++) {
      // from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41 F(j/k) = (k*(k+1)/2)+j.
      size_t index = i*(i+1)/2+j;
      int a1,a2;
      if (i == 0) {a1 = 0;}
      else {a1 = str_record.alleles_to_include.at(i-1);}
      if (j == 0) {a2 = 0;}
      else {a2 = str_record.alleles_to_include.at(j-1);}
      pair<int, int> atype;
      if (a1 < a2) atype = pair<int,int>(a1,a2);
      else atype = pair<int,int>(a2,a1);
      if (str_record.likelihood_grid.find(atype) != str_record.likelihood_grid.end()&&
          str_record.posterior_grid.find(atype) != str_record.posterior_grid.end()) {
        float lik = (str_record.likelihood_grid.at(atype));
        if (lik < NEGATIVE_INFINITY) {lik = NEGATIVE_INFINITY;}
        float post = -10*log10(str_record.posterior_grid.at(atype));
        if (post > POSITIVE_INFINITY) {post = POSITIVE_INFINITY;}
        genotype_likelihoods[index] = lik;
        genotype_posteriors[index] = post;
      }
    }
  }
  stringstream genotype_likelihoods_string;
  stringstream genotype_scaled_likelihoods_string;
  stringstream genotype_posteriors_string;
  if (str_record.coverage == 0) {
    genotype_likelihoods_string << ".";
    genotype_scaled_likelihoods_string << ".";
    genotype_posteriors_string << ".";
  } else {
    for (size_t i = 0; i < genotype_likelihoods.size(); i++) {
      genotype_likelihoods_string << genotype_likelihoods.at(i);
      genotype_scaled_likelihoods_string << static_cast<int>(-10*(genotype_likelihoods.at(i)-str_record.max_log_lik));
      genotype_posteriors_string << genotype_posteriors.at(i);
      if (i != genotype_likelihoods.size()-1) {
        genotype_likelihoods_string << ",";
        genotype_scaled_likelihoods_string << ",";
        genotype_posteriors_string << ",";
      }
    }  
  }
  stringstream gbstring;
  stringstream score1string;
  stringstream score2string;
  stringstream max_partial_string;
  if (str_record.allele1 == MISSING || str_record.coverage == 0) {
    gbstring << ".";
    score1string << ".";
  } else {
    gbstring << str_record.allele1;
    score1string << str_record.allele1_score;
  }
  gbstring << "/";
  if (str_record.allele2 == MISSING || str_record.coverage == 0) {
    gbstring << ".";
    score2string << "."; 
  } else {
    gbstring << str_record.allele2;
    score2string << str_record.allele2_score;
  }
  if (str_record.partial_coverage == 0) {
    max_partial_string << "."; 
  } else {
    max_partial_string << str_record.max_partial_string;
  }
  output_stream << genotype_string.str() << ":"
                << str_record.readstring << ":"
                << str_record.partialreadstring << ":"
                << str_record.conflicting << ":"
                << str_record.coverage << ":"
                << gbstring.str() << ":"
                << genotype_likelihoods_string.str() << ":"
                << genotype_scaled_likelihoods_string.str() << ":"
                << genotype_posteriors_string.str() << ":"
                << max_partial_string.str() << ":"
                << str_record.partial_coverage << ":"
                << str_record.score << ":"
                << score1string.str() << ":"
                << score2string.str() << ":"
                << str_record.num_stitched << ":"
                << str_record.agreeing << endl;
}
