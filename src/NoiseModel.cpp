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

#include <cmath>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "src/common.h"
#include "src/logistic_regression.h"
#include "src/NoiseModel.h"
#include "src/runtime_parameters.h"
#include "src/TextFileReader.h"
#include "src/TextFileWriter.h"

const size_t MIN_TRAIN_COV = 5;
const float MIN_TRAIN_AGREE = 0.5;
const float SMALL_CONST = 1e-4;
const size_t MIN_READS_FOR_TRAINING = 1000;
const size_t NUM_PARAMS = 5;

using namespace std;

/* Same as R dpois */
static float dpois(int step, float mean) {
  if (step < 0) return 0;
  float p = exp(-1*mean);
  for (int i = 0; i < step; i++) {
    p = p*mean;
    p = p/(i+1);
  }
  return p;  
}

/* same as R dgeom */
static float dgeom(int step, float psucc) {
  if (step < 0) return 0;
  float p = psucc;
  for (int i = 0; i < step; i++) {
    p = p*(1-psucc);
  }
  return p;
}

/* Return step size module period */
static float mmod(int step, int period) {
  int astep = abs(step);
  int n1 = astep % period;
  int n2 = abs(astep % period - period);
  return (n1 < n2) ? n1:n2;
}

NoiseModel::NoiseModel(const string& strinfofile,
                       const vector<string>& haploid_chroms_,
		       const string& noise_model) {
  // Read STR info
  ReadSTRInfo(strinfofile);
  // Set noise model filenames
  stutter_model_filename = noise_model + ".stuttermodel";
  stepsize_model_filename = noise_model + ".stepmodel";
  // Determine MIN PERIOD
  if (command == "classify") {
    MIN_PERIOD = DetermineMinPeriod(stepsize_model_filename);
  } else {
    MIN_PERIOD = 1;
  }
  // Set haploid chromosomes
  haploid_chroms = haploid_chroms_;
  // initialize pdf
  for (int i = 0; i < MAX_PERIOD-MIN_PERIOD+1; i++) {
    vector<float> v0(0);
    vector<float> v1(0);
    for (int j = 0; j < MAX_PERIOD*3*2+1; j++) {
      v1.push_back(0);
    }
    pdf_model.push_back(v1);
  }
}

int NoiseModel::DetermineMinPeriod(const string& stepsize_model_filename) {
  TextFileReader modelFile(stepsize_model_filename);
  string line;
  int numlines = 0;
  while (modelFile.GetNextLine(&line)) {
    numlines += 1;
  }
  if (numlines == 17) {
    return 2;
  } else {
    return 1;
  }
}

void NoiseModel::ReadSTRInfo(const string& filename) {
  TextFileReader infoFile(filename);
  string line;
  while (infoFile.GetNextLine(&line)) {
    vector<string> items;
    split(line, '\t', items);
    if (items.size() == 0) break;
    if (items.size() != 6) {
      PrintMessageDieOnError("STR info file has invalid format", ERROR);
    }
    STRINFO info;
    info.chrom = items[0];
    info.start = atoi(items[1].c_str());
    info.end = atoi(items[2].c_str());
    info.score = atof(items[3].c_str());
    info.gc = atof(items[4].c_str());
    info.entropy = atof(items[5].c_str());
    pair<string, int> coord (info.chrom, info.start);
    strInfo.insert(pair< pair<string, int>, STRINFO >
                   (coord, info));
  }
}

void NoiseModel::Train(ReadContainer* read_container) {
  vector<AlignedRead> reads_for_training(0);
  map<int, map <int,int> > step_size_by_period;
  for (map<pair<string, int>, list<AlignedRead> >::iterator
         it = read_container->aligned_str_map_.begin();
       it != read_container->aligned_str_map_.end(); it++) {
    const list<AlignedRead>& aligned_reads = it->second;
    // check if haploid
    bool is_haploid = false;
    if ((find(haploid_chroms.begin(), haploid_chroms.end(),
              aligned_reads.front().chrom) != haploid_chroms.end() ) ||
        find(haploid_chroms.begin(), haploid_chroms.end(), "all")!= haploid_chroms.end()) {
      is_haploid = true;
    }
    if (!is_haploid) continue;
    // check if in our list of STRs
    pair<string, int> coord(aligned_reads.front().chrom,
                            aligned_reads.front().msStart);
    if (strInfo.find(coord) == strInfo.end()) continue;
    // check if has unique mode
    float unique_mode;
    if (HasUniqueMode(aligned_reads, &unique_mode)) {
      for (list<AlignedRead>::const_iterator it2 =
             aligned_reads.begin(); it2 != aligned_reads.end(); it2++) {
        if (it2->mapq != 0) continue;
        AlignedRead aread;
        aread.chrom = it2->chrom; aread.msStart = it2->msStart;
        aread.msEnd = it2->msEnd; aread.read_start = it2->read_start;
        aread.period = it2->period; aread.diffFromRef = it2->diffFromRef;
        aread.stutter = aread.diffFromRef != unique_mode;
        reads_for_training.push_back(aread);
        if (aread.stutter && abs(aread.diffFromRef) <= 3*aread.period) {
          step_size_by_period[aread.period][aread.diffFromRef]++;
        }
      }
    }
  }

  // *** Step 1: model of probability of stuttering *** //
  if (debug) {
    PrintMessageDieOnError("[NoiseModel] Fitting mutation prob", DEBUG);
  }
  FitMutationProb(reads_for_training);
  // *** Step 2: Fitting step size distribution *** //
  if (debug) {
    PrintMessageDieOnError("[NoiseModel] Fitting step size prob", DEBUG);
  }
  FitStepProb(step_size_by_period);
}

bool NoiseModel::HasUniqueMode(const list<AlignedRead>&
                               aligned_reads,
                               float* mode) {
  if (aligned_reads.size() < MIN_TRAIN_COV) {
    return false;
  }
  
  int top_copy_count = 0;
  float top_copy;
  int second_copy_count = 0;
  float second_copy;
  map<float, int> copy_to_count;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    float copy = (*it).diffFromRef;
    copy_to_count[copy]++;
  }
  for (map<float, int>::const_iterator it = copy_to_count.begin();
       it != copy_to_count.end(); it++) {
    int count = it->second;
    float copy = it->first;
    if (count > top_copy_count) {
      second_copy_count = top_copy_count;
      second_copy = top_copy;
      top_copy_count = count;
      top_copy = copy;
    } else if (count > second_copy_count) {
      second_copy = copy;
      second_copy_count = count;
    }
  }
  if (top_copy_count > second_copy_count &&
      static_cast<float>(top_copy_count)/
      static_cast<float>(aligned_reads.size()) > MIN_TRAIN_AGREE) {
    *mode = top_copy;
    return true;
  }
  return false;
}

void NoiseModel::FitMutationProb(const vector<AlignedRead>& reads_for_training) {
  // Set up data for logistic regression
  size_t numreads = reads_for_training.size();
  vector<bool> y;
  y.resize(numreads);
  vector<vector<double> >x;
  x.resize(numreads);
  for (size_t i=0; i<numreads; i++) {
    pair<string, int> coord (reads_for_training.at(i).chrom,
			     reads_for_training.at(i).msStart);
    y.at(i) = reads_for_training.at(i).stutter;
    vector<double> xx;
    xx.resize(4);
    xx.at(0) = reads_for_training.at(i).period;
    xx.at(1) = reads_for_training.at(i).msEnd - reads_for_training.at(i).msStart + reads_for_training.at(i).diffFromRef;
    xx.at(2) = strInfo.at(coord).gc;
    xx.at(3) = strInfo.at(coord).score;
    x.at(i) = xx;
  }
  // Run logistic regression
  vector<double> coeffs;
  coeffs.resize(NUM_PARAMS);
  if (logistic_regression(y, x, y.size(), 4,
			  0.001, 1e-5, 1e-5, 1000, &coeffs) != 0) {
    PrintMessageDieOnError("Logistic regression failed", ERROR);
  }
  // Write output file (for now match old format to keep compatbility)
  TextFileWriter nWriter(stutter_model_filename);
  stringstream ss;
  ss << "solver_type L2R_LR" << endl;
  ss << "nr_class 2" << endl;
  ss << "label -1 1" << endl;
  ss << "nr_feature 4" << endl;
  ss << "bias 0" << endl;
  ss << "w" << endl;
  for (int i=0; i<coeffs.size(); i++) {
    ss << coeffs.at(i) << endl;
  }
  nWriter.Write(ss.str());
  if (my_verbose) {
    stringstream msg;
    msg << "Using " << reads_for_training.size() << " reads for training."
         << " (Required: " << MIN_READS_FOR_TRAINING << ")";
    PrintMessageDieOnError(msg.str(), PROGRESS);
  }
  if (reads_for_training.size() < MIN_READS_FOR_TRAINING) {
    PrintMessageDieOnError("Too few reads for training", ERROR);
  }
}

void NoiseModel::FitStepProb(const map<int, map <int,int> > & step_size_by_period) {
  // Fill pdf for observed data
  vector<vector<float> > pdf_obs;
  for (int i = 0; i < MAX_PERIOD-MIN_PERIOD+1; i++) {
    vector<float> v0(0);
    for (int j = 0; j < MAX_PERIOD*3*2+1; j++) {
      v0.push_back(0);
    }
    pdf_obs.push_back(v0);
  }
  // Get avg. step size mod period and prob of incr/decr
  p_incr = 0;
  int total_reads = 0;
  int num_incr = 0;
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    int total_period_reads = 0;
    int total_steps = 0;
    if ((step_size_by_period).find(period) !=
        (step_size_by_period).end()) {
      for (map<int,int>::const_iterator it =
             (step_size_by_period).at(period).begin();
           it != (step_size_by_period).at(period).end(); it++) {
        const int& step = it->first;
        const int& count = it->second;
        const int& modstep = mmod(step,period);
        total_period_reads += count;
        total_steps += modstep*count;
        if (step > 0) {
          num_incr += count;
        }
        pdf_obs.at(period-MIN_PERIOD).at(step+3*MAX_PERIOD) = count;
      }
    }
    if (total_period_reads > 0 && total_steps > 0) {
      float avg_step = static_cast<float>(total_steps)/
        static_cast<float>(total_period_reads);
      average_step_size.push_back(avg_step);
    } else {
      average_step_size.push_back(SMALL_CONST);
    }
    total_reads += total_period_reads;
  }

  //  average_step_size = avg_step_per_period;
  p_incr = static_cast<float>(num_incr)/
    static_cast<float>(total_reads);
  if (debug) {
    stringstream msg;
    msg << "[NoiseModel] num incr " << num_incr << " total " << total_reads
        << "pinr " << p_incr;
    PrintMessageDieOnError(msg.str(), DEBUG);
  }

  // Fill PDF for each period
  // index to element corresponding to 0 diff from ref
  int zero_index = MAX_PERIOD*3;
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    float total_mass = 0;
    for (int i = MAX_PERIOD*-3; i <= MAX_PERIOD*3; i++) {
      float prob = 0;
      if (i > 0 && i <= period*3) {
        int index_plus = zero_index+i;
        int index_minus = zero_index-i;
        // Set poisson probability
        prob = dpois(i, period);
        if (debug) {
          stringstream msg;
          msg << "[NoiseModel] poisson " << prob;
          PrintMessageDieOnError(msg.str(), DEBUG);
        }
        // Set non-unit probability
        prob = prob*dgeom(mmod(i, period), 1/(average_step_size.at(period-MIN_PERIOD)+1));
        if (debug) {
          stringstream msg;
          msg << "[NoiseModel] "
              << "nonunit prob " << prob << "i " 
              << i << " period " << period 
              << " mmod " << mmod(i, period) << " geom " 
              << dgeom(mmod(i, period), 1/(average_step_size.at(period-MIN_PERIOD)+1))
              << " avg stp size " << average_step_size.at(period-MIN_PERIOD);
          PrintMessageDieOnError(msg.str(), DEBUG);
        }
        // Set incr/decr probability
        float prob_plus = prob * p_incr;
        float prob_minus = prob * (1-p_incr);
        pdf_model.at(period-MIN_PERIOD).at(index_plus) = prob_plus;
        pdf_model.at(period-MIN_PERIOD).at(index_minus) = prob_minus;
        total_mass += prob_plus + prob_minus;
      } else {
        pdf_model.at(period-MIN_PERIOD).at(zero_index+i) = 0;
      }
    }
    // Scale
    for (int j = 1; j <= period*3; j++) {
      pdf_model.at(period-MIN_PERIOD).at(zero_index+j) =
        pdf_model.at(period-MIN_PERIOD).at(zero_index+j)/total_mass;
      pdf_model.at(period-MIN_PERIOD).at(zero_index-j) =
	pdf_model.at(period-MIN_PERIOD).at(zero_index-j)/total_mass;
    }
  }

  // Write to file
  TextFileWriter nWriter(stepsize_model_filename);
  stringstream ss;
  // write mean nonunit step size for each period
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    ss << average_step_size[period-MIN_PERIOD] << endl;
  }
  // write prob incr/decr
  ss << "ProbIncrease=" << p_incr << endl;
  // write pdf for model
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    ss << "Period" << period << "Model";
    for (int i = 0; i < MAX_PERIOD*3*2+1; i++) {
      ss << " " << pdf_model.at(period-MIN_PERIOD).at(i);
    }
    ss << endl;
  }
  // write pdf for obs
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    ss << "Period" << period << "Obs";
    for (int i = 0; i < MAX_PERIOD*3*2+1; i++) {
      ss << " " << pdf_obs.at(period-MIN_PERIOD).at(i);
    }
    ss << endl;
  }
  nWriter.Write(ss.str());
}

float NoiseModel::GetTransitionProb(int a, int b,
                                    int period, int length,
                                    float gc, float score) {
  // Outlier reads likely not from stutter
  if (abs(b-a)> 3*period) {
    return 0;
  }
  vector<double> x;
  x.resize(NUM_PARAMS-1);
  x.at(0) = period;
  x.at(1) = length;
  x.at(2) = gc;
  x.at(3) = score;
  float mutProb = logistic_regression_predict(x, stutter_prob_model);
  float lik;
  if (a == b) {
    lik = (1-mutProb);
  } else {
    float stepProb = pdf_model.at(period-MIN_PERIOD).at(b-a + 3*MAX_PERIOD);
    lik = (mutProb)*stepProb;
  }
  return lik;
}

bool NoiseModel::ReadNoiseModelFromFile() {
  // *** Step 1: Read logistic regression model ***
  if (debug) {
    PrintMessageDieOnError("[NoiseModel] Loading stutter prob from file...", DEBUG);
  }
  stutter_prob_model.resize(NUM_PARAMS);
  TextFileReader sFile(stutter_model_filename.c_str());
  string line;
  for (int i=0; i<6; i++) {sFile.GetNextLine(&line);} // header lines
  for (size_t i=0; i<NUM_PARAMS; i++) {
    sFile.GetNextLine(&line);
    stutter_prob_model.at(i) = atof(line.c_str());
  }

  // *** Step 2: Read step size parameters ***
  if (debug) {
    PrintMessageDieOnError("[NoiseModel] Loading step size params from file...", DEBUG);
  }
  TextFileReader nFile(stepsize_model_filename.c_str());
  if (debug) {
    PrintMessageDieOnError("[NoiseModel] Getting average step size...", DEBUG);
  }
  // Get nonunit step size
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    nFile.GetNextLine(&line);
    average_step_size.push_back(atof(line.c_str()));
  }
  if (debug ) {
    PrintMessageDieOnError("[NoiseModel] Getting prob incr/decr...", DEBUG);
  }
  // Get prob incr/decr
  nFile.GetNextLine(&line);
  p_incr = atof(line.c_str());
  // Get pdf
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    nFile.GetNextLine(&line);
    vector<string> items;
    split(line, ' ', items);
    if (items.size() < 3*MAX_PERIOD*2+1) return false;
    for (size_t i = 1; i < items.size(); i++) {
      pdf_model.at(period-MIN_PERIOD).at(i-1) = atof(items[i].c_str());
    }
  }
  return true;
}

NoiseModel::~NoiseModel() {}

