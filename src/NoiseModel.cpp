/*
 * Author: Melissa Gymrek 2012
 */

#include <sstream>
#include <cmath>

#include "NoiseModel.h"
#include "runtime_parameters.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"

const int MIN_PERIOD = 2;
const int MAX_PERIOD = 6;

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

static float invLogit(float x) {
  float expx = exp(x);
  return expx/(1+expx);
}

static float ppois(int step, float mean) {
  if (step < 0) return 0;
  float p = exp(-1*mean);
  for (int i = 0; i < step; i++) {
    p = p*mean;
    p = p/(i+1);
  }
  return p;
}

NoiseModel::NoiseModel(RInside* _r){
  R = _r;
}

void NoiseModel::Train(ReadContainer* read_container) {
  // Populate training data
  if (debug) {
    cerr << "Populating training data..." << endl;
  }
  for (map<pair<string, int>, list<AlignedRead> >::iterator
	 it = read_container->aligned_str_map_.begin();
       it != read_container->aligned_str_map_.end(); it++) {
    // check if haploid
    const list<AlignedRead>& aligned_reads = it->second;
    bool is_haploid = ((aligned_reads.front().chrom == "chrX" ||
			aligned_reads.front().chrom == "chrY" ));
    if (!is_haploid) continue;
    // check if has unique mode
    float unique_mode;
    if (HasUniqueMode(aligned_reads, &unique_mode)) {
      // add to training data (mode, period, copynum, mutated)
      int period = aligned_reads.front().period;
      for (list<AlignedRead>::const_iterator it2 = 
	     aligned_reads.begin(); it2 != aligned_reads.end(); it2++) {
	training_data.mode.push_back(unique_mode);
	training_data.period.push_back(period);
	training_data.copynum.push_back((*it2).diffFromRef);
	training_data.mutated.push_back((*it2).diffFromRef
					!= unique_mode);
      }
    }
  }

  if (debug) {
    cerr << "Converint to R format" << endl;
  }
  // Convert training data vectors to R vectors
  int td_length = training_data.mode.size();
  r_training_data.mode = Rcpp::NumericVector(td_length);
  r_training_data.period = Rcpp::NumericVector(td_length);
  r_training_data.copynum = Rcpp::NumericVector(td_length);
  r_training_data.mutated = Rcpp::NumericVector(td_length);
  for (int i = 0; i < td_length; i++) {
    r_training_data.mode[i] = training_data.mode.at(i);
    r_training_data.period[i] = training_data.period.at(i);
    r_training_data.copynum[i] = training_data.copynum.at(i);
    r_training_data.mutated[i] = training_data.mutated.at(i);
  }
  (*R)["mode"] = r_training_data.mode;
  (*R)["period"] = r_training_data.period;
  (*R)["copynum"] = r_training_data.copynum;
  (*R)["mutated"] = r_training_data.mutated;

  if (debug) {
    cerr << "Fitting models..." << endl;
  }
  // Fit model
  FitMutationProb();
  FitStepProb();
}

bool NoiseModel::HasUniqueMode(const list<AlignedRead>&
			       aligned_reads,
			       float* mode) {
  if (aligned_reads.size() == 1) {
    *mode = aligned_reads.front().diffFromRef;
    return true;
  }
  int top_copy_count;
  float top_copy;
  int second_copy_count;
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
  if (top_copy_count > second_copy_count) {
    *mode = top_copy;
    return true;
  }
  return false;
}

void NoiseModel::FitMutationProb() {
  if (debug) {
    cerr << "Fitting mutation prob..." << endl;
  }
  Rcpp::NumericVector periods(MAX_PERIOD-MIN_PERIOD+1);
  for (int i = 0; i < MAX_PERIOD-MIN_PERIOD+1; i++ ){
    periods[i] = i+MIN_PERIOD;
  }
  (*R)["all_periods"] = periods;

  // Do logistic regression
  string evalString = "mut_model = glm(mutated~period, " \
    "family='binomial'); \n" \
    "mutIntercept = mut_model$coefficients[[1]]; \n" \
    "mutSlope = mut_model$coefficients[[2]]; \n" \
    "c(mutIntercept, mutSlope)"; 
  if (debug) {
    cerr << evalString << endl;
  }
  SEXP ans = (*R).parseEval(evalString); // return intercept, slope
  Rcpp::NumericVector v(ans);
  mutIntercept = v[0];
  mutSlope = v[1];
}

void NoiseModel::FitStepProb() {
  if (debug) {
    cerr << "Fitting step prob..." << endl;
  }
  string evalString = "data = data.frame(mutated=mutated, period = period, copynum = copynum, mode = mode); \n" \
    "mutated_data = data[data$mutated,];\n" \
    "mutated_data$step = abs(mutated_data$copynum - mutated_data$mode); \n" \
    "step_model = glm(mutated_data$step ~ mutated_data$period, family='poisson'); \n" \
    "poisIntercept = step_model$coefficients[[1]];\n" \
    "poisSlope = step_model$coefficients[[2]];\n" \
    "if (is.na(poisSlope)) {poisSlope = 0};\n" \
    "c(poisIntercept, poisSlope)";
  //    "mutated_data$step = round(mutated_data$copynum/mutated_data$period); \n"
  if (debug) {
    cerr << evalString << endl;
  }
  SEXP ans = (*R).parseEval(evalString); // return intercept, slope
  Rcpp::NumericVector v(ans);
  poisIntercept = v[0];
  poisSlope = v[1];
}

float NoiseModel::GetTransitionProb(float a,
				    float b, int period) {
  float mutProb = invLogit(mutIntercept + mutSlope*period);
  float poisMean = exp(poisIntercept + poisSlope*period);
  if (a == b) {
    return mutProb;
  } else {
    int diff = (int)(abs(b-a));
    return (1-mutProb)*ppois(diff-1, poisMean);
  }
  return 0.0;
}

bool NoiseModel::ReadNoiseModelFromFile(string filename) {
  TextFileReader nFile(filename);
  string line;
  bool mI = false;
  bool mS = false;
  bool pI = false;
  bool pS = false;
  while (nFile.GetNextLine(&line)) {
    vector<string> items;
    split(line,'=', items);
    if (debug) {
      cerr << "reading line " << line << " " << items.size() << endl;
    }
    if (items.size() == 0) break;
    if (items.size() != 2) return false;
    if (items.at(0) == "mutIntercept") {
      mutIntercept = atof(items.at(1).c_str());
      mI = true;
    }
    if (items.at(0) == "mutSlope") {
      mutSlope = atof(items.at(1).c_str());
      mS = true;
    }
    if (items.at(0) == "poisIntercept") {
      poisIntercept = atof(items.at(1).c_str());
      pI = true;
    }
    if (items.at(0) == "poisSlope") {
      poisSlope = atof(items.at(1).c_str());
      pS = true;
    }
  }
  if (! (mI && mS && pI && pS)) return false;
  return true;
}

bool NoiseModel::WriteNoiseModelToFile(string filename) {
  TextFileWriter nWriter(filename);
  stringstream ss;
  ss << "mutIntercept=" << mutIntercept << endl
     << "mutSlope=" << mutSlope << endl
     << "poisIntercept=" << poisIntercept << endl
     << "poisSlope=" << poisSlope << endl;
  nWriter.Write(ss.str());
}

NoiseModel::~NoiseModel(){}

