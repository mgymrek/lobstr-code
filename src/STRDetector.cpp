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

#include <fftw3.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <string>
#include <vector>

#include "src/common.h"
#include "src/EntropyDetection.h"
#include "src/FFT_four_nuc_vectors.h"
#include "src/HammingWindowGenerator.h"
#include "src/ISatellite.h"
#include "src/runtime_parameters.h"
#include "src/STRDetector.h"
#include "src/TukeyWindowGenerator.h"

using namespace std;

/*
   FFT-with-FFTW copied from http://nashruddin.com/fft-with-fftw-example.html
*/

STRDetector::STRDetector() {}

bool STRDetector::ProcessReadPair(ReadPair* read_pair, string* err, string* messages) {
  *err = "Detection-errors-here:";
  *messages = "Detection-notes-here:";
  read_pair->read1_passed_detection = false;
  read_pair->read2_passed_detection = false;
  if (ProcessRead(&read_pair->reads.at(0), err, messages)) {
    read_pair->read1_passed_detection = true;
  }
  // Returns true if at least one read in the pair is detected
  if (read_pair->reads.at(0).paired) {
    if (ProcessRead(&read_pair->reads.at(1), err, messages)) {
      read_pair->read2_passed_detection = true;
    }
    return (read_pair->read1_passed_detection ||
            read_pair->read2_passed_detection);
  } else {
    return read_pair->read1_passed_detection;
  }
}

bool STRDetector::ProcessRead(MSReadRecord* read, string* err, string* messages) {
  // Get the size of the read so we don't keep computing
  size_t read_length = read->nucleotides.size();
  // Preprocessing checks
  if (read_length < (fft_window_size-1) ||
      calculate_N_percentage(read->nucleotides) > percent_N_discard) {
    if (debug) {
      *err += "failed-read-length-check";
    }
    return false;
  }

  vector<double> window_lobe_detection;

  //  Step 1 - sliding window
  size_t nuc_start;
  size_t nuc_end;
  string detected_microsatellite_nucleotides = "";
  size_t next_best_period;
  EntropyDetection ed_filter(read->nucleotides,
                             fft_window_size, fft_window_step);
  if (!ed_filter.EntropyIsAboveThreshold()) {
    if (debug) {
      stringstream msg;
      msg << "failed-entropy-threshold-" << ed_filter.GetMaxEntropy();
      *err += msg.str();
    }
    return false;
  }
  size_t start, end;
  bool rep_end = false;
  ed_filter.FindStartEnd(&start, &end, &rep_end);

  nuc_start = start * fft_window_step;  // + extend_flank;
  nuc_end = (end + 2) * fft_window_step;  // - extend_flank;

  if (nuc_start >= read_length ||
      nuc_start <= 0 ||
      nuc_end-nuc_start + 1 <= 0 ||
      nuc_start >= nuc_end ||
      nuc_end-nuc_start+1 >= read_length ||
      nuc_end >= read_length) {
    if (debug) {
      *err += "failed-STR-region-location-sanity-check";
    }
    return false;
  }

  // allow more nucleotides in detection step
  if ((nuc_start-extend_flank >= 0) &&
      (nuc_start-extend_flank < read_length) &&
      (nuc_end + extend_flank + 1 < read_length)) {
    detected_microsatellite_nucleotides =
      read->nucleotides.substr(nuc_start-extend_flank,
                               nuc_end - nuc_start+1+2*extend_flank);
  } else {
    detected_microsatellite_nucleotides =
      read->nucleotides.substr(nuc_start, nuc_end - nuc_start+1);
  }

  size_t best_period = 0;
  int min_period_to_try = min_period == 1 ? 2 : min_period;

  // if contains > 80% of same nucleotide and has stretch of more than 6, use that
  string over_abundant_nuc =
    OneAbundantNucleotide(detected_microsatellite_nucleotides, 0.9);
  if (over_abundant_nuc != "" &&
      CountAbundantNucRuns(detected_microsatellite_nucleotides, over_abundant_nuc.at(0)) > 6) {
    if (debug) {
      *messages += "Passed-homopolymer-one-abundant-nuc-check;";
    }
    best_period = 1;
  } else {
    // Step 3 - Do FFT on the newly extract microsatellite region
    TukeyWindowGenerator *pTukeyWindowGenerator =
      TukeyWindowGenerator::GetTukeyWindowSingleton();
    const WindowVector *pTukeyWindow =
      pTukeyWindowGenerator->GetWindow
      (detected_microsatellite_nucleotides.length());
    FFT_FOUR_NUC_VECTORS ms;
    ms.resize(1024);
    ms.build_plan(1024);
    ms.set_nucleotides(detected_microsatellite_nucleotides);
    ms.multiply_nuc_matrix_by_vector(*pTukeyWindow);
    ms.pad_zeros_to_end(detected_microsatellite_nucleotides.length());
    ms.execute();
    ms.out_complex_to_magnitude();
    ms.create_summed_matrix();
    ms.destroy_plan();

    // Step 4 Period Detection
    std::vector<double> &fft_vec = ms.summed_matrix;
    size_t lobe_width =
      round(1024.0 / static_cast<double>
            (detected_microsatellite_nucleotides.length()));
    HammingWindowGenerator *pHammingWindowGenerator =
      HammingWindowGenerator::GetHammingWindowSingleton();
    const WindowVector *pHamming_Noise =
      pHammingWindowGenerator->GetWindow(lobe_width);
    double noise_y = 0;
    std::vector<double> noise;
    noise.resize(lobe_width);
    for (size_t i = 0; i < noise.size(); ++i) {
      noise[i] = fft_vec [ rand()%1024 ];  // NOLINT
      noise_y += noise[i] *  pHamming_Noise->at(i);
    }
    std::vector<double> energy;
    const WindowVector *pHamming_roi =
      pHammingWindowGenerator->GetWindow(lobe_width * 2 + 1);
    for (size_t period = min_period_to_try;
         period <= max_period_to_try; period++) {
      double period_energy = 0;
      size_t f = 1;
      double center = 1024 * static_cast<double>(f)/
        static_cast<double>(period);
      double roi_y = 0;
      size_t roi_x_start = round(center - static_cast<double>(lobe_width));
      size_t roi_x_end = round(center + static_cast<double>(lobe_width)+1);
      if (roi_x_end > 1023)
        roi_x_end = 1023;
      for (size_t roi_x = roi_x_start; roi_x < roi_x_end; roi_x++) {
        roi_y += fft_vec[roi_x] * pHamming_roi->at(roi_x-roi_x_start);
      }
      period_energy = (roi_y - noise_y) * period;
      energy.push_back(period_energy);
    }

    // Find the max. period
    std::vector<double>::iterator max_it =
      max_element(energy.begin(), energy.end());
    size_t max_index = distance(energy.begin(), max_it);
    best_period = max_index + min_period_to_try;
    float best_energy = energy.at(max_index);
    float next_best_energy = 0;
    for (size_t newperiod = min_period_to_try; newperiod <= max_period;
         newperiod++) {
      if (newperiod != best_period) {
        if ((energy.at(newperiod - min_period_to_try))>next_best_energy) {
          next_best_energy = energy.at(newperiod - min_period_to_try);
        }
      }
    }

    // check that energy is high enough
    if (best_energy < period_energy_threshold) {
      if (debug) {
        stringstream msg;
        msg << "failed-FFT-energy-threshold-" << best_energy << ";";
        *err += msg.str();
      }
      return false;
    }
    if (next_best_period == 4 && best_period == 2) {
      best_period = 4;
    }
  }

  // Step 6 - Store values in the Read Record
  read->ms_start = nuc_start;
  read->ms_end   = nuc_end;
  read->ms_repeat_best_period = best_period;

  // set indices of left, STR, and right regions
  if ((read->ms_start >=
       static_cast<int>(read_length)) ||
      (read->ms_end+1 >=
       static_cast<int>(read_length))) {
    if (debug) {
      *err += "failed-second-STR-region-location-sanity-check;";
    }
    return false;
  }

  read->left_flank_nuc = (nuc_start>0) ?
    read->nucleotides.substr(0, read->ms_start) : "-";
  read->right_flank_nuc = read->nucleotides.substr(read->ms_end+1);
  read->detected_ms_region_nuc = detected_microsatellite_nucleotides;

  // adjust for max flank region lengths,if repetitive end, don't trim
  read->left_flank_index_from_start = 0;
  read->right_flank_index_from_end = 0;
  if ((read->left_flank_nuc.size() > max_flank_len)) {
    string left_flank = read->left_flank_nuc;
    read->left_flank_nuc = left_flank.substr(left_flank.length()-
                                             max_flank_len, max_flank_len);
    read->left_flank_index_from_start = left_flank.length() -
      max_flank_len;
  }
  if ((read->right_flank_nuc.size() > max_flank_len)) {
    string right_flank = read->right_flank_nuc;
    read->right_flank_nuc = right_flank.substr(0, max_flank_len);
    read->right_flank_index_from_end = right_flank.length() - max_flank_len;
  }

  size_t nuc_len = (read->orig_nucleotides.length() -
                    read->right_flank_index_from_end -
                    read->left_flank_index_from_start);

  if (((read->left_flank_index_from_start+nuc_len) <=
       read_length) &
      ((read->left_flank_index_from_start+nuc_len) >= 0)) {
    read->nucleotides = read->
      orig_nucleotides.substr(read->left_flank_index_from_start, nuc_len);
    read_length = read->nucleotides.size();
    read->quality_scores = read->
      quality_scores.substr(read->left_flank_index_from_start,
                            read_length);
  } else {
    if (debug) {
      *err += "failed-after-trimming-flanking-regions;";
    }
    return false;
  }
  if ( ((read->left_flank_nuc.length() >= min_flank_len) &
        (read->right_flank_nuc.length() >= min_flank_len) &
        (best_period <= max_period) &
        (best_period >= min_period) &
        (read->detected_ms_region_nuc.length() >= MIN_STR_LENGTH))) {
    string repseq = "";
    string repseq_error = "";
    string second_best_repseq = "";
    if (best_period == 1) {
      repseq = getFirstString(over_abundant_nuc, reverseComplement(over_abundant_nuc));
    } else if (!getMSSeq(read->detected_ms_region_nuc,
                         read->ms_repeat_best_period, &repseq, &second_best_repseq, &repseq_error)) {
      if (debug) {
        *err += "failed-STR-when-detecting-motif:" + repseq_error;
      }
      return false;
    }
    read->repseq = repseq;
    return true;
  } else {
    if (debug) {
      *err += "failed-period-and-flanking-length-checks;";
    }
    return false;
  }
}

