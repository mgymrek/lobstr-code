/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

 This file is part of MicroSatelliteDetector.

 MicroSatelliteDetector is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MicroSatelliteDetector is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MicroSatelliteDetector.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "ISatellite.h"
#include "TestFftwSatellite.h"

#include <fftw3.h>

#include "runtime_parameters.h"

#include "MainLobeDetection.h"
#include "TukeyWindowGenerator.h"
#include "HammingWindowGenerator.h"
#include "MicroSatelliteDetection.h"

#include "common.h"
#include "runtime_parameters.h"

using namespace std;


TestFftwSatellite::TestFftwSatellite()
{
}


/*
   FFT-with-FFTW copied from http://nashruddin.com/fft-with-fftw-example.html
 */


bool TestFftwSatellite::ProcessRead(ReadRecord& read)
{

	if (read.nucleotides.length() < (fft_window_size-1)) {
	  return false;
	}

	if (calculate_N_percentage(read.nucleotides) > percent_N_discard) {
	  // cerr << "Discarding read '" << read.ID << "' (too many Ns)" << endl;
	  return false;
	}

	if (read.nucleotides.length() < min_read_length) {
	  return false;
	  cerr << "Discarding read '" << read.ID << "' (read length too short)" << endl;
	}
	// check if meets quality threshold
	//if (avgQuality(read.quality_scores) < quality_control){return false;}

	if (fftw_debug || lobe_debug || microsatellite_detection_debug)
		cerr << "Processing read: " << read.nucleotides << endl;

	vector<double> window_lobe_detection;

	/*
	  Step 1 -
		For each window (of fft_window_size & fft_window_step ):
			1. Calculate FFT on nucleotides matrix
			2. Sum up the matrix
			3. Detect Lobe on summed matrix
			4. Store results in window_lobes_detection vector
	*/
	MainLobeDetection mld;
	mld.calculate_lobe_detection_window(read.nucleotides, window_lobe_detection);

	/*
	   Step 2 - is this read a MicroSatellite ?
	 */
	if (lobe_debug || microsatellite_detection_debug)
		debug_print_matlab_vector(window_lobe_detection, "window_lobe_detection_vector");

	MicroSatelliteDetection ms_detect(window_lobe_detection);

	if (!ms_detect.is_above_threshold()) {
		//Not above threshold, bail out
		if (microsatellite_detection_debug)
			cerr << "MicroSatellite-Detection: window-lobe-vector doesn't cross threshold ("
				<< lobe_detection_threshold << ") - discarding this read." << endl;
		return false;
	}


	size_t lobe_vector_start;
	size_t lobe_vector_end;
	ms_detect.find_start_end(lobe_vector_start,lobe_vector_end);

	size_t nuc_start = lobe_vector_start * fft_window_step ; // add mgymrek to add fft_window_step
	size_t nuc_end   = lobe_vector_end * fft_window_step + fft_window_size ; // add mgymrekt o subtract fft_window_step
	string detected_microsatellite_nucleotides = read.nucleotides.substr(nuc_start,nuc_end - nuc_start+1);

	if (microsatellite_detection_debug) {
		cerr << "MicroSatellite-Detection: " << endl
			<< "Max. Value = " << ms_detect.find_max_lobe_detection_value() << endl
			<< "Lobe-Vector (start = " << lobe_vector_start << "  end = " << lobe_vector_end << ")" << endl
			<< "Nucleotides ( start = " << nuc_start << "  end = " << nuc_end << ") = " << detected_microsatellite_nucleotides << endl ;
	}

	if (ms_detect.is_whole_read_microsatellite(lobe_vector_start, lobe_vector_end)) {
		//The entire read is most likely a microsatellite

		//In DEBUG mode - print warning and continue.
		//In non-debug mode - discard this read
		if (microsatellite_detection_debug)
			cerr << "MicroSatellite-Detection: The entire read is a MicroSatellite" << endl;
		else {
		  //return false;
		}
	}

	/*
	   Step 3 -
		Do FFT on the newly extract microsatellite region
	 */
	TukeyWindowGenerator *pTukeyWindowGenerator = TukeyWindowGenerator::GetTukeyWindowSingleton();
	const WindowVector *pTukeyWindow = pTukeyWindowGenerator->GetWindow( detected_microsatellite_nucleotides.length() ) ;

	FFT_FOUR_NUC_VECTORS ms;
	ms.resize(1024);
	ms.build_plan(1024);
	ms.set_nucleotides(detected_microsatellite_nucleotides);
	if (microsatellite_detection_debug) {
	  cerr << "MicroSatellite-detection, read = " << detected_microsatellite_nucleotides << endl;
	  cerr << "Matrix before Tukey multiplication" << endl;
	  ms.debug_print_input(detected_microsatellite_nucleotides.length());
	}
	if (pTukeyWindow->size() > ms.A.vector_size) return false;
	ms.multiply_nuc_matrix_by_vector(*pTukeyWindow);
	if (microsatellite_detection_debug) {
	  cerr << "Matrix after Tukey multiplication" << endl;
	  ms.debug_print_input(detected_microsatellite_nucleotides.length());
	}
	ms.pad_zeros_to_end(detected_microsatellite_nucleotides.length());
	ms.execute();
	ms.out_complex_to_magnitude();
	ms.create_summed_matrix();

	if (microsatellite_detection_debug) {
	  cerr << "MicroSatellite-detection, read = " << detected_microsatellite_nucleotides << endl;
	  ms.debug_print_input(detected_microsatellite_nucleotides.length());
		ms.debug_print_output_complex(1024);
		ms.debug_print_output_magnitude(1024);
		cerr << "MicroSatellite-detection, FFT on detected region:" << endl;
		debug_print_matlab_vector(ms.summed_matrix,"fft_ms");
	}

	ms.destroy_plan();

	/*
	   Step 4 -
	   Period Detection
	 */
	size_t best_period = 0;
	size_t next_best_period = 0;
	string ms_period_nuc;

	char over_abundtant_nuc = OneAbundantNucleotide(detected_microsatellite_nucleotides);

	if (over_abundtant_nuc != 0) {
		//One nucleotide occupies 80% or more of the microsatellite - discard it

		if ( period_detection_debug ) {
			cerr << "MicroSatellite-detection - Over abundant nucleotide found in detected region:" << endl
			<< "    " << detected_microsatellite_nucleotides << endl
			<< "  returning best_period=1 and seq = " << ms_period_nuc << endl;


		}
		best_period = 1 ;
		ms_period_nuc = over_abundtant_nuc;
	} else {

		/* Do period detection */


		//Same "fft_vec" as in the Matlab "detect_period" function
		std::vector<double> &fft_vec = ms.summed_matrix;

		size_t start_period=1; // mgymrek changed from 1 made these runtime params
		//size_t stop_period =6; // mgymrek changed from 9
		size_t lobe_width = round( 1024.0 / ((double)detected_microsatellite_nucleotides.length()) );
		HammingWindowGenerator *pHammingWindowGenerator = HammingWindowGenerator::GetHammingWindowSingleton();
		const WindowVector *pHamming_Noise = pHammingWindowGenerator->GetWindow( lobe_width );
		double noise_y = 0;
		std::vector<double> noise;
		noise.resize(lobe_width);
		for (size_t i = 0 ; i<noise.size(); ++i) {
			noise[i] = fft_vec [ rand()%1024 ] ;
			noise_y += noise[i] *  pHamming_Noise->at(i);
		}
		if (force_noise_y) {
			noise_y = force_noise_y_value; //debugging - force noise_y
			cerr << "DEBUG!! - forcing noise_y to " << force_noise_y << endl;
		}

		if ( period_detection_debug ) {
			cerr << "Period-Detection:" << endl
				<< "\tlobe-width = " << lobe_width << endl
				<< "\tstart_period = " << start_period << endl
				<< "\tstop_period = " << stop_period << endl
				<< "\tnoise_y  = " << noise_y << endl;
		}

		std::vector<double> energy;

		const WindowVector *pHamming_roi = pHammingWindowGenerator->GetWindow( lobe_width * 2 + 1 );

		// 'period' is 't' in the matlab "detect_period" function
		for (size_t period = start_period; period<=stop_period; period++) {

			double period_energy = 0 ;

			size_t f = 1 ;
			double center = 1024 * (double)f / (double)period;

			double roi_y = 0 ;
			size_t roi_x_start = round(center - (double)lobe_width);
			size_t roi_x_end = round(center + (double)lobe_width+1) ;
			if (roi_x_end>1023)
				roi_x_end=1023;
			size_t roi_x_size = roi_x_end - roi_x_start ;
			if ( period_detection_debug ) {
				cerr << "Period " << period << endl;
				cerr << "roi_x_size = " << roi_x_size << endl;
				cerr << "DEBUG - roi_x_start = " << roi_x_start
					<< " roi_x_end = " << roi_x_end
					<< " roi_x_size = " << roi_x_size
					<< " max size = " << pHamming_roi->size() << endl;
			}
			for ( size_t roi_x = roi_x_start ; roi_x < roi_x_end ; roi_x++ ) {
		//		cerr << "fft_vec[roi_x=" << roi_x<<"] = " << fft_vec[roi_x]
		//			<< " hamming offset = " << roi_x - roi_x_start <<endl;
				roi_y += fft_vec[roi_x] * pHamming_roi->at(roi_x-roi_x_start);
			}

			//In matlab period_detect, lines 32/33
			period_energy = (roi_y - noise_y) * period;

			if ( period_detection_debug ) {
				cerr << "Period-dection (period = " << period << "):" << endl
					<< "\tcenter = " << center << endl
					<< "\troi_x = " << (roi_x_start) << " to " << (roi_x_end) << endl
					<< "\troi_y = " << roi_y << endl
					<< "\tperiod_energy = " << period_energy << endl;
			}
			energy.push_back(period_energy);
		}

		//Find the max. period
		std::vector<double>::iterator max_it = max_element(energy.begin(),energy.end());
		size_t max_index = distance(energy.begin(), max_it);

		best_period = max_index + start_period;
		float best_energy = energy.at(max_index);
		float next_best_energy = 0;
		
		

		
		for (size_t newperiod = start_period; newperiod<=stop_period; newperiod++) {
		  if (newperiod != best_period){
		    if ((energy.at(newperiod - start_period))>next_best_energy){
		      next_best_energy = energy.at(newperiod - start_period);
		      next_best_period = newperiod;
		    }
		  }
		}
		
		if(fabs(next_best_energy-best_energy)/((best_energy+next_best_energy)/2)  <=closeness & next_best_period != 0){
		  //cerr << best_period << " " << next_best_period << " " << next_best_energy << " "<< best_energy << " " << read.nucleotides << endl;
		}else{ next_best_period = 0;}
		
		if ( period_detection_debug ) {
			debug_print_matlab_vector(energy,"energy");
			cerr << "Period-detection results:" << endl
				<< "max_index = " << max_index << endl
				<< "best_period = " << best_period << endl;
		}


	}


	/*
	   Step 6 - Store values in the Read Record
	*/
	read.ms_start = nuc_start;
	read.ms_end   = nuc_end;
	read.ms_repeat_best_period = best_period;
	read.ms_repeat_next_best_period = next_best_period;
	read.left_flank_nuc = (nuc_start>0) ? read.nucleotides.substr(0, read.ms_start) : "-" ;
	read.right_flank_nuc = read.nucleotides.substr(read.ms_end+1); // changed mgymrek
	read.detected_ms_region_nuc = detected_microsatellite_nucleotides;
	read.detected_ms_nuc = ms_period_nuc;
 if((read.left_flank_nuc.length() > min_flank_length) & (read.right_flank_nuc.length() > min_flank_length)){
	//return TRUE to output this read, FALSE to discard it
	return true;
 }else{ return false;}
}

