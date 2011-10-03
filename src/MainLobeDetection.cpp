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
#include <algorithm>
#include <numeric>

#include "runtime_parameters.h"

#include "HammingWindowGenerator.h"

#include "MainLobeDetection.h"

#include "common.h"


using namespace std;

double MainLobeDetection::detect_one_window(const std::vector<double> &v, size_t window_size)
{
  size_t lobe_width = round( ((double)v.size())/((double)window_size) );

	if (lobe_debug)
		cerr << "Main-Lobe-Detection ( input-vector-size = " << v.size()
			<< "  window-size = " << window_size
		     << " lobe-width = " << lobe_width << " )" << endl;


	HammingWindowGenerator *pHammingWindowGenerator = HammingWindowGenerator::GetHammingWindowSingleton();
	const WindowVector *pHammingWindow = pHammingWindowGenerator->GetWindow( lobe_width*2 );

	//TODO:
	//multiply the hamming_window vector by the summed_matrix vector (v)

	vector<double> swapped_v(v);
	//swapped_v.resize(v.size());
	//size_t mid_point = v.size()/2;
	//for (size_t i=0; i<mid_point; ++i) {
	//  swapped_v[i] = v[i+mid_point];
	//  swapped_v[mid_point+i] = v[i];
	//}

	if (lobe_debug) {
	  cerr << "Lobe-Debug, doing ffshift(v)" << endl;
	  cerr << "--- matblab eqv. after_fftshift = fftshift(before_fftshift)" << endl;
	  debug_print_matlab_vector(v,"before_fftshift");
	  debug_print_matlab_vector(swapped_v,"after_fftshift");
	}

	//Copy just the two lobes into a new vector
	vector<double> regions_of_interest;
	regions_of_interest.resize(lobe_width*2);
	for (size_t i=0;i<lobe_width;i++) {
	  regions_of_interest[lobe_width + i] = swapped_v[i];
	  regions_of_interest[i] = swapped_v[swapped_v.size()-lobe_width+i];
	}
	if (lobe_debug) {
	  cerr << "--- matlab eqv.   regions_of_interest=after_fftshift XXXXXXX whatever" << endl;
	  debug_print_matlab_vector(regions_of_interest, "regions_of_interest");
	}

	double result = 0;
	for ( size_t i = 0 ; i < regions_of_interest.size(); i++ ) {
	  result += regions_of_interest[i] * (*pHammingWindow)[i];
	}

	double sum = std::accumulate(swapped_v.begin(),swapped_v.end(),0.0);
	double mean = sum / swapped_v.size();

	double noise=0;
	for ( size_t i = 0 ; i<pHammingWindow->size(); ++i) {
	  noise += (*pHammingWindow)[i] * mean ;
	}

	double energy = result / noise;

	if (lobe_debug) {
	  cerr << "Lobe detection results " << endl
	       << "   result = " << result << "\t----- matlab eqv.    roi_y * hamming(" << lobe_width*2 << ")" << endl
	       << "   sum(input) = " << sum << "\t----- matlab eqv.    sum(c_f)" << endl
	       << "   mean(input)= " << mean << "\t----- matlab eqv.     mean(c_f)" << endl
	       << "   noise = " << noise << "\t----- matlab eqv.    noise = sum(mean(c_f)* hamming(" << lobe_width*2 << "))" << endl
	       << "   energy = " << energy << "\t---- matlab eqv.   energy = result / noise" << endl;
	}

	return energy;
}


void MainLobeDetection::calculate_lobe_detection_window(const std::string& nucleotides,
			std::vector<double>& /*output*/ lobe_detection_window)
{
	data.resize(256);
	data.build_plan(256);

	for (size_t i=0;i<nucleotides.length()-fft_window_size;i+=fft_window_step) {
		string window_nucleotides = nucleotides.substr(i, fft_window_size);

		if (fftw_debug || lobe_debug)
			cerr << "Processing window (offset "
				<< i <<" size " << fft_window_size << "): "
					<< window_nucleotides << endl;

		size_t size = window_nucleotides.length(); //should be equal to fft_window_size

		data.set_nucleotides(window_nucleotides);
		data.pad_zeros_to_end(size);

		if (fftw_debug)
			data.debug_print_input(size);

		data.execute();
		if (fftw_debug)
			data.debug_print_output_complex(256);

		data.out_complex_to_magnitude();
		if (fftw_debug)
			data.debug_print_output_magnitude(256);

		data.create_summed_matrix();
		if (fftw_debug)
			data.debug_print_summed_matrix();

		double lobe_detection = detect_one_window(data.summed_matrix, fft_window_size);

		lobe_detection_window.push_back(lobe_detection);
	}
	data.destroy_plan();
}
