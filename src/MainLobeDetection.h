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
#ifndef __MAIN_LOBE_DETECTION_H__
#define __MAIN_LOBE_DETECTION_H__

#include <vector>

#include "FFT_nuc_vectors.h"
#include "FFT_four_nuc_vectors.h"


class MainLobeDetection
{
private:
	FFT_FOUR_NUC_VECTORS data;
public:
	double detect_one_window(const std::vector<double> &v, size_t window_size);

	void calculate_lobe_detection_window(const std::string& nucleotides,
			std::vector<double>& /*output*/ lobe_detection_window);
};


#endif
