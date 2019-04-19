/*
This file is part of UMAP.  For copyright information see the COPYRIGHT
file in the top level directory, or at
https://github.com/LLNL/umap/blob/master/COPYRIGHT
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License (as published by the Free
Software Foundation) version 2.1 dated February 1999.  This program is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the terms and conditions of the GNU Lesser General Public License
for more details.  You should have received a copy of the GNU Lesser General
Public License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include "fitsio.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utility/commandline.hpp"
#include "../utility/umap_fits_file.hpp"
#include "../utility/time.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include "cube.hpp"

using namespace median;

using pixel_type = float;
constexpr size_t default_num_random_vector = 0;

void map_fits(const std::string &filename,
	size_t *size_x,
	size_t *size_y,
	size_t *size_k,
	pixel_type **image_data) {
	size_t byte_per_element;
	// Map FITS files using UMap
	*image_data = (pixel_type *)utility::umap_fits_file::PerFits_alloc_cube(filename, &byte_per_element,
		size_x, size_y, size_k);

	if (*image_data == nullptr) {
		std::cerr << "Failed to allocate memory for cube" << std::endl;
		std::abort();
	}

	if (sizeof(pixel_type) != byte_per_element) {
		std::cerr << "Pixel type is not float" << std::endl;
		std::abort();
	}
}

std::size_t get_num_vectors() {
	std::size_t num_random_vector = default_num_random_vector;
	const char *buf = std::getenv("NUM_VECTORS");
	if (buf != nullptr) {
		num_random_vector = std::stoll(buf);
	}
	return num_random_vector;
}

std::vector<double> read_exposuretime(const size_t size_k) {
	std::vector<double> exposuretime_list;

	const char *exposuretime_file_name = std::getenv("EXPOSURETIME_FILE");
	if (exposuretime_file_name != nullptr) {
		std::ifstream ifs(exposuretime_file_name);
		if (!ifs.is_open()) {
			std::cerr << "Cannot open " << exposuretime_file_name << std::endl;
			std::abort();
		}
		for (double exposuretime; ifs >> exposuretime;) {
			exposuretime_list.emplace_back(exposuretime);
		}
		if (exposuretime_list.size() != size_k) {
			std::cerr << "#of lines in " << exposuretime_file_name << " is not the same as #of fits files" << std::endl;
			std::abort();
		}
	}
	else {
		// If a list of exposure times is not given, assume that each exposure time is 40 s
		exposuretime_list.resize(size_k);
		for (size_t i = 0; i < size_k; ++i) exposuretime_list[i] = 40;
	}

	return exposuretime_list;
}

// Function to read data info from a csv file
// Reads timestamp, psf fwhm, (ra/dec), and background sky noise 
std::tuple<std::vector<unsigned long>, std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> read_list_csv(const size_t size_k) {
	std::vector<double> psf_list;
	std::vector<unsigned long> timestamp_list;
	std::vector<std::vector<double>> ra_dec_list;
	std::vector<double> noise_list;

	const char *list_file_name = std::getenv("DATA_LIST_FILE");
	if (list_file_name != nullptr) {
		std::ifstream ifs(list_file_name);

		if (!ifs.is_open()) {
			std::cerr << "Cannot open " << list_file_name << std::endl;
			std::abort();
		}

		std::string line;
		std::getline(ifs, line); //skip first row of header info

		int check = 0;
		double psf, mjd, ra, dec, noise;
		char vala, valb, valc, vald, vale;
		std::string frame;
		double mjd_start = 0;

		while (std::getline(ifs, line))
		{
			std::istringstream iss{ line };

			std::getline(iss, frame, ',');
			do {

				if (iss >> psf >> valb >> mjd >> valc >> ra >> vald >> dec >> vale >> noise)
				{

					if (check == 0)
						mjd_start = mjd;
					unsigned long time = std::round((mjd - mjd_start) * 24 * 60 * 60 * 100); // Hundreths of a second
					std::vector<double> ra_dec = { ra, dec };

					if (psf == 0) // For nan rows we set values to the previous row
					{
						time = timestamp_list[(check - 1)];
						psf = psf_list[(check - 1)];
						ra_dec = ra_dec_list[(check - 1)];
						noise = noise_list[(check - 1)];
					}

					timestamp_list.push_back(time);
					psf_list.push_back(psf);
					ra_dec_list.push_back(ra_dec);
					noise_list.push_back(noise);
					++check;
				}
			} while (!iss.eof());
		}

		if (psf_list.size() != size_k) {
			std::cerr << "#of lines in " << list_file_name << " is not the same as #of fits files" << std::endl;
			std::abort();
		}
	}
	else {
		timestamp_list.resize(size_k);
		psf_list.resize(size_k);
		ra_dec_list.resize(size_k);
		noise_list.resize(size_k);
		for (size_t i = 0; i < size_k; ++i) {
			timestamp_list[i] = (double)i*100.0;
			psf_list[i] = 2.0;
			ra_dec_list[i] = { 0.0,0.0 };
			noise_list[i] = 100.0;
		}
	}

	return std::make_tuple(timestamp_list, psf_list, ra_dec_list, noise_list);
}


std::vector<double> pixel_to_equ(std::vector<double> xy_pos) {
	double delta_x = (-0.00007325)*(xy_pos[0] - 78000);
	double delta_y = (0.00007325)*(xy_pos[1] - 78000);
	std::vector<double> radec = { delta_x + 139.489333678528 , delta_y + 15.7289075993474 };

	return radec;
}

std::vector<double> equ_to_pixel(std::vector<double> ra_dec) {
	double delta_x = (xy_pos[0] - 139.489333678528)/(-0.00007325);
	double delta_y = (xy_pos[1] - 15.7289075993474)/(0.00007325);
	std::vector<double> xy_pos = { delta_x + 78000 , delta_y + 78000 };

	return xy_pos;
}

// Function to calculate relevant information about a given vector
// Returns: <SNR, weighted sum, number of frames intersected>
template <typename iterator_type>
std::tuple<double, typename iterator_type::value_type, int>
vector_info(iterator_type iterator_begin, iterator_type iterator_end) {
	using value_type = typename iterator_type::value_type;

	if (iterator_begin == iterator_end)
		return std::tuple<double, value_type, int>(0, 0, 0);

	// DECAM info
	double dark_noise = 0.417; // electrons per pixel per second
	double readout_noise = 7; // electrons

	value_type total_signal = 0;
	double total_B = 0;
	double total_R = 0;
	double total_D = 0;
	int frame_num = 0;

	for (auto iterator(iterator_begin); iterator != iterator_end; ++iterator) {
		std::tuple<pixel_type, int, double, double> snr_info = iterator.snr_info();
		const value_type value = std::get<0>(snr_info);
		int num_pixels = std::get<1>(snr_info);

		if (num_pixels == 0) continue; // For when the vector hits nothing in an image

		total_signal += value;
		++frame_num;

		// SNR calculation
		double B = std::get<3>(snr_info); // pull background noise from list
		double exp_time = std::get<2>(snr_info);

		total_B += B * num_pixels * exp_time;
		total_R += num_pixels * readout_noise*readout_noise;
		total_D += dark_noise * num_pixels * exp_time;
	}

	double SNR = 0;
	if ((total_signal > 0) && (frame_num != 0))
		SNR = total_signal / sqrt(total_signal + total_B + total_D + total_R);

	return std::tuple<double, value_type, int>(SNR, total_signal, frame_num);
}


std::vector<vector_xy> read_vector_catalog(const cube<pixel_type> &cube, const char *catalog_filename, double pixel_scale) {
	
	/*
	There is definitely a way to read in the fits table row by row with a for loop
	and convert to vector_xy on the spot and insert in-place into the output vector,
	but I can't figure it out right now.
	*/

	std::vector<vector_xy> result;

	fitsfile *fptr;

	int status, anynull;
	float floatnull = 0.;
	long nrows;

	status = 0;
	
	if (fits_open_table(&fptr, catalog_filename, READONLY, &status)) {
		fits_report_error(stderr, status);
		exit(status);
	}
	
	fits_get_num_rows(fptr, &nrows, &status); //get number of rows to use as number of elements per column
	
	//std::vector<double> dec(nrows), v_ra(nrows), v_dec(nrows), time(nrows);
	double ra[100000], dec[100000], v_ra[100000], v_dec[100000], time[100000];
	
	fits_read_col(fptr, TDOUBLE, 1 , 1, 1, nrows, &floatnull, &v_ra, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 2 , 1, 1, nrows, &floatnull, &v_dec, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 3 , 1, 1, nrows, &floatnull, &time, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 4 , 1, 1, nrows, &floatnull, &ra, &anynull, &status);
	fits_read_col(fptr, TDOUBLE, 5 , 1, 1, nrows, &floatnull, &dec, &anynull, &status);
	
	for (int ii = 0; ii < nrows; ii++) {
		//NEED FUNCTIONS HERE FOR CONVERSION

		double x_slope, y_slope, x_int, y_int;

		x_slope = -v_ra[ii] / pixel_scale;
		y_slope = v_dec[ii] / pixel_scale;

		std::vector<double> radec = { ra[ii],dec[ii] };

		std::vector<double> position = some_func(radec);

		x_int = position[0] - x_slope * time[ii];
		y_int = position[1] - y_slope * time[ii];

		vector_xy vec{ x_slope,x_int,y_slope,y_int };

		result.push_back(vec);

	}
	
	return result;
}

std::pair<double, std::vector<std::tuple<vector_xy, double, double, int>>> shoot_vector_catalog(const cube<pixel_type> &cube, const std::vector<vector_xy> vector_catalog, const std::size_t num_vec) {
	// Array to store results of the median calculation
	std::vector<std::tuple<vector_xy, double, double, int>> result;

	double total_execution_time = 0.0;


	

		// Shoot random vectors using multiple threads

		for (int ii = 0; ii < ((num_vec == 0) ? (nrows) : (num_vec)); ii++) {
			vector_xy current_vector = vector_catalog[ii];
			cube_iterator_with_vector<pixel_type> begin(cube, current_vector, 0.0);
			cube_iterator_with_vector<pixel_type> end(cube, current_vector);

			// vector info stored as [VECTOR_XY, SNR, SUM, NUMBER OF FRAMES]
			const auto start = utility::elapsed_time_sec();
			std::tuple<double, double, int> v_info = vector_info(begin, end);
			result.push_back(std::make_tuple(current_vector, std::get<0>(v_info), std::get<1>(v_info), std::get<2>(v_info)));
			total_execution_time += utility::elapsed_time_sec(start);
		}
	

	return std::make_pair(total_execution_time, result);
}


// Function to write results to a csv file in the form:
// ID | X_INTERCEPT | Y_INTERCEPT | X_SLOPE | Y_SLOPE | SNR | SUM | NUMBER OF FRAMES HIT
void write_tocsv(std::vector<std::tuple<vector_xy, double, double, int>> &result) {

	double snr_cutoff = 0;
	std::ofstream out("vector_catalog_output.csv");

	out << "ID,X_INTERCEPT,Y_INTERCEPT,X_SLOPE,Y_SLOPE,SNR,SUM,NUMBER_OF_FRAMES_HIT\n";

	long long id = 1;
	for (auto& row : result) {

		if (std::get<1>(row) <= snr_cutoff) {
			++id;
			continue;
		}
		out << id << ',';
		out << std::get<0>(row).x_intercept << ',';
		out << std::get<0>(row).y_intercept << ',';
		out << std::get<0>(row).x_slope << ',';
		out << std::get<0>(row).y_slope << ',';
		out << std::get<1>(row) << ',';
		out << std::get<2>(row) << ',';
		out << std::get<3>(row) << ',';
		out << '\n';
		++id;
	}
}

int main(int argc, char **argv) {
	utility::umt_optstruct_t options;
	umt_getoptions(&options, argc, argv);

#ifdef _OPENMP
	omp_set_num_threads(options.numthreads);
#endif

	size_t size_x; size_t size_y; size_t size_k;
	pixel_type *image_data;
	map_fits(options.filename, &size_x, &size_y, &size_k, &image_data);
	std::tuple<std::vector<unsigned long>, std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> lists = read_list_csv(size_k);
	cube<pixel_type> cube(size_x, size_y, size_k, image_data, std::get<0>(lists), read_exposuretime(size_k), std::get<1>(lists), std::get<2>(lists), std::get<3>(lists));

	const std::size_t num_random_vector = get_num_vectors();


	double pixel_scale = 0.27;
	
	const char *catalog_filename = std::getenv("VECTOR_CATALOG");
	std::vector<vector_xy> vector_catalog = read_vector_catalog(cube, catalog_filename, pixel_scale);

	auto result = shoot_vector_catalog(cube, vector_catalog, num_random_vector);




	//std::cout << "#of vectors = " << vector_catalog.size()
	//	<< "\nexecution time (sec) = " << result.first
	//	<< "\nvectors/sec = " << static_cast<double>(vector_catalog.size()) / result.first << std::endl;

	//write_tocsv(result.second);

	utility::umap_fits_file::PerFits_free_cube(image_data);

	return 0;
}
