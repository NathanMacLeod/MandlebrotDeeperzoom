#pragma once
#include "BigFloat.h"
#include <vector>

template <class float_type>
struct cmplx {
	float_type re;
	float_type im;

	cmplx operator*(cmplx c) { return cmplx{ re * c.re - im * c.im, re * c.im + c.re * im }; }
	cmplx operator+(cmplx c) { return cmplx{ re + c.re, im + c.im }; }
	cmplx operator-(cmplx c) { return cmplx{ re - c.re, im - c.im }; }
	float_type mag_sqr() { return re * re + im * im; }
};

struct rgb {
	uint8_t r, g, b;
};

struct image {
	unsigned int height, width;
	std::vector<rgb> pixels;
};

template<int mantissa_len>
std::vector<cmplx<double>>get_mandlebrot_sequence_of_coord(cmplx<BigFloat<mantissa_len>> c, int max_itr);

template<class float_type>
float get_iterations_to_diverge(cmplx<float_type> x_0, int max_itr, bool smoothed = true);

template<class float_type>
float continue_iterations_to_diverge(cmplx<float_type> x_i, cmplx<float_type> x_0, int max_itr, bool smoothed = true);

float get_smoothed_escape_radius(double escape_radius, cmplx<double> first_escaped_pos, int itr_to_escape);

// Anchor point strategy comes from SUPERFRACTALTHING MATHS (2013) by K.I. Martin https://web.archive.org/web/20160408070057/http://superfractalthing.co.nf/sft_maths.pdf
template<int mantissa_len>
std::vector<float> approx_mandlebrot_itrerations_of_grid(int screen_width, int screen_height, int max_itr, int min_pixel_x, int min_pixel_y, int grid_width, int grid_height, double inv_zoom, cmplx<BigFloat<mantissa_len>> camera_coord, bool smoothed = true);

template<class float_type>
std::vector<float> calculate_mandlebrot_simple(int screen_width, int screen_height, double zoom, cmplx<float_type> camera_coord, int max_itr, bool smoothed = true);

template<int mantissa_len>
std::vector<float> calculate_mandlebrot(int screen_width, int screen_height, double zoom, cmplx<BigFloat<mantissa_len>> camera_coord, int max_itr, int approx_grid_size);

template<int mantissa_len>
std::vector<float> calculate_mandlebrot_dynamic_accuracy(int screen_width, int screen_height, double zoom, cmplx<BigFloat<mantissa_len>> camera_coord, int max_itr, int approx_grid_size);

std::vector<float> downsize_cacluation(const std::vector<float>& calc, int calc_width, int calc_height, int shrink_factor);

image render_calculation_histogram(int image_width, int image_height, const std::vector<float>& itr_to_diverge_per_pixel, rgb min_itr_color, rgb max_itr_color);