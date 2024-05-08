#pragma once
#include "BigFloat.h"
#include <vector>
#include <map>

template <class float_type>
struct cmplx {
	float_type re;
	float_type im;

	cmplx operator*(cmplx c) { return cmplx{ re * c.re - im * c.im, re * c.im + c.re * im }; }
	cmplx operator+(cmplx c) { return cmplx{ re + c.re, im + c.im }; }
	cmplx operator-(cmplx c) { return cmplx{ re - c.re, im - c.im }; }
	float_type mag_sqr() { return re * re + im * im; }
};

static int HASNT_DIVERGED = -1;

struct rgb {
	uint8_t r, g, b;
};

struct image {
	unsigned int height, width;
	std::vector<rgb> pixels;
};

template<int mantissa_len>
static std::vector<cmplx<double>>get_mandlebrot_sequence_of_coord(cmplx<BigFloat<mantissa_len>> c, int max_itr) {
	if (c.mag_sqr().to_double() >= 4) return {};
	std::vector<cmplx<double>> out;
	out.push_back(cmplx<double>{c.re.to_double(), c.im.to_double()});

	cmplx<BigFloat<mantissa_len>> x_i = c;
	for (int i = 1; i < max_itr && x_i.mag_sqr().to_double() < 4; i++) {
		x_i = x_i * x_i + c;
		out.push_back(cmplx<double>{x_i.re.to_double(), x_i.im.to_double()});
	}

	return out;
}

template<class float_type>
static float get_iterations_to_diverge(cmplx<float_type> x_0, int max_itr, bool smoothed) {
	int itr = 0;
	cmplx<float_type> x_i = x_0;

	while (itr < max_itr && x_i.mag_sqr() < 4) {
		itr++;
		x_i = x_i * x_i + x_0;
	}

	if (x_i.mag_sqr() < 4) return HASNT_DIVERGED;
	return smoothed ? get_smoothed_escape_radius(2, x_i, itr) : itr;
}

template<class float_type>
static float continue_iterations_to_diverge(cmplx<float_type> x_i, cmplx<float_type> x_0, int max_itr, bool smoothed) {
	int itr = 0;
	while (itr < max_itr && x_i.mag_sqr() < 4) {
		itr++;
		x_i = x_i * x_i + x_0;
	}

	if (x_i.mag_sqr() < 4) return HASNT_DIVERGED;
	return smoothed ? get_smoothed_escape_radius(2, x_i, itr) : itr;
}

template<int len>
float get_smoothed_escape_radius(double escape_radius, cmplx<BigFloat<len>> first_escaped_pos, int itr_to_escape) {
	double mag = sqrt(first_escaped_pos.mag_sqr().to_double());
	return itr_to_escape - log(log(mag) / log(escape_radius)) / log(2);
}

float get_smoothed_escape_radius(double escape_radius, cmplx<double> first_escaped_pos, int itr_to_escape);

// Anchor point strategy comes from SUPERFRACTALTHING MATHS (2013) by K.I. Martin https://web.archive.org/web/20160408070057/http://superfractalthing.co.nf/sft_maths.pdf
template<int mantissa_len>
static std::vector<float> approx_mandlebrot_itrerations_of_grid(int screen_width, int screen_height, int max_itr, int min_pixel_x, int min_pixel_y, int grid_width, int grid_height, double inv_zoom, cmplx<BigFloat<mantissa_len>> camera_coord, bool smoothed) {
	int anchor_pixel_x = min_pixel_x + grid_width / 2;
	int anchor_pixel_y = min_pixel_y + grid_height / 2;
	int delta_from_center_x = min_pixel_x - screen_width / 2;
	int delta_from_center_y = min_pixel_y - screen_height / 2;

	cmplx<BigFloat<mantissa_len>> anchor_point = camera_coord + cmplx<BigFloat<mantissa_len>>{
		BigFloat<mantissa_len>(delta_from_center_x* inv_zoom),
			BigFloat<mantissa_len>(-delta_from_center_y * inv_zoom) }; //neg sign to make screen up positive imaginary axis

	std::vector<cmplx<double>> anchor_sequence = get_mandlebrot_sequence_of_coord(anchor_point, max_itr);
	int n_cells = grid_width * grid_height;
	std::vector<float> n_itr_till_escape_per_cell(n_cells, HASNT_DIVERGED);

	std::vector<cmplx<double>> last_seq_val_per_cell(n_cells);
	cmplx<double> anchor_point_approx = cmplx<double>{ anchor_point.re.to_double(), anchor_point.im.to_double() };
	for (int grid_x = 0; grid_x < grid_width; grid_x++) {
		for (int grid_y = 0; grid_y < grid_height; grid_y++) {
			int indx = grid_width * grid_y + grid_x;
			cmplx<double> start_seq_val = anchor_point_approx + cmplx<double>{ inv_zoom* (grid_x - grid_width / 2), inv_zoom* (grid_width / 2 - grid_y) };
			last_seq_val_per_cell[indx] = start_seq_val;
			if (start_seq_val.mag_sqr() >= 4) {
				n_itr_till_escape_per_cell[indx] = 0;
			}
		}
	}

	cmplx<double> A_i = { 1, 0 };
	cmplx<double> B_i = { 0, 0 };
	cmplx<double> C_i = { 0, 0 };
	int num_approx_iterations_completed = 0;
	for (int i = 1; i < anchor_sequence.size(); i++) {
		cmplx<double> A_new = cmplx<double>{ 2, 0 } *anchor_sequence[i - 1] * A_i + cmplx<double>{1, 0};
		cmplx<double> B_new = cmplx<double>{ 2, 0 } *anchor_sequence[i - 1] * B_i + A_i * A_i;
		cmplx<double> C_new = cmplx<double>{ 2, 0 } *anchor_sequence[i - 1] * C_i + cmplx<double>{ 2, 0 } *A_i* B_i;

		A_i = A_new;
		B_i = B_new;
		C_i = C_new;

		double k = grid_width * inv_zoom;
		if (C_i.mag_sqr() * k * k * k * 1 > B_i.mag_sqr() * k * k) {
			//no longer within accuracy threshold
			break;
		}
		num_approx_iterations_completed++;

		for (int grid_x = 0; grid_x < grid_width; grid_x++) {
			for (int grid_y = 0; grid_y < grid_height; grid_y++) {

				int indx = grid_width * grid_y + grid_x;
				if (n_itr_till_escape_per_cell[indx] != HASNT_DIVERGED) continue;

				cmplx<double> approx_seq_val;
				cmplx<double> delta = { inv_zoom * (grid_x - grid_width / 2), inv_zoom * (grid_width / 2 - grid_y) };
				cmplx<double> delta_sqr = delta * delta;
				approx_seq_val = anchor_sequence[i] + (A_i * delta + ((B_i * delta_sqr) + (C_i * delta * delta_sqr)));
				last_seq_val_per_cell[indx] = approx_seq_val;

				if (min_pixel_x + grid_x == 63 && min_pixel_y + grid_y == 100) {
					double a_mag = sqrt((A_i * delta).mag_sqr());
					double b_mag = sqrt((B_i * delta_sqr).mag_sqr());
					double c_mag = sqrt((C_i * delta_sqr * delta).mag_sqr());

					int abr = 1 + 2;
				}

				double mag = approx_seq_val.mag_sqr();
				if (approx_seq_val.mag_sqr() >= 4) {
					n_itr_till_escape_per_cell[indx] = smoothed ? get_smoothed_escape_radius(2, approx_seq_val, i) : i;
				}
			}
		}
	}

	//handle calculating the rest if anchor sequence terminated before this sequence
	//doesn't work very well with using double precision at a certain depth for some reason. Can be "fixed" by swtiching to contiuing the computation
	//with a BigFloat type, but that really tanks the performance.
	for (int grid_x = 0; grid_x < grid_width; grid_x++) {
		for (int grid_y = 0; grid_y < grid_height; grid_y++) {
			int indx = grid_width * grid_y + grid_x;
			if (n_itr_till_escape_per_cell[indx] != HASNT_DIVERGED) continue;

			int remaining_iterations = max_itr - num_approx_iterations_completed;
			float itr_to_diverge = continue_iterations_to_diverge(last_seq_val_per_cell[indx], anchor_point_approx + cmplx<double>{ inv_zoom* (grid_x - grid_width / 2), inv_zoom* (grid_width / 2 - grid_y) }, remaining_iterations, smoothed);
			if (itr_to_diverge != HASNT_DIVERGED) {
				n_itr_till_escape_per_cell[indx] = num_approx_iterations_completed + itr_to_diverge;
			}
		}
	}

	return n_itr_till_escape_per_cell;
}

template<class float_type>
std::vector<float> calculate_mandlebrot_simple(int screen_width, int screen_height, double zoom, cmplx<float_type> camera_coord, int max_itr, bool smoothed) {
	std::vector<float> out(screen_width * screen_height);
	double inv_zoom = 1.0 / zoom;
	for (int x = 0; x < screen_width; x++) {
		for (int y = 0; y < screen_height; y++) {

			cmplx<float_type> pixel_delta_from_camera = {
				inv_zoom * (x - screen_width / 2),
				inv_zoom * (screen_height / 2 - y)
			};
			cmplx<float_type> c = camera_coord + pixel_delta_from_camera;
			float itr_to_diverge = get_iterations_to_diverge(c, max_itr, smoothed);

			int indx = x + y * screen_width;
			out[indx] = itr_to_diverge;
		}
	}

	return out;
}

template<int mantissa_len>
std::vector<float> calculate_mandlebrot(int screen_width, int screen_height, double zoom, cmplx<BigFloat<mantissa_len>> camera_coord, int max_itr, int approx_grid_size, bool smoothed) {
	std::vector<float> out(screen_width * screen_height);
	double inv_zoom = 1.0 / zoom;

	for (int x = 0; x < screen_width; x += approx_grid_size) {
		for (int y = 0; y < screen_height; y += approx_grid_size) {

			int grid_width = std::min<int>(approx_grid_size, screen_width - x);
			int grid_height = std::min<int>(approx_grid_size, screen_height - y);
			std::vector<float> itr_vals = approx_mandlebrot_itrerations_of_grid(screen_width, screen_height, max_itr, x, y, grid_width, grid_height, inv_zoom, camera_coord, smoothed);

			for (int u = 0; u < grid_width; u++) {
				for (int w = 0; w < grid_height; w++) {
					int itr_val = itr_vals[u + w * grid_width];

					int indx = x + u + (y + w) * screen_width;
					out[indx] = itr_val;
				}
			}
		}
	}

	return out;
}

double ten_pow_i(int i);

template<int mantissa_len>
std::vector<float> calculate_mandlebrot_dynamic_accuracy(int screen_width, int screen_height, double zoom, cmplx<BigFloat<mantissa_len>> camera_coord, int max_itr, int approx_grid_size, bool smoothed) {
	std::vector<float> computation;
	if (zoom < ten_pow_i(14)) {
		cmplx<double> double_camera = cmplx<double>{ camera_coord.re.to_double(), camera_coord.im.to_double() };
		computation = calculate_mandlebrot_simple(screen_width, screen_height, zoom, double_camera, max_itr, smoothed);
		printf("double");
	}
	else if (zoom < ten_pow_i(17)) {
		cmplx<BigFloat<2>> cam = cmplx<BigFloat<2>>{ BigFloat<2>(camera_coord.re), BigFloat<2>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
		printf("prec 2");
	}
	else if (zoom < ten_pow_i(27)) {
		cmplx<BigFloat<3>> cam = cmplx<BigFloat<3>>{ BigFloat<3>(camera_coord.re), BigFloat<3>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
		printf("prec 3");
	}
	else if (zoom < ten_pow_i(37)) {
		cmplx<BigFloat<4>> cam = cmplx<BigFloat<4>>{ BigFloat<4>(camera_coord.re), BigFloat<4>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
		printf("prec 4");
	}
	else if (zoom < ten_pow_i(47)) {
		cmplx<BigFloat<5>> cam = cmplx<BigFloat<5>>{ BigFloat<5>(camera_coord.re), BigFloat<5>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
		printf("prec 5");
	}
	else if (zoom < ten_pow_i(56)) {
		cmplx<BigFloat<6>> cam = cmplx<BigFloat<6>>{ BigFloat<6>(camera_coord.re), BigFloat<6>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
	}
	else if (zoom < ten_pow_i(66)) {
		cmplx<BigFloat<7>> cam = cmplx<BigFloat<7>>{ BigFloat<7>(camera_coord.re), BigFloat<7>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
	}
	else if (zoom < ten_pow_i(76)) {
		cmplx<BigFloat<8>> cam = cmplx<BigFloat<8>>{ BigFloat<8>(camera_coord.re), BigFloat<8>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
	}
	else if (zoom < ten_pow_i(85)) {
		cmplx<BigFloat<9>> cam = cmplx<BigFloat<9>>{ BigFloat<9>(camera_coord.re), BigFloat<9>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
	}
	else {
		cmplx<BigFloat<10>> cam = cmplx<BigFloat<10>>{ BigFloat<10>(camera_coord.re), BigFloat<10>(camera_coord.im) };
		computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size, smoothed);
	}

	return computation;
}

std::vector<float> downsize_cacluation(const std::vector<float>& calc, int calc_width, int calc_height, int shrink_factor);

image render_calculation_histogram(int image_width, int image_height, const std::vector<float>& itr_to_diverge_per_pixel, rgb min_itr_color, rgb max_itr_color);

void export_iamge(std::string path, image img);