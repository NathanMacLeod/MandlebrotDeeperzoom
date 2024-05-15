#include "../include/Mandlebrot.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

float get_smoothed_escape_radius(double escape_radius, cmplx<double> first_escaped_pos, int itr_to_escape) {
	double mag = sqrt(first_escaped_pos.mag_sqr());
	return itr_to_escape - log(log(mag) / log(escape_radius)) / log(2);
}

double ten_pow_i(int i) {
	double prod = 1;
	for (int j = 0; j < i; j++) prod *= 10;
	return prod;
}

std::vector<float> downsize_cacluation(const std::vector<float>& calc, int calc_width, int calc_height, int shrink_factor) {
	std::vector<float> out;
	for (int y = 0; y < calc_height; y += shrink_factor) {
		for (int x = 0; x < calc_width; x += shrink_factor) {

			float avg = 0;
			int count = 0;

			for (int w = y; w < y + shrink_factor && w < calc_height; w++) {
				for (int u = x; u < x + shrink_factor && u < calc_width; u++) {


					int index = calc_width * w + u;
					avg += calc[index];
					count++;

				}
			}

			avg /= count;
			out.push_back(avg);
		}
	}

	return out;
}

void export_iamge(std::string path, image img) {
	stbi_write_png(path.c_str(), img.width, img.height, 3, img.pixels.data(), sizeof(rgb) * img.width);
}

image render_calculation_histogram(int image_width, int image_height, const std::vector<float>& itr_to_diverge_per_pixel, const std::vector<rgb>& gradient) {
	image out;
	out.width = image_width;
	out.height = image_height;

	if (gradient.size() < 2) {
		printf("Gradient needs to contain at least 2 colors\n");
		return out;
	}

	//using histogram approach described here: https://www.fractalforums.com/programming/newbie-how-to-map-colors-in-the-mandelbrot-set/
	std::map<float, int> num_pixels_at_iteration_count;
	int total_divergent_pixels = 1; //want the numerator in gradient to always be at least 1 smaller than the denominator
	for (float itr_count : itr_to_diverge_per_pixel) {
		if (itr_count == HASNT_DIVERGED) continue;

		total_divergent_pixels++;
		if (num_pixels_at_iteration_count.find(itr_count) == num_pixels_at_iteration_count.end()) {
			num_pixels_at_iteration_count[itr_count] = 1;
		}
		else {
			num_pixels_at_iteration_count[itr_count]++;
		}
	}

	std::map<float, double> itr_count_to_gradient_map;
	int pixel_sum = 0;
	//note iteration will be in sorted order of lowest iteration count to highest iteration count. This is what we want.
	for (auto itr = num_pixels_at_iteration_count.begin(); itr != num_pixels_at_iteration_count.end(); itr++) {
		pixel_sum += itr->second;


		double f = double(pixel_sum) / total_divergent_pixels;
		itr_count_to_gradient_map[itr->first] = pow(double(pixel_sum) / total_divergent_pixels, 5);
	}

	for (int y = 0; y < image_height; y++) {
		for (int x = 0; x < image_width; x++) {

			int indx = x + image_width * y;
			float itr_to_diverge = itr_to_diverge_per_pixel[indx];
			double gradient_value = itr_count_to_gradient_map[itr_to_diverge];

			rgb color;


			if (itr_to_diverge == HASNT_DIVERGED) {
				color = { 0, 0, 0 };
			}
			else {
				double scaled_gradient = gradient_value * (gradient.size() - 1);
				int i = int(scaled_gradient);
				double t = fmod(scaled_gradient, 1.0);

				color.r = gradient[i].r + (gradient[i+1].r - gradient[i].r) * t;
				color.g = gradient[i].g + (gradient[i+1].g - gradient[i].g) * t;
				color.b = gradient[i].b + (gradient[i+1].b - gradient[i].b) * t;
			}

			out.pixels.push_back(color);

			//Draw(x, y, olc::Pixel(color.r, color.b, color.g));
		}
	}

	return out;
}