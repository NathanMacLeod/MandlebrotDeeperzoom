#define OLC_PGE_APPLICATION
#include "../include/olcPixelGameEngine.h"
#include "../include/BigFloat.h"
#include <regex>
#include <cassert>
#include <thread>

#define N_THREADS 12
#define DEEP_ZOOM false
static int HASNT_DIVERGED = -1;
static int BLOCK_SIZE = 8;

class Example : public olc::PixelGameEngine
{
public:
	long double camPos[2] = { 0, 0 };
	long double defaultBounds[4] = { -2.5, 1, -1, 1 };

	/*static void threadDrawRow(void* obj, int y, int interval, long double xItrVal, long double yItrVal, long double xLow, long double yLow, long double invZoom) {
		((Example*)obj)->drawRow(y, interval, xItrVal, yItrVal, xLow, yLow, invZoom);
	}*/

	template <class float_type>
	struct cmplx {
		float_type re;
		float_type im;

		cmplx operator*(cmplx c) { return cmplx{ re * c.re - im * c.im, re * c.im + c.re * im }; }
		cmplx operator+(cmplx c) { return cmplx{ re + c.re, im + c.im }; }
		cmplx operator-(cmplx c) { return cmplx{ re - c.re, im - c.im }; }
		float_type mag_sqr() { return re * re + im * im; }
	};

	/*template<int exp_len, int mantissa_len>
	std::vector<cmplx<double>> get_mandlebrot_sequence_of_coord(cmplx<BigFloat<exp_len, mantissa_len>> c, int max_itr) {
		std::vector<cmplx<double>> out(max_itr + 1);
		out[0] = cmplx<double>{ c.re.to_double(), c.im.to_double() };

		cmplx<BigFloat<exp_len, mantissa_len>> x_i = c;
		for (int i = 1; i <= max_itr; i++) {
			x_i = x_i * x_i + c;
			out[i] = cmplx<double>{ x_i.re.to_double(), x_i.im.to_double() };
		}

		return out;
	}*/

	template<int exp_len , int mantissa_len>
	std::vector<cmplx<double>>get_mandlebrot_sequence_of_coord(cmplx<BigFloat<exp_len, mantissa_len>> c, int max_itr) {
		if (c.mag_sqr().to_double() >= 4) return {};
		std::vector<cmplx<double>> out;
		out.push_back(cmplx<double>{c.re.to_double(), c.im.to_double()});

		cmplx<BigFloat<exp_len, mantissa_len>> x_i = c;
		for (int i = 1; i < max_itr && x_i.mag_sqr().to_double() < 4; i++) {
			x_i = x_i * x_i + c;
			out.push_back(cmplx<double>{x_i.re.to_double(), x_i.im.to_double()});
		}

		return out;
	}

	template<class float_type>
	float get_iterations_to_diverge(cmplx<float_type> x_0, int max_itr, bool smoothed=true) {
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
	float continue_iterations_to_diverge(cmplx<float_type> x_i, cmplx<float_type> x_0, int max_itr, bool smoothed=true) {
		int itr = 0;
		while (itr < max_itr && x_i.mag_sqr() < 4) {
			itr++;
			x_i = x_i * x_i + x_0;
		}

		if (x_i.mag_sqr() < 4) return HASNT_DIVERGED;
		return smoothed ? get_smoothed_escape_radius(2, x_i, itr) : itr;
	}

	float get_smoothed_escape_radius(double escape_radius, cmplx<double> first_escaped_pos, int itr_to_escape) {
		double mag = sqrt(first_escaped_pos.mag_sqr());
		return itr_to_escape - log(log(mag) / log(escape_radius)) / log(2);
	}

	// Anchor point strategy comes from SUPERFRACTALTHING MATHS (2013) by K.I. Martin https://web.archive.org/web/20160408070057/http://superfractalthing.co.nf/sft_maths.pdf
	template<int exp_len, int mantissa_len>
	std::vector<float> approx_mandlebrot_itrerations_of_grid(int screen_width, int screen_height, int max_itr, int min_pixel_x, int min_pixel_y, int grid_width, int grid_height, double inv_zoom, cmplx<BigFloat<exp_len, mantissa_len>> camera_coord, bool smoothed=true) {
		int anchor_pixel_x = min_pixel_x + grid_width / 2;
		int anchor_pixel_y = min_pixel_y + grid_height / 2;
		int delta_from_center_x = min_pixel_x - screen_width / 2;
		int delta_from_center_y = min_pixel_y - screen_height / 2;

		cmplx<BigFloat<exp_len, mantissa_len>> anchor_point = camera_coord + cmplx<BigFloat<exp_len, mantissa_len>>{
																				 BigFloat<exp_len, mantissa_len>(delta_from_center_x* inv_zoom), 
																				 BigFloat<exp_len, mantissa_len>(-delta_from_center_y * inv_zoom) }; //neg sign to make screen up positive imaginary axis

		std::vector<cmplx<double>> anchor_sequence = get_mandlebrot_sequence_of_coord(anchor_point, max_itr);
		int n_cells = grid_width * grid_height;
		std::vector<float> n_itr_till_escape_per_cell(n_cells, HASNT_DIVERGED);

		std::vector<cmplx<double>> last_seq_val_per_cell(n_cells);
		cmplx<double> anchor_point_approx = cmplx<double>{ anchor_point.re.to_double(), anchor_point.im.to_double() };
		for (int grid_x = 0; grid_x < grid_width; grid_x++) {
			for (int grid_y = 0; grid_y < grid_height; grid_y++) {
				int indx = grid_width * grid_y + grid_x;
				cmplx<double> start_seq_val = anchor_point_approx + cmplx<double>{ inv_zoom * (grid_x - grid_width / 2), inv_zoom * (grid_width / 2 - grid_y) };
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
			cmplx<double> A_new = cmplx<double>{ 2, 0 } * anchor_sequence[i - 1] * A_i + cmplx<double>{1, 0};
			cmplx<double> B_new = cmplx<double>{ 2, 0 } * anchor_sequence[i - 1] * B_i + A_i * A_i;
			cmplx<double> C_new = cmplx<double>{ 2, 0 } * anchor_sequence[i - 1] * C_i + cmplx<double>{ 2, 0 } * A_i * B_i;

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

					/*cmplx<BigFloat<exp_len, mantissa_len>> real_p = camera_coord + cmplx<BigFloat<exp_len, mantissa_len>>{
																							    BigFloat<exp_len, mantissa_len>(inv_zoom * delta_from_center_x + delta.re),
																							    BigFloat<exp_len, mantissa_len>(delta.im - delta_from_center_y * inv_zoom) };

					cmplx<BigFloat<exp_len, mantissa_len>> real_seq_val = get_mandlebrot_sequence_of_coord(real_p, max_itr)[i];
					cmplx<double> real_seq_val_double = { real_seq_val.re.to_double(), real_seq_val.im.to_double() };
					int a = 1;*/

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
		for (int grid_x = 0; grid_x < grid_width; grid_x++) {
			for (int grid_y = 0; grid_y < grid_height; grid_y++) {
				int indx = grid_width * grid_y + grid_x;
				if (n_itr_till_escape_per_cell[indx] != HASNT_DIVERGED) continue;

				int remaining_iterations = max_itr - num_approx_iterations_completed;
				/*int itr_to_diverge = continue_iterations_to_diverge(
					cmplx<BigFloat<exp_len, mantissa_len>>{ BigFloat<exp_len, mantissa_len>(last_seq_val_per_cell[indx].re), BigFloat<exp_len, mantissa_len>(last_seq_val_per_cell[indx].im) },
					anchor_point + cmplx<BigFloat<exp_len, mantissa_len>>{ inv_zoom* (grid_x - grid_width / 2), inv_zoom* (grid_width / 2 - grid_y)
				}, remaining_iterations);*/
				int itr_to_diverge = continue_iterations_to_diverge(last_seq_val_per_cell[indx], anchor_point_approx + cmplx<double>{ inv_zoom* (grid_x - grid_width / 2), inv_zoom* (grid_width / 2 - grid_y) }, remaining_iterations);
				if (itr_to_diverge != HASNT_DIVERGED) {
					n_itr_till_escape_per_cell[indx] = num_approx_iterations_completed + itr_to_diverge;
				}
			}
		}

		return n_itr_till_escape_per_cell;
	}

	template<class float_type>
	std::vector<float> calculate_mandlebrot_simple(int screen_width, int screen_height, double zoom, cmplx<float_type> camera_coord, int max_itr, bool smoothed =true) {
		std::vector<float> out(screen_width * screen_height);
		double inv_zoom = 1.0 / zoom;
		for (int x = 0; x < screen_width; x++) {
			for (int y = 0; y < screen_height; y++) {

				cmplx<float_type> pixel_delta_from_camera = {
					inv_zoom * (x - screen_width / 2),
					inv_zoom * (screen_height / 2 - y)
				};
				cmplx<float_type> c = camera_coord + pixel_delta_from_camera;
				float itr_to_diverge = get_iterations_to_diverge(c, max_itr);

				int indx = x + y * screen_width;
				out[indx] = itr_to_diverge;

				/*double r, g, b;
				double q = double(itr_to_diverge) / max_itr;

				if (itr_to_diverge == HASNT_DIVERGED) {
					r = g = b = 0;
				}
				else if (q > 0.5) {
					r = 1;
					g = q * inv_zoom;
					b = 1;
				}
				else {
					r = sqrt(q);
					g = 0;
					b = sqrt(q);

				}

				Draw(x, y, olc::Pixel(r * 255, g * 255, b * 255));*/
			}
		}

		return out;
	}

	template<int exp_len, int mantissa_len>
	std::vector<float> calculate_mandlebrot(int screen_width, int screen_height, double zoom, cmplx<BigFloat<exp_len, mantissa_len>> camera_coord, int max_itr, int approx_grid_size) {
		std::vector<float> out(screen_width * screen_height);
		double inv_zoom = 1.0 / zoom;

		for (int x = 0; x < screen_width; x += approx_grid_size) {
			for (int y = 0; y < screen_height; y += approx_grid_size) {

				int grid_width = std::min<int>(approx_grid_size, screen_width - x);
				int grid_height = std::min<int>(approx_grid_size, screen_height - y);
				std::vector<float> itr_vals = approx_mandlebrot_itrerations_of_grid(screen_width, screen_height, max_itr, x, y, grid_width, grid_height, inv_zoom, camera_coord);

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

	double ten_pow_i(int i) {
		double prod = 1;
		for (int j = 0; j < i; j++) prod *= 10;
		return prod;
	}

	template<int exp_len, int mantissa_len>
	std::vector<float> calculate_mandlebrot_dynamic_accuracy(int screen_width, int screen_height, double zoom, cmplx<BigFloat<exp_len, mantissa_len>> camera_coord, int max_itr, int approx_grid_size) {
		std::vector<float> computation;
		if (zoom < ten_pow_i(14)) {
			cmplx<double> double_camera = cmplx<double>{ camera_coord.re.to_double(), camera_coord.im.to_double() };
			computation = calculate_mandlebrot_simple(screen_width, screen_height, zoom, double_camera, max_itr);
			printf("double");
		}
		else if (zoom < ten_pow_i(17)) {
			cmplx<BigFloat<1, 2>> cam = cmplx<BigFloat<1, 2>>{ BigFloat<1, 2>(camera_coord.re), BigFloat<1, 2>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
			printf("prec 2");
		}
		else if (zoom < ten_pow_i(27)) {
			cmplx<BigFloat<1, 3>> cam = cmplx<BigFloat<1, 3>>{ BigFloat<1, 3>(camera_coord.re), BigFloat<1, 3>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
			printf("prec 3");
		}
		else if (zoom < ten_pow_i(37)) {
			cmplx<BigFloat<1, 4>> cam = cmplx<BigFloat<1, 4>>{ BigFloat<1, 4>(camera_coord.re), BigFloat<1, 4>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
			printf("prec 4");
		}
		else if (zoom < ten_pow_i(47)) {
			cmplx<BigFloat<1, 5>> cam = cmplx<BigFloat<1, 5>>{ BigFloat<1, 5>(camera_coord.re), BigFloat<1, 5>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
			printf("prec 5");
		}
		else if (zoom < ten_pow_i(56)) {
			cmplx<BigFloat<1, 6>> cam = cmplx<BigFloat<1, 6>>{ BigFloat<1, 6>(camera_coord.re), BigFloat<1, 6>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
		}
		else if (zoom < ten_pow_i(66)) {
			cmplx<BigFloat<1, 7>> cam = cmplx<BigFloat<1, 7>>{ BigFloat<1, 7>(camera_coord.re), BigFloat<1, 7>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
		}
		else if (zoom < ten_pow_i(76)) {
			cmplx<BigFloat<1, 8>> cam = cmplx<BigFloat<1, 8>>{ BigFloat<1, 8>(camera_coord.re), BigFloat<1, 8>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
		}
		else if (zoom < ten_pow_i(85)) {
			cmplx<BigFloat<1, 9>> cam = cmplx<BigFloat<1, 9>>{ BigFloat<1, 9>(camera_coord.re), BigFloat<1, 9>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
		}
		else {
			cmplx<BigFloat<1, 10>> cam = cmplx<BigFloat<1, 10>>{ BigFloat<1, 10>(camera_coord.re), BigFloat<1, 10>(camera_coord.im) };
			computation = calculate_mandlebrot(screen_width, screen_height, zoom, cam, max_itr, approx_grid_size);
		}

		return computation;
	}

	/*void drawRow(int y, int interval, long double xItrVal, long double yItrVal, long double xLow, long double yLow, long double invZoom) {
		while (y < ScreenHeight()) {
			long double y0 = yLow + y * yItrVal;
			for (int x = 0; x < ScreenWidth(); x++) {
				long double x0 = xLow + x * xItrVal;
				drawMandelBrotPixel(x, y, x0, y0, invZoom);
			}
			y += interval;
		}
	}*/

	/*void drawMandelBrotPixel(int x, int y, BigFloat<2, 3> x0, BigFloat<2, 3> y0, long double invZoom) {
		BigFloat<2, 3> a = x0;
		BigFloat<2, 3> b = y0;
		double maxItr = 1500;
		int i = 0;
		BigFloat<2, 3> divBound = 4;

		while ((a * a + b * b - divBound).get_is_negative() && i < maxItr) {
			i++;
			BigFloat<2, 3> temp = a * a - b * b;
			b = a.getDouble() * b + y0;
			a = temp + x0;
		}

		double quotient = i / maxItr;

		double r;
		double g;
		double blu;

		if (i == maxItr) {
			r = g = blu = 0;
		}
		else if (quotient > 0.5) {
			r = 1;
			g = quotient * invZoom;
			blu = 1;
		}
		else {
			r = sqrt(quotient);
			g = 0;
			blu = sqrt(quotient);

		}

		Draw(x, y, olc::Pixel(r * 255, g * 255, blu * 255));
	}*/

	//void drawMandlebrot(long double camPos[2], long double zoom) {
	//	long double invZoom = 1 / zoom;

	//	long double xLow = camPos[0] + defaultBounds[0] * invZoom;
	//	long double xHigh = camPos[0] + defaultBounds[1] * invZoom;
	//	long double yLow = camPos[1] + defaultBounds[2] * invZoom;
	//	long double yHigh = camPos[1] + defaultBounds[3] * invZoom;

	//	long double xItrVal = (xHigh - xLow) / ScreenWidth();
	//	long double yItrVal = (yHigh - yLow) / ScreenHeight();

	//	//drawRow(0, 1, xItrVal, yItrVal, xLow, yLow);
	//	std::thread* threads[N_THREADS];
	//	for (int i = 0; i < N_THREADS; i++) {
	//		threads[i] = new std::thread(threadDrawRow, this, i, N_THREADS, xItrVal, yItrVal, xLow, yLow, invZoom);
	//	}

	//	for (int i = 0; i < N_THREADS; i++) {
	//		threads[i]->join();
	//		delete threads[i];
	//	}
	//}

	struct rgb {
		double r, g, b;
	};

	void render_calculation_histogram(const std::vector<float>& itr_to_diverge_per_pixel, rgb min_itr_color, rgb max_itr_color) {
		assert(itr_to_diverge_per_pixel.size() == ScreenWidth() * ScreenHeight());

		//using histogram approach described here: https://www.fractalforums.com/programming/newbie-how-to-map-colors-in-the-mandelbrot-set/
		std::map<float, int> num_pixels_at_iteration_count; 
		int total_divergent_pixels = 0;
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
			itr_count_to_gradient_map[itr->first] = double(pixel_sum) / total_divergent_pixels;
		}

		for (int x = 0; x < ScreenWidth(); x++) {
			for (int y = 0; y < ScreenHeight(); y++) {
				int indx = x + ScreenWidth() * y;
				float itr_to_diverge = itr_to_diverge_per_pixel[indx];
				double gradient_value = itr_count_to_gradient_map[itr_to_diverge];

				double r, g, b;
				

				if (itr_to_diverge == HASNT_DIVERGED) {
					r = g = b = 0;
				}
				else {
					r = min_itr_color.r + (max_itr_color.r - min_itr_color.r) * gradient_value;
					g = min_itr_color.g + (max_itr_color.g - min_itr_color.g) * gradient_value;
					b = min_itr_color.b + (max_itr_color.b - min_itr_color.b) * gradient_value;
				}

				Draw(x, y, olc::Pixel(r * 255, g * 255, b * 255));
			}
		}
	}

public:
	Example()
	{
		sAppName = "Mandlebrot";
	}

	//very cool, zooms into spiral and then mini mandlebrot
	//static long double zoom = 199671284473369792;
	//static long double camPos[2]{ -0.8032523144678346, 0.1780442509067320 };

	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		static long double zoom = 1706878340403585575169294336.0;
		//static long double camPos[2]{ -0.7092299821065109, -0.3497044921822606 };
		//static long double camPos[2]{ 0, 0 };
		cmplx<BigFloat<1, 4>> camera_coord = { 
			BigFloat<1, 4>("-5.2331947713954048804324105129573794751258630835645843424439668068239369870002149206876172497521879823953138281122398137766573137334493739068752571341747170663438737392425537109375e-1"),
			BigFloat<1, 4>("-6.881333022784251045352075710453609692937405306722069606536530036149330119054353256161330336414040677611968943591786095878097525763809479037633708153887113212476833723485469818115234375e-1")
		};
		//cmplx<double> camera_coord_doub = { camPos[0], camPos[1] };
		//render_mandlebrot(zoom, camera_coord, 20000, BLOCK_SIZE);
		printf("%s, %s\n", camera_coord.re.to_scientific_notation().c_str(), camera_coord.im.to_scientific_notation().c_str());
		std::vector<float> calc = calculate_mandlebrot(ScreenWidth(), ScreenHeight(), zoom, camera_coord, 41154, BLOCK_SIZE);
		render_calculation_histogram(calc, {0.0f, 0.0f, 0.0f}, { 0.7, 0, 1.0 });
		//render_mandlebrot(ScreenWidth() / 4, cmplx<BigFloat<1, 1>>{0, 0}, 1000, 16);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		static long double zoom = ScreenWidth() / 4;/*149849401626652258854764544.0;*///
		static BigFloat<1, 6> camPos[2]{ -0.7092299821065109, 0.3497044921822606 };
		bool lclick = GetMouse(0).bPressed;
		bool rclick = GetMouse(1).bPressed;
		if (DEEP_ZOOM) {
			zoom += fElapsedTime * zoom / 10;
			int max_itr = 700 * log(zoom / 50);
			cmplx<BigFloat<1, 6>> camera_coord = { camPos[0], camPos[1] };
			std::vector<float> calc = calculate_mandlebrot_dynamic_accuracy(ScreenWidth(), ScreenHeight(), zoom, camera_coord, max_itr, BLOCK_SIZE);
			render_calculation_histogram(calc, { 0.0f, 0.0f, 0.0f }, { 0.7, 0, 1.0 });
			//printf("%.16f, %.16fi, %f, %d\n", camPos[0], camPos[1], zoom, max_itr);
		}
		if (lclick || rclick) {

			if (lclick) {
				camPos[0] += BigFloat<1, 3>(GetMouseX() - ScreenWidth() / 2) / zoom;
				camPos[1] += BigFloat<1, 3>(ScreenHeight() / 2 - GetMouseY()) / zoom;
				zoom *= 1.5;
			}
			else {
				zoom /= 1.5;
			}

			int max_itr = 700 * log(zoom / 50);
			cmplx<BigFloat<1, 6>> camera_coord = { camPos[0], camPos[1] };
			std::vector<float> calc = calculate_mandlebrot_dynamic_accuracy(ScreenWidth(), ScreenHeight(), zoom, camera_coord, max_itr, BLOCK_SIZE);
			render_calculation_histogram(calc, { 0.0f, 0.0f, 0.0f }, { 0.7, 0, 1.0 });
			printf("(%s, %s), %f, %d\n", camPos[0].to_scientific_notation().c_str(), camPos[1].to_scientific_notation().c_str(), zoom, max_itr);
		}
		// called once per frame
		return true;
	}
};


int main()
{
	Example demo;
	if (demo.Construct(300, 200, 3, 3))
		demo.Start();

	return 0;
}
