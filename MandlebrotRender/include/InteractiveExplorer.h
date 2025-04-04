#pragma once
#include "../include/olcPixelGameEngine.h"
#include "../include/Mandlebrot.h"
#include <regex>
#include <cassert>
#include <thread>

class InteractiveExplorer : public olc::PixelGameEngine
{
public:

	void draw_image(image img) {
		for (int x = 0; x < ScreenWidth(); x++) {
			for (int y = 0; y < ScreenHeight(); y++) {
				int indx = x + ScreenWidth() * y;
				rgb color = img.pixels[indx];
				Draw(x, y, olc::Pixel(color.r, color.g, color.b));
			}
		}
	}

	void draw_mandlebrot(cmplx<BigFloat<4>> cam_pos, double zoom, int max_iterations, int supersample_factor, bool smoothing_enabled, int approximation_block_size, const std::vector<rgb>& gradient) {
		int calc_width = ScreenWidth() * supersample_factor;
		int calc_height = ScreenHeight() * supersample_factor;

		std::vector<float> calc;
		if (approximation_block_size > 1) calc = calculate_mandlebrot_dynamic_accuracy(calc_width, calc_height, zoom * supersample_factor, cam_pos, max_iterations, approximation_block_size * supersample_factor, smoothing_enabled);
		else                              calc = calculate_mandlebrot_simple(calc_width, calc_height, zoom * supersample_factor, cam_pos, max_iterations, smoothing_enabled);

		if (supersample_factor != 1) calc = downsize_cacluation(calc, calc_width, calc_height, supersample_factor);
		image img = render_calculation_histogram(ScreenWidth(), ScreenHeight(), calc, gradient);
		draw_image(img);
	}

	int get_dynamic_iteration_count(double zoom) {
		return 700 * log(zoom / 50);
	}

	InteractiveExplorer(render_settings settings)
		: settings(settings)
	{
		sAppName = "Mandlebrot";
	}

	bool OnUserCreate() override
	{
		zoom = ScreenWidth() / 4;
		cam_pos = { 0, 0 };

		draw_mandlebrot(cam_pos, zoom, settings.max_iterations, settings.supersample_factor, settings.smoothing_enabled, settings.approximation_block_size, settings.gradient);
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		bool lclick = GetMouse(0).bPressed;
		bool rclick = GetMouse(1).bPressed;
		if (lclick || rclick) {

			if (lclick) {
				cam_pos.re += BigFloat<3>(GetMouseX() - ScreenWidth() / 2) / zoom;
				cam_pos.im += BigFloat<3>(ScreenHeight() / 2 - GetMouseY()) / zoom;
				zoom *= 1.5;
			}
			else {
				zoom /= 1.5;
			}

			draw_mandlebrot(cam_pos, zoom, settings.max_iterations, settings.supersample_factor, settings.smoothing_enabled, settings.approximation_block_size, settings.gradient);
			printf("(%s, %s), %f\n", cam_pos.re.to_scientific_notation().c_str(), cam_pos.im.to_scientific_notation().c_str(), zoom);
		}
		return true;
	}

private:
	static const int BLOCK_SIZE = 8;
	cmplx<BigFloat<4>> cam_pos;
	double zoom;

	render_settings settings;
};