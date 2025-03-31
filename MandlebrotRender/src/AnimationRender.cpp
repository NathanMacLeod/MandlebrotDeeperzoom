#include "../include/AnimationRender.h"
#include <thread>
#include <atomic>
#include <format>

int calculate_animation_framecount(double target_zoom, double starting_zoom, int fps, double zoom_speed) {
	//zoom(t) = starting_zoom * (zoom_speed + 1)^time
	double time_to_zoom = log(target_zoom / starting_zoom) / log(zoom_speed + 1);
	return time_to_zoom * fps * time_to_zoom;
}

double get_zoom_of_frame(int frame_number, int fps, double starting_zoom, double zoom_speed) {
	return starting_zoom + pow(zoom_speed + 1, double(frame_number) / fps);
}

void render_frame(int frame_number, const render_settings& settings) {
	double starting_zoom = settings.image_width / 4;
	double zoom = get_zoom_of_frame(frame_number, settings.fps, starting_zoom, settings.zoom_speed);
	int calc_width = settings.image_width * settings.supersample_factor;
	int calc_height = settings.image_height * settings.supersample_factor;

	std::vector<float> calc;
	if (settings.approximation_block_size > 1) {
		calc = calculate_mandlebrot_dynamic_accuracy(
			calc_width, calc_height, zoom * settings.supersample_factor, settings.cam_pos, settings.max_iterations,
			settings.approximation_block_size * settings.supersample_factor, settings.smoothing_enabled
		);
	}
	else {
		calc = calculate_mandlebrot_simple(
			calc_width, calc_height, zoom * settings.supersample_factor, settings.cam_pos, settings.max_iterations, settings.smoothing_enabled
		);
	}

	if (settings.supersample_factor != 1) calc = downsize_cacluation(calc, calc_width, calc_height, settings.supersample_factor);
	image img = render_calculation_histogram(settings.image_width, settings.image_height, calc, settings.gradient);
	std::string path = std::format("{}/{}.png", settings.write_directory_path, frame_number);
	export_iamge(path, img);
}

void multithreaded_render(const render_settings& settings)
{
	double starting_zoom = settings.image_width / 4;
	int final_frame_number = std::min<int>(settings.start_frame + settings.max_frames_to_render, calculate_animation_framecount(settings.target_zoom, starting_zoom, settings.fps, settings.zoom_speed));
	std::atomic_int next_frame = settings.start_frame;

	if (settings.thread_count <= 1) {
		for (int i = 0; i < final_frame_number; i++) {
			render_frame(i, settings);
		}
	}
	else {
		std::vector<std::thread> threads(settings.thread_count);
		for (int i = 0; i < settings.thread_count; i++) {
			threads[i] = std::thread([=, &next_frame]() {
				int my_frame = -1;
				while ((my_frame = std::atomic_fetch_add(&next_frame, 1)) < final_frame_number) {
					printf("Start frame %d/%d\n", my_frame - settings.start_frame, final_frame_number - settings.start_frame);
					render_frame(my_frame, settings);
					printf("Finish frame %d/%d\n", my_frame - settings.start_frame, final_frame_number - settings.start_frame);
				}
				});
		}

		for (std::thread& thread : threads) thread.join();
	}
}