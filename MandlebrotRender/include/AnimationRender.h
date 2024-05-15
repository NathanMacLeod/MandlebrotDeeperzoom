#pragma once
#include "Mandlebrot.h"

int calculate_animation_framecount(double target_zoom, double starting_zoom, int fps, double zoom_speed);

double get_zoom_of_frame(int frame_number, int fps, double starting_zoom, double zoom_speed);

void render_frame(int frame_number, const render_settings& settings);

void multithreaded_render(const render_settings& settings);