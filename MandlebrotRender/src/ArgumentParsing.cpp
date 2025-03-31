#include "../include/ArgumentParsing.h"

bool arg_is_long_option(std::string arg) {
	return arg.size() >= 2 && arg[0] == '-' && arg[1] == '-';
}

bool arg_is_option(std::string arg) {
	return arg.size() >= 1 && arg[0] == '-';
}

option_parse_response parse_next_option_and_args(const std::vector<option>& options, int argc, char* argv[], int pos) {
	option_parse_response response = { "", {}, pos, true }; //note response invaid is set to true

	if (pos >= argc) {
		printf("err: no more args left to parse\n");
		return response;
	}
	std::string arg = argv[pos];
	if (!arg_is_option(arg)) {
		printf("err: options arg must begin with -\n");
	}

	bool is_long_option = arg_is_long_option(arg);
	int num_dashes = is_long_option ? 2 : 1;
	arg = arg.substr(num_dashes, arg.size() - num_dashes);
	std::replace(arg.begin(), arg.end(), '-', '_');

	int expected_arg_count = 0;
	for (const option& o : options) {
		if ((is_long_option && o.long_name == arg) || (!is_long_option && o.short_name == arg)) {
			response.option_long_name = o.long_name;
			response.invalid = false;
			expected_arg_count = o.n_args;
		}
	}

	if (response.option_long_name == "") {
		printf("err: option %s not recognized\n", arg.c_str());
		return response;
	}

	pos++;
	for (int i = 0; i < expected_arg_count; i++) {
		if (pos >= argc) {
			printf("err: option %s expected %d arguments, only had %d\n", response.option_long_name.c_str(), expected_arg_count, i);
			response.invalid = true;
			return response;
		}
		response.args.push_back(std::string(argv[pos]));
		pos++;
	}

	response.new_pos = pos;
	return response;
}

static std::vector<option> options = {
	{"help", "h", "Print this message and exit", 0},
	{"real_coordinate", "re", "(Only significant for animation render) real coordinate to center the zoom on.", 1},
	{"imaginary_coordinate", "im", "(Only significant for animation render) imaginary coordinate to center the zoom on.", 1},
	{"target_zoom", "z", "(Only significant for animation render) how deep the animation will zoom to.", 1},
	{"fps", "f", "(Only significant for animation render) framerate for the generated animation.", 1},
	{"output_directory", "o", "(Only significant for animation render) framerate for the generated animation.", 1},
	{"thread_count", "t", "(Only significant for animation render) number of threads to use.", 1},
	{"zoom_speed", "zs", "(Only significant for animation render) How quickly the animation zooms in.", 1},
	{"start_frame", "t", "(Only significant for animation render) Skip rendering the first number of frames in the animation.", 1},
	{"max_frames_rendered", "m", "(Only significant for animation render) Set an upper limit on the total numer of frames to render", 1},
	{"approximation_block_size", "a", "For deeper zooms, uses an approximation method where few BigFloat samples are taken, and then nearby samples within block size are approximated from that sample.", 1},
	{"max_iterations", "i", "The max number of iterations tested when determining if a sequence diverges or not. Higher iterations results in more detail, but longer computation time", 1},
	{"image_width", "w", "Width of the render", 1},
	{"image_height", "h", "Height of the render", 1},
	{"smoothing_enabled", "se", "Enable smooth coloring", 0},
	{"smoothing_disabled", "sd", "disable smooth coloring", 0},
	{"supersample_factor", "ss", "Increases the number of samples per pixel by the square of the factor", 1},
	{"interactive_mode", "im", "Open an interactive explorer. Controls are left-click to zoom, centered at the point clicked, right-click to zoom out. The coordinates and zoom are printed after every click.", 0},
	{"animation_mode", "am", "Render animation frames to output directory specified by --output-directory", 0}
};

run_info parse_run_info_from_args(int argc, char** argv) {
	run_info out;
	out.dont_run = false;

	std::string real_coordinate = "-0.8032523144678346";
	std::string imag_coordinate = "0.1780442509067320";

	int arg_pos = 1;
	while (arg_pos < argc) {
		option_parse_response resp = parse_next_option_and_args(options, argc, argv, arg_pos);

		if (resp.invalid || resp.option_long_name == "help") {
			printf("\n==========MANDLEBROT==RENDERER==========\n");
			printf("Run either an interactive explorer, or render a zoom animation to a specific location.\n");
			printf("Example usage:\nLaunch interactive explorer, rendering with 1000 max iterations:\t./mandlebrot --interactive-mode --max-iterations 1000\n");
			printf("Generate frames for an animation zooming to a specific coordinate:\t./mandlebrot --animation-mode --real-coordinate 1.423212e-2 --imaginary-coordinate 0.456 --target-zoom 1500000000 --image-width 960 --image-height 520 --output-directory ./render_dir\n");

			printf("\n");
			for (option o : options) {
				std::string long_name_flag_form = o.long_name;
				std::replace(long_name_flag_form.begin(), long_name_flag_form.end(), '_', '-');
				printf("--%s; -%s:\t\t%d args. %s\n", long_name_flag_form.c_str(), o.short_name.c_str(), o.n_args, o.description.c_str());
			}
			out.dont_run = true;
			return out;
		}
		else if (resp.option_long_name == "real_coordinate") real_coordinate = resp.args[0];
		else if (resp.option_long_name == "imaginary_coordinate") imag_coordinate = resp.args[0];
		else if (resp.option_long_name == "target_zoom") out.settings.target_zoom = stof(resp.args[0]);
		else if (resp.option_long_name == "fps") out.settings.fps = stoi(resp.args[0]);
		else if (resp.option_long_name == "output_directory") out.settings.write_directory_path = resp.args[0];
		else if (resp.option_long_name == "thread_count") out.settings.thread_count = stoi(resp.args[0]);
		else if (resp.option_long_name == "zoom_speed") out.settings.zoom_speed = stof(resp.args[0]);
		else if (resp.option_long_name == "start_frame") out.settings.start_frame = stoi(resp.args[0]);
		else if (resp.option_long_name == "max_frames_rendered") out.settings.max_frames_to_render = stoi(resp.args[0]);
		else if (resp.option_long_name == "approximation_block_size") out.settings.approximation_block_size = stoi(resp.args[0]);
		else if (resp.option_long_name == "max_iterations") out.settings.max_iterations = stoi(resp.args[0]);
		else if (resp.option_long_name == "image_width") out.settings.image_width = stoi(resp.args[0]);
		else if (resp.option_long_name == "image_height") out.settings.image_height = stoi(resp.args[0]);
		else if (resp.option_long_name == "smoothing_enabled") out.settings.smoothing_enabled = true;
		else if (resp.option_long_name == "smoothing_disabled") out.settings.smoothing_enabled = false;
		else if (resp.option_long_name == "supersample_factor") out.settings.supersample_factor = stoi(resp.args[0]);
		else if (resp.option_long_name == "interactive_mode") out.use_interactive_mode = true;
		else if (resp.option_long_name == "animation_mode") out.use_interactive_mode = false;

		arg_pos = resp.new_pos;
	}

	if (!BigFloat<5>::scientific_notation_number_valid(real_coordinate)) {
		printf("Error: real coodinate %s is not valid notation. Example of a valid number: \"43.5e+15\"\n", real_coordinate);
		out.dont_run = true;
		return out;
	}
	if (!BigFloat<5>::scientific_notation_number_valid(imag_coordinate)) {
		printf("Error: imaginary coodinate %s is not valid notation. Example of a valid number: \"43.5e+15\"\n", imag_coordinate);
		out.dont_run = true;
		return out;
	}

	out.settings.cam_pos = { real_coordinate, imag_coordinate };
	return out;
}