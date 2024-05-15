#pragma once
#include <string>
#include <vector>
#include "./Mandlebrot.h"

struct option {
	std::string long_name;
	std::string short_name;
	std::string description;
	int n_args;
};

struct option_parse_response {
	std::string option_long_name;
	std::vector<std::string> args;
	int new_pos;
	bool invalid;
};

struct run_info {
	bool use_interactive_mode = true;
	bool parse_error = false;
	render_settings settings;
};

bool arg_is_long_option(std::string arg);

bool arg_is_option(std::string arg);

option_parse_response parse_next_option_and_args(const std::vector<option>& options, int argc, char* argv[], int pos);

run_info parse_run_info_from_args(int argc, char** argv);