#include "../include/Mandlebrot.h"
#define OLC_PGE_APPLICATION
#include "../include/InteractiveExplorer.h"
#include "../include/AnimationRender.h"
#include "../include/ArgumentParsing.h"

int main(int argc, char* argv[])
 {
	run_info info = parse_run_info_from_args(argc, argv);
	if (info.parse_error) return 1;

	if (info.use_interactive_mode) {
		InteractiveExplorer demo(info.settings);
		if (demo.Construct(300, 200, 3, 3))
			demo.Start();
	}
	else {
		multithreaded_render(info.settings);
	}
}