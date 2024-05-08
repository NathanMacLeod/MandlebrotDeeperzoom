#include "../include/Mandlebrot.h"
#define OLC_PGE_APPLICATION
#include "../include/InteractiveExplorer.h"
#include "../include/optionparser.h"

enum option_index { UNKNOWN, HELP, INTERACTIVE, RENDER, POS_REAL };

int main() {
	InteractiveExplorer demo;
	if (demo.Construct(300, 200, 3, 3))
		demo.Start();
}