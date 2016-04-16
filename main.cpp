#include <iostream>
#include <vector>
#include <functional>
#include "complex_numbers.h"
#include "functions.h"
#include "grapher.h"
#include "parser.h"
#include "SDL/SDL.h"

using namespace std;
using namespace std::placeholders;

SDL_Surface * graph(const string, double = -3, double = 3, double = -3, double = 3); // Default range is -3 <= x, y <= 3

int main(int argc, char **argv) {
	freopen("CON", "w", stdout); // Prevents redirect of output
	freopen("CON", "w", stderr); // and errors...
	if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
		cerr << "Failed to initialize environment\n";
		return -1;
	}

	bool flag = true;
	while (flag) {
		try {
			flag = false;
			string expr;
			cout << "Enter an expression: ";
			cin >> expr;

			SDL_Surface *image = graph(expr); // Get the graph, as a surface.
			display_image(image, 2000); // Display for 5000 milliseconds
		} catch (const exception &e) {
			flag = true;
			cout << "Invalid expression. ";
		}
	}

	SDL_Quit();
	return 0;
}

SDL_Surface * graph(const string s, double x_min, double x_max, double y_min, double y_max) {
	ExpressionTreeNode *expr_tree = parse(s);
	function<ComplexNumber(ComplexNumber)> f = bind(evaluate_tree, expr_tree, _1);
	return map_function_to_pixels(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, x_min, x_max, y_min, y_max, f);
}