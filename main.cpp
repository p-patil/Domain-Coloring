#include <iostream>
#include <ctime>
#include "complex_numbers.h"
#include "functions.h"
#include "grapher.h"
#include "SDL/SDL.h"

using namespace std;

ComplexNumber identity(ComplexNumber z) { return z; }

ComplexNumber square(ComplexNumber z) { return z * z; }

ComplexNumber exponential(ComplexNumber z) {
	Exponential exp (1, 1);
	return exp.eval(z);
}

ComplexNumber inverse(ComplexNumber z) { return 1 / z; }

ComplexNumber sine(ComplexNumber z) {
	Sine s (1, 1, 0);
	return s.eval(z);
}

ComplexNumber polynomial(ComplexNumber z) { return (z * z * z - 1) * (z - 1); }

int main(int argc, char **argv) {
	srand(time(NULL));

	freopen("CON", "w", stdout); // Prevents redirect of output
	freopen("CON", "w", stderr); // and errors...
	if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
		cerr << "Failed to initialize environment\n";
		return -1;
	}

	SDL_Surface *image = map_function_to_pixels(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, -3, 3, -3, 3, identity);
	display_image(image, 3000);

	SDL_FreeSurface(image);
	SDL_Quit();
	return 0;
}