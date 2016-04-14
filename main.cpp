#include <iostream>
#include <ctime>
#include "complex_numbers.h"
#include "functions.h"
#include "grapher.h"
#include "SDL/SDL.h"

using namespace std;

ComplexNumber square(ComplexNumber);
ComplexNumber identity(ComplexNumber);

int main(int argc, char **argv) {
	srand(time(NULL));

	freopen("CON", "w", stdout); // Prevents redirect of output
	freopen("CON", "w", stderr); // and errors...
	if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
		cerr << "Failed to initialize environment\n";
		return -1;
	}

	SDL_Surface *image = map_function_to_pixels(SCREEN_WIDTH, SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_HEIGHT, SCREEN_BPP, identity);
	// SDL_Surface *image = random_image(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP);
	display_image(image, 2000);

	SDL_FreeSurface(image);
	SDL_Quit();
}

ComplexNumber square(ComplexNumber z) { return z * z; }

ComplexNumber identity(ComplexNumber z) { return z; }