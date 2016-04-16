#ifndef GRAPHER_H
#define GRAPHER_H

#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define SCREEN_BPP 32 // Bits per pixel
#define PIXEL_WIDTH 8
#define CHANNELS 3 // 3 channels for RGB

#include "SDL/SDL.h"
#include "complex_numbers.h"

using namespace std;

class RGBPixel {
	int bits_per_pixel;
	int red, green, blue;

	public:
		RGBPixel(int, int, int, int);

		RGBPixel(int, int);

		int get_bpp(void) const;

		int get_red(void) const;

		int get_green(void) const;

		int get_blue(void) const;

		void set_red(int red);

		void set_green(int green);

		void set_blue(int blue);

		int as_int(void) const;

};

ostream& operator <<(ostream &, const RGBPixel &);

SDL_Surface * random_image(int, int, int);

SDL_Surface * map_function_to_pixels(int, int, int, double, double, double, double, function<ComplexNumber(ComplexNumber)>);

void display_image(SDL_Surface *, long);

// HELPER FUNCTIONS

RGBPixel get_pixel32(const SDL_Surface *, int, int);

void set_pixel32(SDL_Surface *, int, int, RGBPixel);

RGBPixel numberToPixel(ComplexNumber, double);

double * hsl_to_rgb(double h, double s, double l);

#endif