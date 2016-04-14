#include <iostream>
#include <stdexcept>
#include <SDL/SDL.h>
#include <cstdlib>
#include <limits>
#include <cmath>
#include "grapher.h"
#include "complex_numbers.h"

using namespace std;

RGBPixel::RGBPixel(int bpp, int r, int g, int b) {
	if (bpp < 0 || bpp > 32) {
		cerr << "Invalid bits per pixel\n";
		exit(1);
	}

	this->bits_per_pixel = bpp;
	this->red = r;
	this->green = g;
	this->blue = b;
}

RGBPixel::RGBPixel(int bpp, int pixel_as_int) {
	if (32 % bpp != 0 || bpp > 32 / 3) {
		cerr << "Invalid bits per pixel: " << bpp;
		exit(1);
	}

	int blue_mask = ~((-1) << bits_per_pixel);
	int green_mask = blue_mask << bits_per_pixel;
	int red_mask = blue_mask << (bits_per_pixel * 2);

	this->bits_per_pixel = bpp;
	this->red = pixel_as_int & red_mask;
	this->green = pixel_as_int & green_mask;
	this->blue = pixel_as_int & blue_mask;
}

int RGBPixel::get_bpp(void) const { return this->bits_per_pixel; }

int RGBPixel::get_red(void) const { return this->red; }

int RGBPixel::get_green(void) const {	return this->green; }

int RGBPixel::get_blue(void) const { return this->blue; }

void RGBPixel::set_red(int red) { this->red = red; }

void RGBPixel::set_green(int green) { this->green = green; }

void RGBPixel::set_blue(int blue) { this->blue = blue; }

int RGBPixel::as_int(void) const {
	if (32 % bits_per_pixel != 0) {
		return -1;
	}

	return (red << (bits_per_pixel * 2)) | (green << bits_per_pixel) | blue;
}

RGBPixel numberToPixel(ComplexNumber z, double maximum) {
	RGBPixel p (8, 0, 0, 0);
	double hue, saturation, lightness;

	hue = z.angle() / (2 * PI);
	saturation = 0.5;
	lightness = z.mag() / maximum;

	int *rgb = hslToRgb(hue, saturation, lightness);

	p.set_red(rgb[0]);
	p.set_green(rgb[1]);
	p.set_blue(rgb[2]);

	free(rgb);
	return p;
}

// Override
ostream& operator <<(ostream &stream, const RGBPixel &p) {
	stream << "(" << p.get_red() << ", " << p.get_green() << ", " << p.get_blue() << ")";
	return stream;
}

// Returns the pixel at coordinate (i, j) on the given surface. Assumes width, height, depth > 0 and depth <= 4.
RGBPixel get_pixel32(const SDL_Surface *surface, int x, int y) {
	// Retrieve pixel, as a 32-bit int.
	int pixel_as_int = ((int *) surface->pixels)[x + (y * surface->w)];

	// Create and return an RGBPixel.
	RGBPixel p (PIXEL_WIDTH, pixel_as_int);
	return p;
}

// Sets the pixel at coordinate (i, j) on the given surface to the given pixel.
void set_pixel32(SDL_Surface *surface, int x, int y, RGBPixel pixel) {
	// Get the pixel array.
	int *pixels = (int *) surface->pixels;

	// Write the given RGBPixel, as a 32-bit int, to the array.
	pixels[x + (y * surface->w)] = pixel.as_int();
}

// Randomly initializes pixels. Useful for testing purposes.
SDL_Surface * random_image(int width, int height, int bpp) {
	int blue_mask = ~((-1) << bpp); // Red mask
	int green_mask = blue_mask << bpp; // Green mask
	int red_mask = blue_mask << (bpp * 2); // Blue mask

	SDL_Surface *image = SDL_CreateRGBSurface(SDL_SWSURFACE, // Flag specifying to create image in memory
											  width, // Width of image
											  height, // Height of image
											  bpp, // Bits per pixel
											  red_mask,
											  green_mask,
											  blue_mask,
											  0); // No alpha mask

	if (SDL_MUSTLOCK(image)) {
		SDL_LockSurface(image);
	}

	int rgb_limit = 1 << PIXEL_WIDTH;
	RGBPixel p (PIXEL_WIDTH, 0, 0, 0);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			p.set_red(rand() % rgb_limit);
			p.set_green(rand() % rgb_limit);
			p.set_blue(rand() % rgb_limit);

			set_pixel32(image, i, j, p);
		}
	}

	if (SDL_MUSTLOCK(image)) {
		SDL_UnlockSurface(image);
	}

	return image;
}

// w_scale is the number of pixels per unit on the real axis, h_scale is the number of pixels per unit on the imaginary axis.
SDL_Surface * map_function_to_pixels(int width, int w_scale, int height, int h_scale, int bpp, ComplexNumber (*f)(ComplexNumber)) {
	// Create the surface to return, initially empty.
	int blue_mask = ~((-1) << bpp); // Red mask
	int green_mask = blue_mask << bpp; // Green mask
	int red_mask = blue_mask << (bpp * 2); // Blue mask

	SDL_Surface *coordinate_plane = SDL_CreateRGBSurface(SDL_SWSURFACE, // Flag specifying to create image in memory
													     width, // Width of image
													     height, // Height of image
													     bpp, // Bits per pixel
													     red_mask,
													     green_mask,
													     blue_mask,
													     0); // No alpha mask

	// Lock threads if necessary.
	if (SDL_MUSTLOCK(coordinate_plane)) {
		SDL_LockSurface(coordinate_plane);
	}

	// Compute the range of the given function over the given domain, and find the maximum (for use in pixel conversion).
	ComplexNumber **values = (ComplexNumber **) malloc(sizeof(ComplexNumber *) * width);
	ComplexNumber z (0, 0);
	double maximum = - numeric_limits<double>::infinity();
	double x, y;
	for (int i = 0; i < width; i++) {
		values[i] = (ComplexNumber *) malloc(sizeof(ComplexNumber) * height);
		for (int j = 0, y = - ((double) height) / (2 * h_scale); j < height; j++) {
			// Center the coordinate axes, scale the x- and y-coordinates appropriately.
			x = (i - ((double) width) / 2) / w_scale;
			y = (j - ((double) height) / 2) / h_scale;

			z.set(x, y); // Build the complex number
			values[i][j] = f(z);

			// Update maximum.
			if (values[i][j].mag() > maximum) {
				maximum = values[i][j].mag();
			}
		}
	}

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			set_pixel32(coordinate_plane, i, j, numberToPixel(values[i][j], maximum)); // Set the pixel to its value

			// Implement logic to handle poles (division by zero)
		}
	}

	// Unlock threads if necessary.
	if (SDL_MUSTLOCK(coordinate_plane)) {
		SDL_UnlockSurface(coordinate_plane);
	}

	// Free memory and return.
	for (int i = 0; i < width; i++) {
		free(values[i]);
	}
	free(values);

	return coordinate_plane;
}

// Displays the given surface for a given amount of time, in milliseconds.
void display_image(SDL_Surface *image, long display_time) {
	SDL_Surface *screen = SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, SDL_SWSURFACE); // Creates a 640 px by 480 px screen with 32 bits per pixel, in memory

	if (screen == NULL) {
		cerr << "Failed to initialize screen\n";
		exit(1);
	}

	if (SDL_BlitSurface(image, NULL, screen, NULL) != 0) {
		cerr << "failed: " << SDL_GetError() << endl;
		exit(1);
	}

	if (SDL_Flip(screen) == -1) {
		cerr << "Failed to update screen\n";
		exit(1);
	}
	SDL_Delay(display_time); // Display screen for 2000 milliseconds
}

// Converts HSL to RGB. Assumes 0 <= h, s, l <= 1.
int * hslToRgb(double h, double s, double l) {
	double r, g, b;

	if (s == 0) {
		r = g = b = l;
	} else {
		double q = (l < 0.5) ? (l * (1 + s)) : (l + s - l * s);
		double p = 2 * l - q;

		r = hslToRgbHelper(p, q, h + 1.0 / 3);
		g = hslToRgbHelper(p, q, h);
		b = hslToRgbHelper(p, q, h - 1.0 / 3);
	}

	int *ret = (int *) malloc(sizeof(int) * 3);
	ret[0] = (int) (r * 255 + 0.5);
	ret[1] = (int) (g * 255 + 0.5);
	ret[2] = (int) (b * 255 + 0.5);

	return ret;
}

double hslToRgbHelper(double p, double q, double t) {
	if (t < 0) {
		t++;
	} else if (t > 1) {
		t--;
	}

	if (t < 1.0 / 6) {
		return p + (q - p) * 6 * t;
	} else if (t < 1.0 / 2) {
		return q;
	} else if (t < 2.0 / 3) {
		return p + (q - p) * (2/3 - t) * 6;
	} else {
	    return p;		
	}
}