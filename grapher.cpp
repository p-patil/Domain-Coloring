#include <SDL/SDL.h>
#include <iostream>
#include <limits>
#include <functional>

#include "grapher.h"
#include "complex_numbers.h"
#include "omp.h"

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

int RGBPixel::get_green(void) const {    return this->green; }

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

// Override
ostream& operator <<(ostream &stream, const RGBPixel &p) {
    stream << "(" << p.get_red() << ", " << p.get_green() << ", " << p.get_blue() << ")";
    return stream;
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

// TODO: Make this function compatible with GPUs using CUDA (dramatically speeds up computation)
// Primary function that returns a surface which displays the graph of the given complex function in the given bounds. w_scale is the number of pixels 
// per unit on the real axis, h_scale is the number of pixels per unit on the imaginary axis.
SDL_Surface * map_function_to_pixels(int width, int height, int bpp, double x_min, double x_max, double y_min, double y_max, function<ComplexNumber(ComplexNumber)> f) {
    // Create the surface to return, initially empty.
    int blue_mask = ~((-1) << bpp); // Blue mask
    int green_mask = blue_mask << bpp; // Green mask
    int red_mask = blue_mask << (bpp * 2); // Red mask

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
    double *maximum_arr = (double *) malloc(sizeof(double) * omp_get_num_threads());
    for (int i = 0; i < omp_get_num_threads(); i++) {
        maximum_arr[i] = - numeric_limits<double>::infinity();
    }

    ComplexNumber z (0, 0);
    double x, y;
    #pragma omp parallel for private(z, x, y) // Parallelize computation.
    for (int i = 0; i < width; i++) {
        values[i] = (ComplexNumber *) malloc(sizeof(ComplexNumber) * height);

        for (int j = 0; j < height; j++) {
            // Center the coordinate axes, scale the x- and y-coordinates appropriately.
            x = x_min + (((double) i) / width) * (x_max - x_min);
            y = y_min + (((double) j) / height) * (y_max - y_min);

            z.set(x, y); // Build the complex number

            try { // Store the image of f at z
                values[i][j] = f(z);
            } catch (ComplexDivisionByZeroException &e) {
                // Caught division by zero exception, so z is a pole of f; set to infinity
                z.set(numeric_limits<double>::infinity(), numeric_limits<double>::infinity());
                values[i][j] = z;
            }

            // Update maximum.
            if (values[i][j].mag() > maximum_arr[omp_get_thread_num()] && z.mag() != numeric_limits<double>::infinity()) {
                maximum_arr[omp_get_thread_num()] = values[i][j].mag();
            }
        }
    }

    // Join all maxima found by threads
    double maximum = - numeric_limits<double>::infinity();
    for (int i = 0; i < omp_get_num_threads(); i++) {
        if (maximum_arr[i] > maximum) {
            maximum = maximum_arr[i];
        }
    }
    free(maximum_arr);

    z.set(numeric_limits<double>::infinity(), numeric_limits<double>::infinity()); // Represents infinity
    RGBPixel white (8, 255, 255, 255); // White pixel, for coloring infinite values

    #pragma omp parallel for
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (values[i][j] != z) {
                set_pixel32(coordinate_plane, i, j, numberToPixel(values[i][j], maximum)); // Set the pixel to its value
            } else {
                set_pixel32(coordinate_plane, i, j, white);
            }
        }

        free(values[i]); // Free memory after use
    }

    // Unlock threads if necessary.
    if (SDL_MUSTLOCK(coordinate_plane)) {
        SDL_UnlockSurface(coordinate_plane);
    }

    // Free memory and return.
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

// HELPER FUNCTIONS

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

RGBPixel numberToPixel(ComplexNumber z, double maximum) {
    RGBPixel p (8, 0, 0, 0);
    double hue, saturation, lightness;

    hue = (z.angle() + PI) / (2 * PI); // Add PI to transform range from [- PI, PI] to [0, 2 * PI]
    saturation = 1;
    lightness = z.mag() / maximum;

    double *rgb = hsl_to_rgb(hue, saturation, lightness);

    // Add rounded versions of RGB to pixel, then return.
    p.set_red((int) (rgb[0] + 0.5));
    p.set_green((int) (rgb[1] + 0.5));
    p.set_blue((int) (rgb[2] + 0.5));

    free(rgb);
    return p;
}

// Converts HSL to RGB. Assumes 0 <= h, s, l <= 1.
double * hsl_to_rgb(double h, double s, double l) {
    double mid1, mid2, v, m;
    double r, g, b;

    r = g = b = l;
    v = (l <= 0.5) ? (l * (1.0 + s)) : (l + s - l * s);

    if (v > 0) {
        int sextant = (int) (h * 6);
        m = l + l - v;
        mid1 = m + ((v - m) * ((h * 6) - sextant));
        mid2 = v - ((v - m) * ((h * 6) - sextant));

        switch (sextant) {
            case 0:
                r = v;
                g = mid1;
                b = m;
                break;
            case 1:
                r = mid2;
                g = v;
                b = m;
                break;
            case 2:
                r = m;
                g = v;
                b = mid1;
                break;
            case 3:
                r = m;
                g = mid2;
                b = v;
                break;
            case 4:
                r = mid1;
                g = m;
                b = v;
                break;
            case 5:
                r = v;
                g = m;
                b = mid2;
                break;
        }
    }
    
    double * rgb = (double *) malloc(sizeof(double) * 3);
    rgb[0] = r * 255;
    rgb[1] = g * 255;
    rgb[2] = b * 255;

    return rgb;        
}