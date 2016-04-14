#include <cmath>
#include <stdexcept>
#include <iostream>
#include "complex_numbers.h"

using namespace std;

ComplexNumber::ComplexNumber() {};

ComplexNumber::ComplexNumber(double real, double imag) : re(real), im(imag) { initialize_polar(real, imag); }

ComplexNumber::ComplexNumber(double x) : re(x), im(0) {
	r = abs(x);

	if (x < 0) {
		theta = PI;
	} else {
		theta = 0;
	}
}

void ComplexNumber::initialize_polar(double real, double imag) {
	r = sqrt(real * real + imag * imag);

	if (real == 0 && imag == 0) {
		theta = 0;
	} else if (imag == 0) {
		if (real > 0) {
			theta = 0;
		} else {
			theta = PI;
		}
	} else if (real == 0) {
		if (imag > 0) {
			theta = PI / 2;
		} else {
			theta = - PI / 2;
		}
	} else {
		theta = atan(abs(imag) / abs(real));

		if (real < 0 && imag > 0) {
			theta = PI - theta;
		} else if (real < 0 && imag < 0) {
			theta = -(PI - theta);
		} else if (real > 0 && imag < 0) {
			theta = -theta;
		}
	}
}

/**
 * Returns real part this complex number.
 */
double ComplexNumber::real(void) const { return re; }

/**
 * Returns the imaginary part of this complex number.
 */
double ComplexNumber::imag(void) const { return im; }

/**
 * Magnitude of the complex number.
 */
double ComplexNumber::mag(void) const { return r; }

/**
 * Polar angle of the complex number, in radians.
 */
double ComplexNumber::angle(void) const { return theta; }

/**
 * Complex conjugate.
 */
ComplexNumber ComplexNumber::conjugate(void) const {
	ComplexNumber w (re, -im);
	return w;
}

bool ComplexNumber::is_real(void) { return im == 0; }

bool ComplexNumber::is_imaginary(void) { return re == 0 && im != 0; }

void ComplexNumber::set(double real, double imag) {
	re = real;
	im = imag;

	initialize_polar(real, imag);
}

/**
 * Assignment operator.
 */
ComplexNumber ComplexNumber::operator =(const ComplexNumber &z) {
	if (this != &z) { // Guard against self-assignment.
		re = z.real();
		im = z.imag();
		r = z.mag();
		theta = z.angle();
	}

	return *this;
}

bool ComplexNumber::operator ==(const ComplexNumber &z) const {
	if (re == z.real() && im == z.imag()) {
		return true;
	} else {
		return false;
	}
}

bool ComplexNumber::operator !=(const ComplexNumber &z) const { return !(*this == z); }

/**
 * Complex addition.
 */
ComplexNumber ComplexNumber::operator +(const ComplexNumber &z) const {
	double sum_re = re + z.real();
	double sum_im = im + z.imag();

	ComplexNumber w (sum_re, sum_im);
	return w;
}

ComplexNumber operator +(double x, const ComplexNumber &z) {
	ComplexNumber x_promoted (x, 0);
	return x_promoted + z;
}

/** 
 * Comlex subtraction.
 */
ComplexNumber ComplexNumber::operator -(const ComplexNumber &z) const {
	double sub_re = re - z.real();
	double sub_im = im - z.imag();

	ComplexNumber w (sub_re, sub_im);
	return w;
}

ComplexNumber operator -(double x, const ComplexNumber &z) {
	ComplexNumber x_promoted (x, 0);
	return x_promoted - z;
}

/**
 * Complex multiplication.
 */
ComplexNumber ComplexNumber::operator *(const ComplexNumber &z) const {
	double prod_re = re * z.real() - im * z.imag();
	double prod_im = re * z.imag() + im * z.real();

	ComplexNumber w (prod_re, prod_im);
	return w;
}

ComplexNumber operator *(double x, const ComplexNumber &z) {
	ComplexNumber x_promoted (x, 0);
	return x_promoted * z;
}

/**
 * Complex division.
 */
ComplexNumber ComplexNumber::operator /(const ComplexNumber &z) const {
	if (z.mag() == 0) {
		throw logic_error("Division by zero");
	}

	double div_re = (re * z.real() + im * z.imag()) / (z.mag() * z.mag());
	double div_im = (im * z.real() - re * z.imag()) / (z.mag() * z.mag());

	ComplexNumber w (div_re, div_im);
	return w;
}

ComplexNumber operator /(double x, const ComplexNumber &z) {
	ComplexNumber x_promoted (x, 0);
	return x_promoted / z;
}

ComplexNumber ComplexNumber::operator -() const {
	ComplexNumber zero (0, 0);
	return zero - *this;
}

ComplexNumber ComplexNumber::operator +=(const ComplexNumber &z) {
	*this = *this + z;
	return *this;
}

ComplexNumber ComplexNumber::operator -=(const ComplexNumber &z) {
	*this = *this - z;
	return *this;
}

ComplexNumber ComplexNumber::operator *=(const ComplexNumber &z) {
	*this = *this * z;
	return *this;
}

ComplexNumber ComplexNumber::operator /=(const ComplexNumber &z) {
	*this = *this / z;
	return *this;
}


/**
 * Control how a complex number is printed.
 */
ostream& operator <<(ostream &stream, const ComplexNumber &z) {
	if (z.real() == 0 && z.imag() == 0) {
		stream << "0";
	} else if (z.real() == 0) {
		if (z.imag() == 1) {
			stream << "i";
		} else {
			stream << z.imag() << "i";
		}
	} else if (z.imag() == 0) {
		stream << z.real();
	} else {
		if (z.imag() == 1) {
			stream << z.real() << " + i";
		} else if (z.imag() == -1) {
			stream << z.real() << " - i";
		} else if (z.imag() > 0) {
			stream << z.real() << " + " << z.imag() << "i";
		} else {
			stream << z.real() << " - " << (-z.imag()) << "i";					
		}
	}

	return stream;
}