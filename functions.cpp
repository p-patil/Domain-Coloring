#include <iostream>
#include <cmath>
#include "functions.h"
#include "complex_numbers.h"

using namespace std;

// Begin polynomial class.

Polynomial::Polynomial(ComplexNumber coeffs[], int length) {
	degree = length - 1;
	coefficients = new ComplexNumber[degree + 1];
	for (int i = 0; i < degree + 1; i++) {
		coefficients[i] = coeffs[length - i - 1];
	}
}

int Polynomial::get_degree(void) const {
	return degree;
}

ComplexNumber * Polynomial::get_coefficients(void) const {
	return coefficients;
}

// Override
ComplexNumber Polynomial::eval(const ComplexNumber z) const {
	ComplexNumber result (0.0, 0.0);
	ComplexNumber pow_of_z (1, 0);
	for (int n = 0; n < degree + 1; n ++) {
		if (coefficients[n] != 0) {
			result += coefficients[n] * pow_of_z;
		}

		pow_of_z *= z;
	}

	return result;
}

// Destructor.
Polynomial::~Polynomial(void) {
	delete[] coefficients;
}

Polynomial Polynomial::operator =(const Polynomial &p) {
	if (this != &p) {
		degree = p.get_degree();

		delete[] coefficients;
		coefficients = new ComplexNumber[degree + 1];
		for (int i = 0; i < p.get_degree() + 1; i++) {
			coefficients[i] = p.get_coefficients()[i];
		}
	}

	return *this;
}

bool Polynomial::operator ==(const Polynomial &z) const {
	if (degree != z.get_degree()) {
		return false;
	}

	for (int i = 0; i < degree + 1; i++) {
		if (coefficients[i] != z.get_coefficients()[i]) {
			return false;
		}
	}

	return true;
}

bool Polynomial::operator !=(const Polynomial &z) const {
	return !(*this == z);
}

/**
 * Control how a polynomial is printed.
 */
ostream& operator <<(ostream &stream, const Polynomial &p) {
	ComplexNumber one (1, 0);
	ComplexNumber minus_one (-1, 0);
	ComplexNumber zero (0 ,0);

	if (p.get_coefficients()[p.get_degree()] != zero) {
		if (p.get_coefficients()[p.get_degree()] == one) {
			stream << "z^" << p.get_degree();
		} else if (p.get_coefficients()[p.get_degree()] == minus_one) {
			stream << "-z^" << p.get_degree();
		} else if (p.get_coefficients()[p.get_degree()].is_real()) {
			stream << p.get_coefficients()[p.get_degree()].real() << "z^" << p.get_degree();
		} else if (p.get_coefficients()[p.get_degree()].is_imaginary()) {
			stream << p.get_coefficients()[p.get_degree()].imag() << "iz^" << p.get_degree();
		} else {
			stream << "(" << p.get_coefficients()[p.get_degree()] << ")z^" << p.get_degree();
		}
	}

	for (int n = p.get_degree() - 1; n > 1; n--) {
		if (p.get_coefficients()[n] != zero) {	
			if (p.get_coefficients()[n] == one) {
				stream << " + z^" << n;
			} else if (p.get_coefficients()[n] == minus_one) {
				stream << " - z^" << n;
			} else if (p.get_coefficients()[n].is_real()) {
				if (p.get_coefficients()[n].real() > 0) {
					stream << " + " << p.get_coefficients()[n].real() << "z^" << n;
				} else {
					stream << " - " << (-p.get_coefficients()[n].real()) << "z^" << n;					
				}
			} else if (p.get_coefficients()[n].is_imaginary()) {
				if (p.get_coefficients()[n].imag() > 0) {
					stream << " + " << p.get_coefficients()[n].imag() << "iz^" << n;
				} else {
					stream << " - " << (-p.get_coefficients()[n].imag()) << "iz^" << n;
				}
			} else {
				stream << " + (" << p.get_coefficients()[n] << ")z^" << n;
			}
		}
	}

	// The cases n = 0, 1 require special formatting.
	if (p.get_coefficients()[1] != zero) {
		if (p.get_coefficients()[1] == one) {
			stream << " + z";
		} else if (p.get_coefficients()[1] == minus_one) {
			stream << "- z";
		} else if (p.get_coefficients()[1].is_real()) {
			if (p.get_coefficients()[1].real() > 0) {
				stream << " + " << p.get_coefficients()[1].real() << "z";
			} else {
				stream << " - " << (-p.get_coefficients()[1].real()) << "z";				
			}
		} else if (p.get_coefficients()[1].is_imaginary()) {
			if (p.get_coefficients()[1].imag() > 0) {
				stream << " + " << p.get_coefficients()[1].imag() << "iz";
			} else {
				stream << " - " << (-p.get_coefficients()[1].imag()) << "iz";				
			}
		} else {
			stream << " + (" << p.get_coefficients()[1] << ")z";
		}
	}

	if (p.get_coefficients()[0] != zero) {
		if (p.get_coefficients()[0].is_real()) {
			if (p.get_coefficients()[0].real() > 0) {
				stream << " + " << p.get_coefficients()[0].real();
			} else {
				stream << " - " << (-p.get_coefficients()[0].real());				
			}
		} else if (p.get_coefficients()[0].is_imaginary()) {
			if (p.get_coefficients()[0].imag() > 0) {
				stream << " + " << p.get_coefficients()[0].imag() << "i";
			} else {
				stream << " - " << (-p.get_coefficients()[0].imag()) << "i";				
			}
		} else {
			stream << " + (" << p.get_coefficients()[0] << ")";
		}
	}


	return stream;
}

// End polynomial class.

// Begin power class.

Power::Power(ComplexNumber exponent) {
	this->exponent = exponent;
}

ComplexNumber Power::get_exponent(void) {
	return exponent;
}

// Override.
// (r * exp(i * theta))^(a + bi) = r^(a + bi) * exp(i * theta * a - theta * b)
//	= (r^a * exp(-theta * b)) * (r^b * exp(theta * a))^i
ComplexNumber Power::eval(const ComplexNumber z) const {
	double r = z.mag(), theta = z.angle();
	double a = exponent.real(), b = exponent.imag();

	double coeff = pow(r, a) * exp(-theta * b);
	ComplexNumber w = iPower(pow(r, b) * exp(theta * a));

	ComplexNumber v (coeff * w.real(), coeff * w.imag());
	return v;
}

// Returns a^i.
// a^i = exp(i * ln(a)) = cos(ln(a)) + i * sin(ln(a))
ComplexNumber Power::iPower(double a) const {
	double pow_re = cos(log(a)), pow_im = sin(log(a));

	ComplexNumber w (pow_re, pow_im);
	return w;
}

// End extended polynomial class.

// Begin N-th root class.

NthRoot::NthRoot(int n) {
	this->n = n;
}

int NthRoot::get_n(void) {
	return n;
}

// Override.
ComplexNumber NthRoot::eval(const ComplexNumber z) const {
	Power power (1 / n);
	return power.eval(z);
}

// End N-th root class.

// Begin exponential class.

Exponential::Exponential(ComplexNumber coeff, ComplexNumber l) {
	this->coefficient = coeff;
	this->lambda = l;
}

ComplexNumber Exponential::get_coefficient(void) const {
	return coefficient;
}

ComplexNumber Exponential::get_lambda(void) const {
	return lambda;
}

ComplexNumber Exponential::eval(const ComplexNumber z) const {
	ComplexNumber exponent = lambda * z;

	double exp_re = exp(exponent.real()) * cos(exponent.imag());
	double exp_im = exp(exponent.real()) * sin(exponent.imag());

	ComplexNumber w (exp_re, exp_im);
	return coefficient * w;
}

Exponential Exponential::operator =(const Exponential &z) {
	if (this != &z) {
		coefficient = z.get_coefficient();
		lambda = z.get_lambda();
	}

	return *this;
}

bool Exponential::operator ==(const Exponential &z) const {
	return (coefficient == z.get_coefficient()) && (lambda == z.get_lambda());
}

bool Exponential::operator !=(const Exponential &z) const {
	return !(*this == z);
}

/**
 * Control how an exponential is printed.
 */
ostream& operator <<(ostream &stream, const Exponential &exp) {
	ComplexNumber one (1, 0);

	if (exp.get_coefficient() != one) {
		if (exp.get_coefficient().is_real() || exp.get_coefficient().is_imaginary()) {
			stream << exp.get_coefficient();			
		} else {
			stream << "(" << exp.get_coefficient() << ")";
		}
	}

	stream << "e^";

	if (exp.get_lambda() != one) {
		if (exp.get_coefficient().is_real() || exp.get_coefficient().is_imaginary()) {
			stream << "(" << exp.get_lambda() << "z)";
		} else {
			stream << "((" << exp.get_lambda() << ")z)";
		}
	} else {
		stream << "z";
	}

	return stream;
}

// End exponential class.

// Begin logarithm class.

// If no base is provided, Euler's number e is taken as the base by default.
Logarithm::Logarithm(ComplexNumber coefficient, ComplexNumber lambda) {
	this->coefficient = coefficient;
	this->lambda = lambda;
	base = exp(1);	
}

Logarithm::Logarithm(ComplexNumber coefficient, ComplexNumber lambda, ComplexNumber base) {
	this->coefficient = coefficient;
	this->lambda = lambda;
	this->base = base;
}

ComplexNumber Logarithm::get_coefficient(void) const {
	return coefficient;
}

ComplexNumber Logarithm::get_lambda(void) const {
	return lambda;
}

ComplexNumber Logarithm::get_base(void) const {
	return base;
}

// Override.
// ln(z) = ln(|r|) + i * theta, where z = r * e^(i * theta)
// log_z(w) = ln(w) * log_z(e) = ln(w) * 1 / (ln(r) + i * theta), where z = r * e^(i * theta).
ComplexNumber Logarithm::eval(const ComplexNumber z) const {
	ComplexNumber argument = z + lambda;
	ComplexNumber ln_z (log(argument.mag()), argument.angle());

	if (base == exp(1)) {
		return coefficient * ln_z;
	} else {
		ComplexNumber one (1, 0);
		ComplexNumber log_base (log(base.mag()), base.angle());
		return coefficient * ln_z * (one / log_base);
	}
}

Logarithm Logarithm::operator =(const Logarithm &z) {
	if (this != &z) {
		base = z.get_base();
		coefficient = z.get_coefficient();
		lambda = z.get_lambda();
	}

	return *this;
}

bool Logarithm::operator ==(const Logarithm &z) const {
	return (base == z.get_base()) && (coefficient == z.get_coefficient()) && (lambda == z.get_lambda());
}

bool Logarithm::operator !=(const Logarithm &z) const {
	return !(*this == z);
}

/**
 * Controls how a logarithm is printed.
 */
ostream& operator <<(ostream &stream, const Logarithm &logarithm) {
	ComplexNumber one (1, 0);
	ComplexNumber minus_one (-1, 0);

	if (logarithm.get_coefficient().is_real() || logarithm.get_coefficient().is_imaginary()) {
		if (logarithm.get_coefficient() == minus_one) {
			stream << "-";
		} else if (logarithm.get_coefficient() != one) {
			stream << logarithm.get_coefficient() << " * ";
		}
	} else {
		stream << "(" << logarithm.get_coefficient() << ") * ";
	}
	
	if (logarithm.get_base() == exp(1)) {
		stream << "ln(";
	} else {
		stream << "log_(" << logarithm.get_base() << ") (";
	}

	stream << "z";

	if (logarithm.get_lambda().is_real()) {
		if (logarithm.get_lambda().real() > 0) {
			stream << " + " << logarithm.get_lambda().real() << ")";
		} else if (logarithm.get_lambda().real() < 0) {
			stream << " - " << (-logarithm.get_lambda().real()) << ")";
		} else {
			stream << ")";
		}
	} else if (logarithm.get_lambda().is_imaginary()) {
		if (logarithm.get_lambda().imag() > 0) {
			if (logarithm.get_lambda().imag() == 1) {
				stream << " + i)";
			} else {
				stream << " + " << logarithm.get_lambda().imag() << "i)";
			}
		} else if (logarithm.get_lambda().imag() < 0) {
			if (logarithm.get_lambda().imag() == -1) {
				stream << " - i)";
			} else {
				stream << " - " << (-logarithm.get_lambda().imag()) << "i)";			
			}
		} else {
			stream << ")";
		}
	} else {
		if (logarithm.get_lambda().real() > 0) {
			stream << " + " << logarithm.get_lambda() << ")";
		} else if (logarithm.get_lambda().real() < 0) {
			ComplexNumber new_lambda (-logarithm.get_lambda().real(), logarithm.get_lambda().imag());
			stream << " - " << new_lambda << ")";
		}
	}

	return stream;
}

// End logarithm class.

// Begin sine class.

Sine::Sine() {}

Sine::Sine(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Sine::getA(void) const {
	return a;
}

ComplexNumber Sine::getB(void) const {
	return b;
}

ComplexNumber Sine::getC(void) const {
	return c;
}

ComplexNumber Sine::eval(const ComplexNumber z) const {
	ComplexNumber i (0, 1);
	ComplexNumber arg = b * z + c;

	Exponential e (1, 1);
	return a * (e.eval(i * arg) - e.eval(-i * arg)) / (2 * i);
}

/**
 * Control how the sine function is printed.
 */
/**
 * Control how the sine function is printed.
 */
ostream& operator <<(ostream &stream, const Sine &s) {
	if (s.getA() == 0 || (s.getB() == 0 && s.getC() == 0)) {
		return stream << "0";
	}

	if (s.getA() != 1) {
		if (s.getA().real() == 0 || s.getA().imag() == 0) {
			stream << s.getA();
		} else {
			stream << "(" << s.getA() << ")";
		}

		stream << " ";
	}

	stream << "sin(";

	if (s.getB() != 1) {
		if (s.getB().real() == 0 || s.getB().imag() == 0) {
			stream << s.getB();
		} else {
			stream << "(" << s.getB() << ")";
		}
	}

	stream << "z";

	if (s.getC() != 0) {
		if (s.getC().is_real()) {
			if (s.getC().real() < 0) {
				stream << " - " << -s.getC().real();
			} else {
				stream << " + " << s.getC().real();
			}
		} else if (s.getC().is_imaginary()) {
			if (s.getC().imag() < 0) {
				stream << " - " << -s.getC().imag();
			} else {
				stream << " + " << s.getC().imag();
			}
		} else {
			stream << s.getC();
		}
	}

	return stream << ")";
}

// End sine class.

// Begin cosine class.

Cosine::Cosine() {}

Cosine::Cosine(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Cosine::getA(void) const {
	return a;
}

ComplexNumber Cosine::getB(void) const {
	return b;
}

ComplexNumber Cosine::getC(void) const {
	return c;
}

ComplexNumber Cosine::eval(const ComplexNumber z) const {
	ComplexNumber i (0, 1);
	ComplexNumber arg = b * z + c;

	Exponential e (1, 1);
	return a * (e.eval(i * arg) + e.eval(-i * arg)) / 2;
}

/**
 * Control how the cosine function is printed.
 */
ostream& operator <<(ostream &stream, const Cosine &s) {
	if (s.getA() == 0 || (s.getB() == 0 && s.getC() == 0)) {
		return stream << "0";
	}

	if (s.getA() != 1) {
		if (s.getA().real() == 0 || s.getA().imag() == 0) {
			stream << s.getA();
		} else {
			stream << "(" << s.getA() << ")";
		}

		stream << " ";
	}

	stream << "cos(";

	if (s.getB() != 1) {
		if (s.getB().real() == 0 || s.getB().imag() == 0) {
			stream << s.getB();
		} else {
			stream << "(" << s.getB() << ")";
		}
	}

	stream << "z";

	if (s.getC() != 0) {
		if (s.getC().is_real()) {
			if (s.getC().real() < 0) {
				stream << " - " << -s.getC().real();
			} else {
				stream << " + " << s.getC().real();
			}
		} else if (s.getC().is_imaginary()) {
			if (s.getC().imag() < 0) {
				stream << " - " << -s.getC().imag();
			} else {
				stream << " + " << s.getC().imag();
			}
		} else {
			stream << s.getC();
		}
	}

	return stream << ")";
}

// End cosine class.

// Begin tangent class.

Tangent::Tangent(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Tangent::getA(void) const {
	return a;
}

ComplexNumber Tangent::getB(void) const {
	return b;
}

ComplexNumber Tangent::getC(void) const {
	return c;
}

ComplexNumber Tangent::eval(const ComplexNumber z) const {
	Sine sine (a, b, c);
	Cosine cosine (a, b, c);

	return sine.eval(z) / cosine.eval(z);
}

/**
 * Control how the tangent function is printed.
 */
ostream& operator <<(ostream &stream, const Tangent &s) {
	if (s.getA() == 0 || (s.getB() == 0 && s.getC() == 0)) {
		return stream << "0";
	}

	if (s.getA() != 1) {
		if (s.getA().real() == 0 || s.getA().imag() == 0) {
			stream << s.getA();
		} else {
			stream << "(" << s.getA() << ")";
		}

		stream << " ";
	}

	stream << "tan(";

	if (s.getB() != 1) {
		if (s.getB().real() == 0 || s.getB().imag() == 0) {
			stream << s.getB();
		} else {
			stream << "(" << s.getB() << ")";
		}
	}

	stream << "z";

	if (s.getC() != 0) {
		if (s.getC().is_real()) {
			if (s.getC().real() < 0) {
				stream << " - " << -s.getC().real();
			} else {
				stream << " + " << s.getC().real();
			}
		} else if (s.getC().is_imaginary()) {
			if (s.getC().imag() < 0) {
				stream << " - " << -s.getC().imag();
			} else {
				stream << " + " << s.getC().imag();
			}
		} else {
			stream << s.getC();
		}
	}

	return stream << ")";
}

// End tangent class.