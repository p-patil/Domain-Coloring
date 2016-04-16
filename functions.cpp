#include <cmath>
#include "functions.h"
#include "complex_numbers.h"

using namespace std;

void Function::print_to_stream(ostream &stream) const {
	stream << "MEOW";
}

// Base class << operator to allow for printing.
ostream& operator <<(ostream &stream, const Function &f) {
	f.print_to_stream(stream << "f(z) = ");
	return stream;
}

// Begin constant class.

Constant::Constant(const ComplexNumber z) { this->constant = z; }

ComplexNumber Constant::get_constant(void) const { return this->constant; }

ComplexNumber Constant::eval(const ComplexNumber z) const { return this->constant; }

void Constant::print_to_stream(ostream &stream) const { stream << this->constant; }

// End constant class.

// Begin polynomial class.

Polynomial::Polynomial(ComplexNumber coeffs[], int length) {
	degree = length - 1;
	coefficients = new ComplexNumber[degree + 1];
	for (int i = 0; i < degree + 1; i++) {
		coefficients[i] = coeffs[length - i - 1];
	}
}

int Polynomial::get_degree(void) const { return degree; }

ComplexNumber * Polynomial::get_coefficients(void) const { return coefficients; }

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
Polynomial::~Polynomial(void) { delete[] coefficients; }

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

bool Polynomial::operator !=(const Polynomial &z) const { return !(*this == z); }

void Polynomial::print_to_stream(ostream &stream) const {
	ComplexNumber one (1, 0);
	ComplexNumber minus_one (-1, 0);
	ComplexNumber zero (0 ,0);

	if (this->coefficients[this->degree] != zero) {
		if (this->coefficients[this->degree] == one) {
			stream << "z^" << this->degree;
		} else if (this->coefficients[this->degree] == minus_one) {
			stream << "-z^" << this->degree;
		} else if (this->coefficients[this->degree].is_real()) {
			stream << this->coefficients[this->degree].real() << "z^" << this->degree;
		} else if (this->coefficients[this->degree].is_imaginary()) {
			stream << this->coefficients[this->degree].imag() << "iz^" << this->degree;
		} else {
			stream << "(" << this->coefficients[this->degree] << ")z^" << this->degree;
		}
	}

	for (int n = this->degree - 1; n > 1; n--) {
		if (this->coefficients[n] != zero) {	
			if (this->coefficients[n] == one) {
				stream << " + z^" << n;
			} else if (this->coefficients[n] == minus_one) {
				stream << " - z^" << n;
			} else if (this->coefficients[n].is_real()) {
				if (this->coefficients[n].real() > 0) {
					stream << " + " << this->coefficients[n].real() << "z^" << n;
				} else {
					stream << " - " << (-this->coefficients[n].real()) << "z^" << n;					
				}
			} else if (this->coefficients[n].is_imaginary()) {
				if (this->coefficients[n].imag() > 0) {
					stream << " + " << this->coefficients[n].imag() << "iz^" << n;
				} else {
					stream << " - " << (-this->coefficients[n].imag()) << "iz^" << n;
				}
			} else {
				stream << " + (" << this->coefficients[n] << ")z^" << n;
			}
		}
	}

	// The cases n = 0, 1 require special formatting.
	if (this->coefficients[1] != zero) {
		if (this->coefficients[1] == one) {
			stream << " + z";
		} else if (this->coefficients[1] == minus_one) {
			stream << "- z";
		} else if (this->coefficients[1].is_real()) {
			if (this->coefficients[1].real() > 0) {
				stream << " + " << this->coefficients[1].real() << "z";
			} else {
				stream << " - " << (-this->coefficients[1].real()) << "z";				
			}
		} else if (this->coefficients[1].is_imaginary()) {
			if (this->coefficients[1].imag() > 0) {
				stream << " + " << this->coefficients[1].imag() << "iz";
			} else {
				stream << " - " << (-this->coefficients[1].imag()) << "iz";				
			}
		} else {
			stream << " + (" << this->coefficients[1] << ")z";
		}
	}

	if (this->coefficients[0] != zero) {
		if (this->coefficients[0].is_real()) {
			if (this->coefficients[0].real() > 0) {
				stream << " + " << this->coefficients[0].real();
			} else {
				stream << " - " << (-this->coefficients[0].real());				
			}
		} else if (this->coefficients[0].is_imaginary()) {
			if (this->coefficients[0].imag() > 0) {
				stream << " + " << this->coefficients[0].imag() << "i";
			} else {
				stream << " - " << (-this->coefficients[0].imag()) << "i";				
			}
		} else {
			stream << " + (" << this->coefficients[0] << ")";
		}
	}
}

// End polynomial class.

// Begin power class.

Power::Power(ComplexNumber exponent) { this->exponent = exponent; }

ComplexNumber Power::get_exponent(void) const { return exponent; }

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

void Power::print_to_stream(ostream &stream) const {
	ComplexNumber zero (0, 0);
	ComplexNumber one (1, 0);
	
	if (this->exponent == zero) {
		stream << 1;
	} else if (this->exponent == one) {
		stream << "z";
	} else {
		stream << "z^";
		if (this->exponent.is_real()) {
			if (this->exponent.real() < 0) {
				stream << "(" << this->exponent << ")";
			} else {
				stream << this->exponent;
			}
		} else {
			stream << "(" << this->exponent << ")";
		}
	}
}

// End power class.

// Begin n-th root class.

NthRoot::NthRoot(int n) { this->n = n; }

int NthRoot::get_n(void) const { return n; }

// Override.
ComplexNumber NthRoot::eval(const ComplexNumber z) const {
	Power power (1 / n);
	return power.eval(z);
}

void NthRoot::print_to_stream(ostream &stream) const {
	stream << "z";
	if (this->n != 1) {
		stream << "^(1 / ";

		if (this->n < 0) {
			stream << "(" << this->n << ")";
		}

		stream << ")";
	}
}

// End n-th root class.

// Begin exponential class.

Exponential::Exponential(ComplexNumber coeff, ComplexNumber l) {
	this->coefficient = coeff;
	this->lambda = l;
}

ComplexNumber Exponential::get_coefficient(void) const { return coefficient; }

ComplexNumber Exponential::get_lambda(void) const { return lambda; }

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

bool Exponential::operator !=(const Exponential &z) const { return !(*this == z); }

void Exponential::print_to_stream(ostream &stream) const {
	ComplexNumber one (1, 0);

	if (this->coefficient != one) {
		if (this->coefficient.is_real() || this->coefficient.is_imaginary()) {
			stream << this->coefficient;			
		} else {
			stream << "(" << this->coefficient << ")";
		}
	}

	stream << "e^";

	if (this->lambda != one) {
		if (this->coefficient.is_real() || this->coefficient.is_imaginary()) {
			stream << "(" << this->lambda << "z)";
		} else {
			stream << "((" << this->lambda << ")z)";
		}
	} else {
		stream << "z";
	}
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

ComplexNumber Logarithm::get_coefficient(void) const { return coefficient; }

ComplexNumber Logarithm::get_lambda(void) const { return lambda; }

ComplexNumber Logarithm::get_base(void) const { return base; }

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

bool Logarithm::operator !=(const Logarithm &z) const { return !(*this == z); }

void Logarithm::print_to_stream(ostream &stream) const {
	ComplexNumber one (1, 0);
	ComplexNumber minus_one (-1, 0);

	if (this->coefficient.is_real() || this->coefficient.is_imaginary()) {
		if (this->coefficient == minus_one) {
			stream << "-";
		} else if (this->coefficient != one) {
			stream << this->coefficient << " * ";
		}
	} else {
		stream << "(" << this->coefficient << ") * ";
	}
	
	if (this->base == exp(1)) {
		stream << "ln(";
	} else {
		stream << "log_(" << this->base << ") (";
	}

	stream << "z";

	if (this->lambda.is_real()) {
		if (this->lambda.real() > 0) {
			stream << " + " << this->lambda.real() << ")";
		} else if (this->lambda.real() < 0) {
			stream << " - " << (-this->lambda.real()) << ")";
		} else {
			stream << ")";
		}
	} else if (this->lambda.is_imaginary()) {
		if (this->lambda.imag() > 0) {
			if (this->lambda.imag() == 1) {
				stream << " + i)";
			} else {
				stream << " + " << this->lambda.imag() << "i)";
			}
		} else if (this->lambda.imag() < 0) {
			if (this->lambda.imag() == -1) {
				stream << " - i)";
			} else {
				stream << " - " << (-this->lambda.imag()) << "i)";			
			}
		} else {
			stream << ")";
		}
	} else {
		if (this->lambda.real() > 0) {
			stream << " + " << this->lambda << ")";
		} else if (this->lambda.real() < 0) {
			ComplexNumber new_lambda (-this->lambda.real(), this->lambda.imag());
			stream << " - " << new_lambda << ")";
		}
	}
}

// End logarithm class.

// Begin sine class.

Sine::Sine() {}

Sine::Sine(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Sine::getA(void) const { return a; }

ComplexNumber Sine::getB(void) const { return b; }

ComplexNumber Sine::getC(void) const { return c; }

ComplexNumber Sine::eval(const ComplexNumber z) const {
	ComplexNumber i (0, 1);
	ComplexNumber arg = b * z + c;

	Exponential e (1, 1);
	return a * (e.eval(i * arg) - e.eval(-i * arg)) / (2 * i);
}

void Sine::print_to_stream(ostream &stream) const {
	if (this->a == 0 || (this->b == 0 && this->c == 0)) {
		stream << "0";
	}

	if (this->a != 1) {
		if (this->a.real() == 0 || this->a.imag() == 0) {
			stream << this->a;
		} else {
			stream << "(" << this->a << ")";
		}

		stream << " ";
	}

	stream << "sin(";

	if (this->b != 1) {
		if (this->b.real() == 0 || this->b.imag() == 0) {
			stream << this->b;
		} else {
			stream << "(" << this->b << ")";
		}
	}

	stream << "z";

	if (this->c != 0) {
		if (this->c.is_real()) {
			if (this->c.real() < 0) {
				stream << " - " << -this->c.real();
			} else {
				stream << " + " << this->c.real();
			}
		} else if (this->c.is_imaginary()) {
			if (this->c.imag() < 0) {
				stream << " - " << -this->c.imag();
			} else {
				stream << " + " << this->c.imag();
			}
		} else {
			stream << this->c;
		}
	}

	stream << ")";
}

// End sine class.

// Begin cosine class.

Cosine::Cosine() {}

Cosine::Cosine(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Cosine::getA(void) const { return a; }

ComplexNumber Cosine::getB(void) const { return b; }

ComplexNumber Cosine::getC(void) const { return c; }

ComplexNumber Cosine::eval(const ComplexNumber z) const {
	ComplexNumber i (0, 1);
	ComplexNumber arg = b * z + c;

	Exponential e (1, 1);
	return a * (e.eval(i * arg) + e.eval(-i * arg)) / 2;
}

void Cosine::print_to_stream(ostream &stream) const {
	if (this->a == 0 || (this->b == 0 && this->c == 0)) {
		stream << "0";
	}

	if (this->a != 1) {
		if (this->a.real() == 0 || this->a.imag() == 0) {
			stream << this->a;
		} else {
			stream << "(" << this->a << ")";
		}

		stream << " ";
	}

	stream << "cos(";

	if (this->b != 1) {
		if (this->b.real() == 0 || this->b.imag() == 0) {
			stream << this->b;
		} else {
			stream << "(" << this->b << ")";
		}
	}

	stream << "z";

	if (this->c != 0) {
		if (this->c.is_real()) {
			if (this->c.real() < 0) {
				stream << " - " << -this->c.real();
			} else {
				stream << " + " << this->c.real();
			}
		} else if (this->c.is_imaginary()) {
			if (this->c.imag() < 0) {
				stream << " - " << -this->c.imag();
			} else {
				stream << " + " << this->c.imag();
			}
		} else {
			stream << this->c;
		}
	}

	stream << ")";
}

// // End cosine class.

// Begin tangent class.

Tangent::Tangent(ComplexNumber a, ComplexNumber b, ComplexNumber c) {
	this->a = a;
	this->b = b;
	this->c = c;
}

ComplexNumber Tangent::getA(void) const { return a; }

ComplexNumber Tangent::getB(void) const { return b; }

ComplexNumber Tangent::getC(void) const { return c; }

ComplexNumber Tangent::eval(const ComplexNumber z) const {
	Sine sine (a, b, c);
	Cosine cosine (a, b, c);

	return sine.eval(z) / cosine.eval(z);
}

void Tangent::print_to_stream(ostream &stream) const {
	if (this->a == 0 || (this->b == 0 && this->c == 0)) {
		stream << "0";
	}

	if (this->a != 1) {
		if (this->a.real() == 0 || this->a.imag() == 0) {
			stream << this->a;
		} else {
			stream << "(" << this->a << ")";
		}

		stream << " ";
	}

	stream << "tan(";

	if (this->b != 1) {
		if (this->b.real() == 0 || this->b.imag() == 0) {
			stream << this->b;
		} else {
			stream << "(" << this->b << ")";
		}
	}

	stream << "z";

	if (this->c != 0) {
		if (this->c.is_real()) {
			if (this->c.real() < 0) {
				stream << " - " << -this->c.real();
			} else {
				stream << " + " << this->c.real();
			}
		} else if (this->c.is_imaginary()) {
			if (this->c.imag() < 0) {
				stream << " - " << -this->c.imag();
			} else {
				stream << " + " << this->c.imag();
			}
		} else {
			stream << this->c;
		}
	}

	stream << ")";
}

// End tangent class.