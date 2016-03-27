#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include "complex_numbers.h"

using namespace std;

class Function {
	public:
		virtual ComplexNumber eval(const ComplexNumber) const = 0;
};

// Polynomials.
class Polynomial : public Function {
	int degree;
	ComplexNumber *coefficients; // Index n corresponds to coefficient of x^n

	public:
		Polynomial(ComplexNumber[], int);

		int get_degree(void) const;

		ComplexNumber * get_coefficients(void) const;

		// Override
		ComplexNumber eval(const ComplexNumber) const;

		// Destructor.
		~Polynomial(void);

		Polynomial operator =(const Polynomial &);

		bool operator ==(const Polynomial &) const;

		bool operator !=(const Polynomial &) const;
};

/**
 * Control how a polynomial is printed.
 */
ostream& operator <<(ostream &, const Polynomial &);

// Class representing functions of the form f(z) = z^a, for complex number a.
class Power : public Function {
	private:
		ComplexNumber iPower(double) const;

	public:
		ComplexNumber exponent;

		Power(ComplexNumber);

		ComplexNumber get_exponent(void);

		// Override.
		ComplexNumber eval(const ComplexNumber) const;
};

class NthRoot : public Function {
	public:
		int n;

		NthRoot(int);

		int get_n(void);

		// Override.
		ComplexNumber eval(const ComplexNumber) const;
};

// Exponential function and its inverse, the logarithm.
class Exponential : public Function {
	// Models the function f(z) = coefficient * e^(lambda * z)
	ComplexNumber coefficient, lambda;

	public:
		Exponential(ComplexNumber, ComplexNumber);

		ComplexNumber get_coefficient(void) const;

		ComplexNumber get_lambda(void) const ;

		// Override.
		ComplexNumber eval(const ComplexNumber) const;

		Exponential operator =(const Exponential &);

		bool operator ==(const Exponential &) const;

		bool operator !=(const Exponential &) const;
};

/**
 * Control how an exponential is printed.
 */
ostream& operator <<(ostream &, const Exponential &);

// Logarithms. Because the class represents the complex logarithm, which is multivalued, the principal log is taken.
class Logarithm : public Function {
	ComplexNumber base, coefficient, lambda; // Represents coefficient * log_(base) (z + lambda)

	public:

		Logarithm(ComplexNumber, ComplexNumber);
		
		Logarithm(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber get_coefficient(void) const;

		ComplexNumber get_lambda(void) const;

		ComplexNumber get_base(void) const;

		// Override.
		ComplexNumber eval(const ComplexNumber) const;

		Logarithm operator =(const Logarithm &);

		bool operator ==(const Logarithm &) const;

		bool operator !=(const Logarithm &) const;
};

/**
 * Control how a logarithm is printed.
 */
ostream& operator <<(ostream &stream, const Logarithm &);


// Trig functions and their inverses.
class Sine : public Function {
	ComplexNumber a, b, c; // Represents a * sin(b * z + c)

	public:
		Sine();

		Sine(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;

		ComplexNumber getB(void) const;

		ComplexNumber getC(void) const;

		// Override
		ComplexNumber eval(const ComplexNumber) const;
};

/**
 * Control how the sine function is printed.
 */
ostream& operator <<(ostream &stream, const Sine &);

class Cosine : public Function {
	ComplexNumber a, b, c; // Represents a * cos(b * z + c)

	public:
		Cosine();

		Cosine(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;

		ComplexNumber getB(void) const;

		ComplexNumber getC(void) const;

		ComplexNumber eval(const ComplexNumber) const;
};

/**
 * Control how the cosine function is printed.
 */
ostream& operator <<(ostream &stream, const Cosine &);

class Tangent : public Function {
	ComplexNumber a, b, c; // Represents a * tan(b * z + c)

	public:
		Tangent(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;
		
		ComplexNumber getB(void) const;
		
		ComplexNumber getC(void) const;

		// Override.
		ComplexNumber eval(const ComplexNumber) const;
};

/**
 * Control how the tangent function is printed.
 */
ostream& operator <<(ostream &stream, const Tangent &);

class Cosecant : public Function {

};

class Secant : public Function {

};

class Cotangent : public Function {

};

class Arcsine : public Function {

};

class Arccosine : public Function {

};

class Arctangent : public Function {

};

class Arccosecant : public Function {

};

class Arcsecant : public Function {

};

class HyperbolicSine : public Function {

};

class HyperbolicCosine : public Function {

};

class HyperbolicTangent : public Function {

};

class HyperbolicCosecant : public Function {

};

class HyperbolicSecant : public Function {

};

class HyperbolicCotangent : public Function {

};

class InverseHyperbolicSine : public Function {

};

class InverseHyperbolicCosine : public Function {

};

class InverseHyperbolicTangent : public Function {

};

class InverseHyperbolicCosecant : public Function {

};

class InverseHyperbolicSecant : public Function {

};

class InverseHyperbolicCotangent : public Function {

};

#endif