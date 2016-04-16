#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include "complex_numbers.h"

using namespace std;

class Function {
	protected:
		// Virtual function used to allow for polymorphic overload of << operator.
		virtual void print_to_stream(ostream &) const;

	public:
		virtual ComplexNumber eval(const ComplexNumber) const = 0;
	
	friend ostream& operator <<(ostream &, const Function &);

};

ostream& operator <<(ostream &, const Function &);

class Constant : public Function {
	ComplexNumber constant;

	public:
		Constant(const ComplexNumber);

		ComplexNumber get_constant(void) const;

		ComplexNumber eval(const ComplexNumber) const override;

		void print_to_stream(ostream &) const override;
};

// Polynomials.
class Polynomial : public Function {
	int degree;
	ComplexNumber *coefficients; // Index n corresponds to coefficient of x^n

	public:
		Polynomial(ComplexNumber[], int);

		// Destructor.
		~Polynomial(void);
		
		int get_degree(void) const;

		ComplexNumber * get_coefficients(void) const;

		ComplexNumber eval(const ComplexNumber) const override;

		Polynomial operator =(const Polynomial &);

		bool operator ==(const Polynomial &) const;

		bool operator !=(const Polynomial &) const;

		void print_to_stream(ostream &) const override;
};

// Class representing functions of the form f(z) = z^a, for complex number a.
class Power : public Function {
	private:
		ComplexNumber exponent;
		ComplexNumber iPower(double) const;

	public:

		Power(ComplexNumber);

		ComplexNumber get_exponent(void) const;

		ComplexNumber eval(const ComplexNumber) const override;

		void print_to_stream(ostream &) const override;
};

class NthRoot : public Function {
	public:
		int n;

		NthRoot(int);

		int get_n(void) const;

		ComplexNumber eval(const ComplexNumber) const override;

		void print_to_stream(ostream &) const override;
};

// Exponential function and its inverse, the logarithm.
class Exponential : public Function {
	// Models the function f(z) = coefficient * e^(lambda * z)
	ComplexNumber coefficient, lambda;

	public:
		Exponential(ComplexNumber, ComplexNumber);

		ComplexNumber get_coefficient(void) const;

		ComplexNumber get_lambda(void) const ;

		ComplexNumber eval(const ComplexNumber) const override;

		Exponential operator =(const Exponential &);

		bool operator ==(const Exponential &) const;

		bool operator !=(const Exponential &) const;
		
		void print_to_stream(ostream &) const override;
};

// Logarithms. Because the class represents the complex logarithm, which is multivalued, the principal log is taken.
class Logarithm : public Function {
	ComplexNumber base, coefficient, lambda; // Represents coefficient * log_(base) (z + lambda)

	public:

		Logarithm(ComplexNumber, ComplexNumber);
		
		Logarithm(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber get_coefficient(void) const;

		ComplexNumber get_lambda(void) const;

		ComplexNumber get_base(void) const;

		ComplexNumber eval(const ComplexNumber) const override;

		Logarithm operator =(const Logarithm &);

		bool operator ==(const Logarithm &) const;

		bool operator !=(const Logarithm &) const;
		
		void print_to_stream(ostream &) const override;
};

// Trig functions and their inverses.
class Sine : public Function {
	ComplexNumber a, b, c; // Represents a * sin(b * z + c)

	public:
		Sine();

		Sine(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;

		ComplexNumber getB(void) const;

		ComplexNumber getC(void) const;

		ComplexNumber eval(const ComplexNumber) const override;
		
		void print_to_stream(ostream &) const override;
};

class Cosine : public Function {
	ComplexNumber a, b, c; // Represents a * cos(b * z + c)

	public:
		Cosine();

		Cosine(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;

		ComplexNumber getB(void) const;

		ComplexNumber getC(void) const;

		ComplexNumber eval(const ComplexNumber) const;
		
		void print_to_stream(ostream &) const override;
};

class Tangent : public Function {
	ComplexNumber a, b, c; // Represents a * tan(b * z + c)

	public:
		Tangent(ComplexNumber, ComplexNumber, ComplexNumber);

		ComplexNumber getA(void) const;
		
		ComplexNumber getB(void) const;
		
		ComplexNumber getC(void) const;

		ComplexNumber eval(const ComplexNumber) const override;
		
		void print_to_stream(ostream &) const override;
};

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

class Arccotangent : public Function {

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