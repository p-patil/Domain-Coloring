#ifndef COMPLEX_NUMBERS_H
#define COMPLEX_NUMBERS_H

// Mathematical constants
#define PI 3.141592653589793238462643383279502884197169399375
#define E 2.718281828459045235360287471352662497757247093699
#define PHI 1.618033988749894848204586834365638117720309179805

using namespace std;

class ComplexNumber {
	double re, im; // Real and imaginary components.
	double r, theta; // Corresponding polar coordinates.

	private:
		void initialize_polar(double, double);

	public:
		ComplexNumber();

		ComplexNumber(double, double);

		/**
		 * Constructor used to implicitly promote double to ComplexNumber.
		 */
		ComplexNumber(double);

		/**
		 * Returns real part this complex number.
		 */
		double real(void) const;

		/**
		 * Returns the imaginary part of this complex number.
		 */
		double imag(void) const;

		/**
		 * Magnitude of the complex number.
		 */
		double mag(void) const;

		/**
		 * Polar angle of the complex number, in radians.
		 */		 
		double angle(void) const;

		/**
		 * Complex conjugate.
		 */
		ComplexNumber conjugate(void) const;

		bool is_real(void) const;

		bool is_imaginary(void) const;

		void set(double, double);

		/**
		 * Assignment operator.
		 */
		ComplexNumber operator =(const ComplexNumber &);

		bool operator ==(const ComplexNumber &) const;

		bool operator !=(const ComplexNumber &) const;

		/**
		 * Complex addition.
		 */
		ComplexNumber operator +(const ComplexNumber &) const;

		friend ComplexNumber operator +(double, const ComplexNumber &);

		/** 
		 * Comlex subtraction.
		 */
		ComplexNumber operator -(const ComplexNumber &) const;
		
		friend ComplexNumber operator -(double, const ComplexNumber &);

		/**
		 * Complex multiplication.
		 */
		ComplexNumber operator *(const ComplexNumber &) const;

		friend ComplexNumber operator *(double, const ComplexNumber &);

		/**
		 * Complex division.
		 */
		ComplexNumber operator /(const ComplexNumber &) const;

		friend ComplexNumber operator /(double, const ComplexNumber &);

		ComplexNumber operator -() const;

		ComplexNumber operator +=(const ComplexNumber &);

		ComplexNumber operator -=(const ComplexNumber &);

		ComplexNumber operator *=(const ComplexNumber &);

		ComplexNumber operator /=(const ComplexNumber &);
};

/**
 * Control how a complex number is printed.
 */
ostream& operator <<(ostream &, const ComplexNumber &);

// Division by zero exception class, so that the poles of complex functions can be found by catching this exception.
class ComplexDivisionByZeroException : public exception {
	public:
		const char * what() const throw();
};

#endif