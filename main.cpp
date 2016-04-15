#include <iostream>
#include <vector>
#include "complex_numbers.h"
#include "functions.h"
#include "grapher.h"
#include "parser.h"
#include "SDL/SDL.h"

using namespace std;

ComplexNumber identity(ComplexNumber z) { return z; }

ComplexNumber square(ComplexNumber z) { return z * z; }

ComplexNumber exponential(ComplexNumber z) {
	Exponential exp (1, 1);
	return exp.eval(z);
}

ComplexNumber inverse(ComplexNumber z) { return 1 / z; }

ComplexNumber sine(ComplexNumber z) {
	Sine s (1, 1, 0);
	return s.eval(z);
}

ComplexNumber polynomial(ComplexNumber z) { return (z * z * z - 1) * (z - 1); }

int main(int argc, char **argv) {
	freopen("CON", "w", stdout); // Prevents redirect of output
	freopen("CON", "w", stderr); // and errors...
	if (SDL_Init(SDL_INIT_EVERYTHING) == -1) {
		cerr << "Failed to initialize environment\n";
		return -1;
	}

	// ExpressionTreeBinaryOp *expr = (ExpressionTreeBinaryOp *) parse("log(5 * z) + e^10");
	// cout << expr->get_operator() << endl; // expect +
	// cout << ((ExpressionTreeFunction *) expr->get_left())->get_function() << endl; // expect sin
	// cout << ((ExpressionTreeBinaryOp *) ((ExpressionTreeFunction *) expr->get_left())->get_argument())->get_operator() << endl; // expect *
	// cout << ((ExpressionTreeLeaf *) ((ExpressionTreeBinaryOp *) ((ExpressionTreeFunction *) expr->get_left())->get_argument())->get_left())->get_val() << endl; // expect 5
	// cout << ((ExpressionTreeVariable *) ((ExpressionTreeBinaryOp *) ((ExpressionTreeFunction *) expr->get_left())->get_argument())->get_right())->is_var_node() << endl; // expect 1 (true)
	// cout << ((ExpressionTreeBinaryOp *) expr->get_right())->get_operator() << endl; // expect ^
	// cout << ((ExpressionTreeConstant *) ((ExpressionTreeBinaryOp *) expr->get_right())->get_left())->get_constant() << endl; // expect e
	// cout << ((ExpressionTreeLeaf *) ((ExpressionTreeBinaryOp *) expr->get_right())->get_right())->get_val() << endl; // expect 10
	// delete expr;

	vector<Function *> arr;
	
	ComplexNumber z (1, 1);
	Constant c (z);
	arr.push_back(&c);

	Exponential expo (1, 1);
	arr.push_back(&expo);

	Logarithm loga (1, 0);
	arr.push_back(&loga);


	for (int i = 0; i < 3; i++) {
		cout << *(arr[i]) << endl;
	}

	SDL_Quit();
	return 0;
}