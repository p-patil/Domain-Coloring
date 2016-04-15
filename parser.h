#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include "complex_numbers.h"
#include "functions.h"

using namespace std;

// Identifiers for operators, functions, mathematical constants
const string VAR = "z";
const string PLUS = "+";
const string MINUS = "-";
const string TIMES = "*";
const string DIVIDE = "/";
const string EXP = "^";
const string LOG = "log";
const string SIN = "sin";
const string COS = "cos";
const string TAN = "tan";
const string e = "e";
const string pi = "pi";
const string phi = "phi";

const vector<string> OPERATORS = {PLUS, MINUS, TIMES, DIVIDE, EXP}; // Set of operators, ordered (ascending) by precedence
const unordered_set<string> FUNCTIONS = {LOG, SIN, COS, TAN};
const unordered_set<string> CONSTANTS = {e, pi, phi};

class ExpressionTreeNode {
	public:
		ExpressionTreeNode(void);

		bool is_op_node(void) const;

		bool is_func_node(void) const;

		bool is_leaf_node(void) const;

		bool is_var_node(void) const;

		bool is_const_node(void) const;
};

class ExpressionTreeBinaryOp : public ExpressionTreeNode {
	ExpressionTreeNode *left, *right;
	string op;

	public:
		ExpressionTreeBinaryOp(void);

		ExpressionTreeBinaryOp(string);

		ExpressionTreeBinaryOp(string, ExpressionTreeNode *, ExpressionTreeNode *);

		~ExpressionTreeBinaryOp();

		string get_operator(void) const;

		void set_operator(string);

		ExpressionTreeNode * get_left(void) const;

		ExpressionTreeNode * get_right(void) const;

		void set_left(ExpressionTreeNode *);

		void set_right(ExpressionTreeNode *);

		bool is_op_node(void) const;
};

class ExpressionTreeFunction : public ExpressionTreeNode {
	ExpressionTreeNode *argument;
	string function;

	public:
		ExpressionTreeFunction(void);
		
		ExpressionTreeFunction(string);

		ExpressionTreeFunction(ExpressionTreeNode *);
		
		ExpressionTreeFunction(ExpressionTreeNode *, string);

		~ExpressionTreeFunction();

		ExpressionTreeNode * get_argument(void) const;

		string get_function(void) const;

		void set_argument(ExpressionTreeNode *);

		void set_function(const string);

		bool is_func_node(void) const;
};

class ExpressionTreeLeaf : public ExpressionTreeNode {
	ComplexNumber val;

	public:
		ExpressionTreeLeaf(void);

		ExpressionTreeLeaf(ComplexNumber);

		~ExpressionTreeLeaf();

		ComplexNumber get_val(void);

		void set_val(ComplexNumber);

		bool is_leaf_node(void) const;
};

class ExpressionTreeConstant : public ExpressionTreeNode {
	string constant;

	public:
		ExpressionTreeConstant(string);

		string get_constant(void) const;

		void set_constant(const string);

		bool is_const_node(void) const;
};

class ExpressionTreeVariable : public ExpressionTreeNode {
	public:
		ExpressionTreeVariable(void);

		bool is_var_node(void) const;
};

vector<Function *> to_elementary_composition(ExpressionTreeNode *);

ExpressionTreeNode * parse(string);

ExpressionTreeNode * simplify(ExpressionTreeNode *);

// HELPER FUNCTIONS

void to_elementary_composition_helper(ExpressionTreeNode *, vector<Function *>);

ExpressionTreeNode * parse_helper(vector<string>);

vector<string> tokenize_shallow(string);

string remove_whitespace(string);

bool is_base_expression(const string);

bool is_number(const string);

double to_number(const string);

bool is_complex_number(const string);

ComplexNumber to_complex_number(const string);

int operator_precedence(string, string);