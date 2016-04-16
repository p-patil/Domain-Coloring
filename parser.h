#include <unordered_set>
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
const string CSC = "csc";
const string SEC = "sec";
const string COT = "cot";
const string ARCSIN = "arcsin";
const string ARCCOS = "arccos";
const string ARCTAN = "arctan";
const string ARCCSC = "arccsc";
const string ARCSEC = "arcsec";
const string ARCCOT = "arccot";
const string SINH = "sinh";
const string COSH = "cosh";
const string TANH = "tanh";
const string CSCH = "csch";
const string SECH = "sech";
const string COTH = "coth";
const string ARCSINH = "arcsinh";
const string ARCCOSH = "arccosh";
const string ARCTANH = "arctanh";
const string ARCCSCH = "arccsch";
const string ARCSECH = "arcsech";
const string ARCCOTH = "arccoth";

const string e = "e";
const string pi = "pi";
const string phi = "phi";

const vector<string> OPERATORS = {PLUS, MINUS, TIMES, DIVIDE, EXP}; // Set of operators, ordered (ascending) by precedence
const unordered_set<string> FUNCTIONS = {LOG, SIN, COS, TAN, CSC, SEC, COT, ARCSIN, ARCCOS, ARCTAN, ARCCSC, ARCSEC, 
									 	 ARCCOT, SINH, COSH, TANH, CSCH, TANH, CSCH, SECH, COTH, ARCSINH, ARCCOSH, 
									 	 ARCTANH, ARCCSCH, ARCSECH, ARCCOTH};
const unordered_set<string> CONSTANTS = {e, pi, phi};

class ExpressionTreeNode {
	public:
		ExpressionTreeNode(void);

		virtual bool is_op_node(void) const = 0;

		virtual bool is_func_node(void) const = 0;

		virtual bool is_leaf_node(void) const = 0;

		virtual bool is_var_node(void) const = 0;

		virtual bool is_const_node(void) const = 0;
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

		virtual bool is_op_node(void) const override;

		virtual bool is_func_node(void) const override;

		virtual bool is_leaf_node(void) const override;

		virtual bool is_var_node(void) const override;

		virtual bool is_const_node(void) const override;
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

		virtual bool is_op_node(void) const override;

		virtual bool is_func_node(void) const override;

		virtual bool is_leaf_node(void) const override;

		virtual bool is_var_node(void) const override;

		virtual bool is_const_node(void) const override;
};

class ExpressionTreeLeaf : public ExpressionTreeNode {
	ComplexNumber val;

	public:
		ExpressionTreeLeaf(void);

		ExpressionTreeLeaf(ComplexNumber);

		~ExpressionTreeLeaf();

		ComplexNumber get_val(void);

		void set_val(ComplexNumber);

		virtual bool is_op_node(void) const override;

		virtual bool is_func_node(void) const override;

		virtual bool is_leaf_node(void) const override;

		virtual bool is_var_node(void) const override;

		virtual bool is_const_node(void) const override;
};

class ExpressionTreeConstant : public ExpressionTreeNode {
	string constant;

	public:
		ExpressionTreeConstant(string);

		string get_constant(void) const;

		void set_constant(const string);

		virtual bool is_op_node(void) const override;

		virtual bool is_func_node(void) const override;

		virtual bool is_leaf_node(void) const override;

		virtual bool is_var_node(void) const override;

		virtual bool is_const_node(void) const override;
};

class ExpressionTreeVariable : public ExpressionTreeNode {
	public:
		ExpressionTreeVariable(void);

		virtual bool is_op_node(void) const override;

		virtual bool is_func_node(void) const override;

		virtual bool is_leaf_node(void) const override;

		virtual bool is_var_node(void) const override;

		virtual bool is_const_node(void) const override;
};

ComplexNumber evaluate_tree(const ExpressionTreeNode *, const ComplexNumber);

ExpressionTreeNode * parse(string);

ExpressionTreeNode * simplify(ExpressionTreeNode *);

// HELPER FUNCTIONS

ExpressionTreeNode * parse_helper(vector<string>);

vector<string> tokenize_shallow(string);

string remove_whitespace(string);

bool is_base_expression(const string);

bool is_number(const string);

double to_number(const string);

bool is_complex_number(const string);

ComplexNumber to_complex_number(const string);

int operator_precedence(string, string);