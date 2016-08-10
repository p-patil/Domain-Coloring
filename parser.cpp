#include <algorithm>
#include <sstream>
#include "parser.h"
#include "complex_numbers.h"
#include "functions.h"

using namespace std;

// Begin parent expression tree node class.

ExpressionTreeNode::ExpressionTreeNode() {}

// End parent expression tree node class.

// Begin parent expression tree node class.

ExpressionTreeBinaryOp::ExpressionTreeBinaryOp(void) {}

ExpressionTreeBinaryOp::ExpressionTreeBinaryOp(string op) {    this->op = op; }

ExpressionTreeBinaryOp::ExpressionTreeBinaryOp(string op, ExpressionTreeNode *l, ExpressionTreeNode *r) {
    this->op = op;
    this->left = l;
    this->right = r;
}

ExpressionTreeBinaryOp::~ExpressionTreeBinaryOp() {
    delete left;
    delete right;
}

string ExpressionTreeBinaryOp::get_operator(void) const { return this->op; }

void ExpressionTreeBinaryOp::set_operator(string op) { this->op = op; }

ExpressionTreeNode * ExpressionTreeBinaryOp::get_left(void) const { return this->left; }

ExpressionTreeNode * ExpressionTreeBinaryOp::get_right(void) const { return this->right; }

void ExpressionTreeBinaryOp::set_left(ExpressionTreeNode *l) { this->left = l; }

void ExpressionTreeBinaryOp::set_right(ExpressionTreeNode *r) { this->right = r; }

bool ExpressionTreeBinaryOp::is_leaf_node(void) const { return false; }

bool ExpressionTreeBinaryOp::is_op_node(void) const { return true; }

bool ExpressionTreeBinaryOp::is_func_node(void) const { return false; }

bool ExpressionTreeBinaryOp::is_var_node(void) const { return false; }

bool ExpressionTreeBinaryOp::is_const_node(void) const { return false; }

// End parent expression tree node class.

// Begin parent expression tree node class.

ExpressionTreeFunction::ExpressionTreeFunction(void) {}

ExpressionTreeFunction::ExpressionTreeFunction(string f) { this->function = f; }

ExpressionTreeFunction::ExpressionTreeFunction(ExpressionTreeNode *a) { this->argument = a; }

ExpressionTreeFunction::ExpressionTreeFunction(ExpressionTreeNode *a, string f) {
    this->argument = a;
    this->function = f;
}

ExpressionTreeFunction::~ExpressionTreeFunction() { delete argument; }

ExpressionTreeNode * ExpressionTreeFunction::get_argument(void) const { return this->argument; }

string ExpressionTreeFunction::get_function(void) const { return this->function; }

void ExpressionTreeFunction::set_argument(ExpressionTreeNode *a) { this->argument = a; }

void ExpressionTreeFunction::set_function(const string f) { this->function = f; }

bool ExpressionTreeFunction::is_leaf_node(void) const { return false; }

bool ExpressionTreeFunction::is_op_node(void) const { return false; }

bool ExpressionTreeFunction::is_func_node(void) const { return true; }

bool ExpressionTreeFunction::is_var_node(void) const { return false; }

bool ExpressionTreeFunction::is_const_node(void) const { return false; }

// End parent expression tree node class.

// Begin parent expression tree node class.

ExpressionTreeLeaf::ExpressionTreeLeaf(void) {}

ExpressionTreeLeaf::ExpressionTreeLeaf(ComplexNumber v) { this->val = v; }

ExpressionTreeLeaf::~ExpressionTreeLeaf() { this->val.~ComplexNumber(); }

ComplexNumber ExpressionTreeLeaf::get_val(void) { return this->val; }

void ExpressionTreeLeaf::set_val(ComplexNumber v) { this->val = v; }

bool ExpressionTreeLeaf::is_leaf_node(void) const { return true; }

bool ExpressionTreeLeaf::is_op_node(void) const { return false; }

bool ExpressionTreeLeaf::is_func_node(void) const { return false; }

bool ExpressionTreeLeaf::is_var_node(void) const { return false; }

bool ExpressionTreeLeaf::is_const_node(void) const { return false; }

// End parent expression tree node class.

// Begin constant expression tree class.

ExpressionTreeConstant::ExpressionTreeConstant(string c) { 
    if (!CONSTANTS.count(c)) {
        throw invalid_argument("Must be known mathematical constant");
    }

    this->constant = c;
}

string ExpressionTreeConstant::get_constant(void) const { return this->constant; }

void ExpressionTreeConstant::set_constant(const string c) { this->constant = c; }

bool ExpressionTreeConstant::is_leaf_node(void) const { return false; }

bool ExpressionTreeConstant::is_op_node(void) const { return false; }

bool ExpressionTreeConstant::is_func_node(void) const { return false; }

bool ExpressionTreeConstant::is_var_node(void) const { return false; }

bool ExpressionTreeConstant::is_const_node(void) const { return true; }

// End constant expression tree class.

// Begin var expression tree node class.

ExpressionTreeVariable::ExpressionTreeVariable(void) {}

bool ExpressionTreeVariable::is_leaf_node(void) const { return false; }

bool ExpressionTreeVariable::is_op_node(void) const { return false; }

bool ExpressionTreeVariable::is_func_node(void) const { return false; }

bool ExpressionTreeVariable::is_var_node(void) const { return true; }

bool ExpressionTreeVariable::is_const_node(void) const { return false; }

// End var expression tree node class.

// TODO: Expand to allow for more functions.
// TODO: When parsing e^(f(z)), use Exponential class instead of Power class.
// Evaluates the given expression tree on the given variable value.
ComplexNumber evaluate_tree(const ExpressionTreeNode *root, const ComplexNumber z) {
    if (root->is_op_node()) {
        ExpressionTreeBinaryOp *temp = (ExpressionTreeBinaryOp *) root;
        ComplexNumber left = evaluate_tree(temp->get_left(), z);
        ComplexNumber right = evaluate_tree(temp->get_right(), z);

        // Find the appropriate operator, apply, then return.
        if (temp->get_operator() == PLUS) {
            return left + right;
        } else if (temp->get_operator() == MINUS) {
            return left - right;
        } else if (temp->get_operator() == TIMES) {
            return left * right;
        } else if (temp->get_operator() == DIVIDE) {
            return left / right;
        } else { // temp->get_operator() == EXP
            Power p (right);
            return p.eval(left);
        }
    } else if (root->is_func_node()) {
        ExpressionTreeFunction *temp = (ExpressionTreeFunction *) root;
        ComplexNumber arg = evaluate_tree(temp->get_argument(), z);

        // Find the appropriate function, evaluate, then return.
        if (temp->get_function() == LOG) {
            Logarithm l (1, 0);
            return l.eval(arg);
        } else if (temp->get_function() == SIN) {
            Sine s (1, 1, 0);
            return s.eval(arg);
        } else if (temp->get_function() == COS) {
            Cosine c (1, 1, 0);
            return c.eval(arg);
        } else { // temp->get_function() == TAN
            Tangent t (1, 1, 0);
            return t.eval(arg);
        }
    } else if (root->is_leaf_node()) {
        ExpressionTreeLeaf *temp = (ExpressionTreeLeaf *) root;
        return temp->get_val(); // Return the value itself
    } else if (root->is_const_node()) {
        ExpressionTreeConstant *temp = (ExpressionTreeConstant *) root;

        // Return the corresponding mathematical constant
        if (temp->get_constant() == e) {
            return E;
        } else if (temp->get_constant() == pi) {
            return PI;
        } else { // temp->get_constant() == phi
            return PHI;
        }
    } else { // root->is_var_node()
        return z;
    }
}

// TODO: Allow for logarithms with arbitrary bases, by parsing for '_' character after 'log' token.
// Primary function that parses the given mathematical expression. Assumes the expression is a single-variable with symbol 'z'.
ExpressionTreeNode * parse(string expr) { return parse_helper(tokenize_shallow(remove_whitespace(expr))); }

// TODO: Write this function
// Simplifies the given expression tree as much as possible using standard rules of algebra. "Simplifying" in this context means
// returning an equivalent expression tree that is as small (in the number of nodes) as possible.
ExpressionTreeNode * simplify(ExpressionTreeNode *root) {

}

// HELPER FUNCTIONS

ExpressionTreeNode * parse_helper(vector<string> shallow_tokens) {
    if (shallow_tokens.size() == 1 && is_base_expression(shallow_tokens[0])) {
        if (is_complex_number(shallow_tokens[0])) {
            ExpressionTreeLeaf *root = new ExpressionTreeLeaf(to_complex_number(shallow_tokens[0]));
            return root;
        } else if (CONSTANTS.count(shallow_tokens[0])) {
            ExpressionTreeConstant *root = new ExpressionTreeConstant(shallow_tokens[0]);
            return root;
        } else { // shallow_tokens[0] == VAR
            ExpressionTreeVariable *root = new ExpressionTreeVariable();
            return root;
        }
    }

    int index = -1; // Index of lowest-precedence operator

    // Iterate through the operators in ascending order by precedence
    for (int i = 0; i < OPERATORS.size(); i++) {
        for (int j = 0; j < shallow_tokens.size(); j++) {
            if (OPERATORS[i] == shallow_tokens[j]) {
                index = j;
                break;
            }
        }

        if (index != -1) {
            break;
        }
    }

    if (index == -1) { // All tokens are either non-shallow or functions - first token must be function
        // Recursively tokenize the non-shallow function argument, set as child, and return.
        ExpressionTreeFunction *root = new ExpressionTreeFunction(shallow_tokens[0]);
        ExpressionTreeNode *arg = parse_helper(tokenize_shallow(shallow_tokens[1]));
        root->set_argument(arg);
        return root;
    } else {
        ExpressionTreeBinaryOp *root = new ExpressionTreeBinaryOp(shallow_tokens[index]);

        // Build left and right token lists
        vector<string> left, right;
        
        for (int i = 0; i < index; i++) {
            left.push_back(shallow_tokens[i]);
        }

        for (int i = index + 1; i < shallow_tokens.size(); i++) { // Skip the operator
            right.push_back(shallow_tokens[i]);
        }

        // Recursively resolve binary parameters, then return
        ExpressionTreeNode *left_child = parse_helper(left), *right_child = parse_helper(right);
        root->set_left(left_child);
        root->set_right(right_child);
        return root;
    }
}

// Returns a vector of tokens from the given expression. The returned vector is shallow in that expressions grouped by 
// parentheses (including function arguments) are considered entire tokens. Assumes expr has no whitespace.
// Since the returned tokens are shallow, pairs of tokens are returned, where each token is bundled with a boolean 
// specifying if the corresponding token is shallow or not.
vector<string> tokenize_shallow(string expr) {
    vector<string> tokens; // Stores tokens
    string s; // Temporary buffer string
    int j = 0, parentheses_count;
    for (int i = 0; i < expr.length(); j++) {
        s = expr.substr(i, 1);
        if (find(OPERATORS.begin(), OPERATORS.end(), s) != OPERATORS.end()) { // Check if operator (use substr since operator is unordered_set<string>)
            tokens.push_back(s);
            i++;
        } else if (s == VAR) { // Check if variable
            tokens.push_back(s);
            i++;
        } else if (isdigit(expr[i])) { // Check if number
            bool decimal_flag = false;
            s = "";
            while (isdigit(expr[i])) {
                s += expr[i];
                i++;
            }

            if (expr[i] == '.') {
                decimal_flag = true;
                i++;
            }

            while (isdigit(expr[i])) {
                s += expr[i];
                i++;
            }

            tokens.push_back(s);
        } else if (expr[i] == '(') { // If parenthetical expression, group as token
            i++; // Skip first '('
            s = "";
            parentheses_count = 1;
            while (parentheses_count > 0) {
                s += expr[i];

                if (expr[i] == '(') {
                    parentheses_count++;
                } else if (expr[i] == ')') {
                    parentheses_count--;
                }

                i++;
            }

            s = s.substr(0, s.length() - 1); // Discard the final ')' at the end
            tokens.push_back(s);
        } else {
            // Check if function
            bool flag = false;
            for (auto itr = FUNCTIONS.begin(); itr != FUNCTIONS.end(); itr++) {
                s = expr.substr(i, (*itr).length());
                if (s == *itr) {
                    // Found a match
                    flag = true;
                    tokens.push_back(s);
                    i += s.length();
                    break;
                }
            }

            // If the token is a function, skip the entire argument as only shallow tokens are considered.
            if (flag) {
                if (expr[i] != '(') {
                    ostringstream stream;
                    stream << "Invalid expression - expected '(' at index " << i << "in \"" << expr << "\"";
                    throw invalid_argument(stream.str());
                }

                i++; // Skip first '('
                j++;
                s = "";
                parentheses_count = 1;
                while (parentheses_count > 0) {
                    s += expr[i];

                    if (expr[i] == '(') {
                        parentheses_count++;
                    } else if (expr[i] == ')') {
                        parentheses_count--;
                    }

                    i++;
                }

                s = s.substr(0, s.length() - 1); // Discard the final ')' at the end
                tokens.push_back(s);
            } else { // Check if mathematical constant
                for (auto itr = CONSTANTS.begin(); itr != CONSTANTS.end(); itr++) {
                    s = expr.substr(i, (*itr).length());
                    if (s == *itr) {
                        // Found a match
                        flag = true;
                        tokens.push_back(s);
                        i += s.length();
                        break;
                    }
                }
            }

            if (!flag) { // Invalid token
                ostringstream stream;
                stream << "Invalid expression - failed at index " << i << " ('" << expr[i] << "') in \"" << expr << "\"";
                throw invalid_argument(stream.str());
            }
        }
    }

    return tokens;
}

string remove_whitespace(const string s) {
    string t (s);
    t.erase(remove_if(t.begin(), t.end(), ::isspace), t.end());
    return t;
}

// Returns if the given expression contains any sub-expressions or not.
bool is_base_expression(const string expr) { return is_number(expr) || is_complex_number(expr) || expr == VAR || CONSTANTS.count(expr); }

// Returns if the given string is a (decimal) number.
bool is_number(const string s) {
    bool decimal_point = false;
    for (int i = (s[0] == '-') ? 1 : 0; i < s.length(); i++) {
        if (s[i] == '.') {
            if (decimal_point) {
                return false;
            } else {
                decimal_point = true;
            }
        } else if (!isdigit(s[i])) {
            return false;
        }
    }
    
    return true;
}

// Converts the given string to a double.
double to_number(const string s) {
    int sign = (s[0] == '-') ? -1 : 1;
    int i = (s[0] == '-') ? 1 : 0;

    double pow_ten = 1;
    for (int j = i; j < s.length() && s[j] != '.'; j++) {
        pow_ten *= 10;
    }

    double n = 0;
    for (pow_ten = pow_ten / 10; i < s.length(); i++, pow_ten /= 10) {
        if (s[i] == '.') {
            pow_ten *= 10;
            continue;
        }

        n += (s[i] - '0') * pow_ten;
    }

    return n;
}

// Returns if the given string is a complex number.
bool is_complex_number(const string s) {
    if (is_number(s)) {
        return true;
    }

    string t = remove_whitespace(s);

    int index = -1;
    for (int i = 0; i < t.length(); i++) {
        if (t[i] == '+' || t[i] == '-') {
            index = i;
            break;
        }
    }

    if (index == -1) {
        return false;
    } else {
        string real = "", imaginary = "";

        for (int i = 0; i < index; i++) {
            real += t[i];
        }

        for (int i = index; i < t.length(); i++) {
            imaginary += t[i];
        }

        if (real.length() == 0 && imaginary.length() == 0) {
            return false;
        } else if (real.length() != 0) {
            return is_number(real);
        } else { // imaginary.length() != 0
            return is_number(t.substr(0, t.length() - 1)) && (t[t.length() - 1] == 'i');
        }
    }
}

// Converts the given string representation of a complex number to a complex number.
ComplexNumber to_complex_number(const string s) {
    if (is_number(s)) {
        ComplexNumber z(to_number(s), 0.0);
        return z;
    }

    string t = remove_whitespace(s);

    int index = -1;
    for (int i = 0; i < t.length(); i++) {
        if (t[i] == '+' || t[i] == '-') {
            index = i;
            break;
        }
    }

    string real = "", imaginary = "";

    for (int i = 0; i < index; i++) {
        real += t[i];
    }

    for (int i = index + 1; i < t.length() - 1; i++) { // Skip the + / - and the i at the end
        imaginary += t[i];
    }

    if (real.length() == 0) {
        ComplexNumber z (0.0, to_number(imaginary));
        return z;
    } else if (imaginary.length() == 0) {
        ComplexNumber z (to_number(real), 0.0);
        return z;
    } else {
        ComplexNumber z (to_number(real), to_number(imaginary));
        return z;
    }
}

// Given two operators, returns the operator precedence: 1 if op1 has higher precedence, 0 if both are same, and -1 if op2
// has higher precedence.
int operator_precedence(string op1, string op2) {
    int op1_hash, op2_hash;

    // Assign op1 a value equivalent to its position in the precedence hierarchy
    // if (op1 == LOG || op1 == SIN || op1 == COS || op1 == TAN) { // Functions have lowest precedence
    if (FUNCTIONS.count(op1)) { // Functions have lowest precedence
        op1_hash = 0;
    } else if (op1 == PLUS || op1 == MINUS) { // Then addition and subtraction
        op1_hash = 1;
    } else if (op1 == TIMES || op1 == DIVIDE) { // Then multiplication and division
        op1_hash = 2;
    } else if (op1 == EXP) { // Finally, exponentiation.
        op1_hash = 3;
    } else {
        throw new invalid_argument("Argument must be an operator or function");
    }

    // Repeat for op2
    // if (op2 == LOG || op2 == SIN || op2 == COS || op2 == TAN) { // Functions have lowest precedence
    if (FUNCTIONS.count(op2)) { // Functions have lowest precedence
        op2_hash = 0;
    } else if (op2 == PLUS || op2 == MINUS) { // Then addition and subtraction
        op2_hash = 1;
    } else if (op2 == TIMES || op2 == DIVIDE) { // Then multiplication and division
        op2_hash = 2;
    } else if (op2 == EXP) { // Finally, exponentiation.
        op1_hash = 3;
    } else {
        throw new invalid_argument("Argument must be an operator or function");
    }

    // Return
    if (op1_hash > op2_hash) {
        return 1;
    } else if (op1_hash == op2_hash) {
        return 0;
    } else {
        return -1;
    }
}