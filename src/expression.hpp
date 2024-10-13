#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <variant>
#include <unordered_map>
#include <string>
#include <optional>

struct Expression {
    // Enum to represent the operation types
    enum class OperationType {
        ADD = 4,
        SUBTRACT = 5,
        MULTIPLY = 6,
        NEGATIVE = 7,
        AND = 0,
        OR = 1,
        XOR = 2,
        NOT = 3,
    };

    // Forward declaration(s)
    struct Node;
    struct Operation;

    /* BEGIN EXPRESSION CLASS */

    std::shared_ptr<Node> root;  // The root node of the expression
    std::vector<std::string> variables;  // Variable IDs used in the expression
    std::unordered_map<std::string, int> variableMap; // Maps variable names to indices
    size_t n_variables = 0;

    // Constructor to initialize the root node
    Expression(std::shared_ptr<Node> root) : root(root) {}
    Expression() {}

    // Add a variable to the expression
    void addVariable(const std::string& name) {
        variableMap[name] = n_variables;
        variables.push_back(name);
        n_variables++;
    }

    // Print the expression starting from the root node
    void print(bool isPolishNotation = false) {
        printNode(root, this, isPolishNotation);
    }

    /* END EXPRESSION CLASS */

    struct Operation {
        OperationType type; // Type of the operation
        std::shared_ptr<Node> left;  // Left operand 
        std::shared_ptr<Node> right; // Right operand (null for unary operations)
        bool isBooleanArithmetic; // Determines if the operation is boolean arithmetic
        bool isLinearMBA; // Determines if the operation is a linear MBA
        bool isConstant;

        Operation(OperationType type, std::shared_ptr<Node> left, std::shared_ptr<Node> right = nullptr)
            : type(type), left(left), right(right), 
            isBooleanArithmetic(calculateBooleanArithmetic()),
            isLinearMBA(calculateLinearMBA()),
            isConstant(calculateConstant()) {}

    private:
        bool calculateBooleanArithmetic() {
            // Check the type of the operation
            switch (type) {
                case OperationType::AND:
                case OperationType::OR:
                case OperationType::XOR:
                case OperationType::NOT:
                    if (left && !left->isBooleanArithmetic) return false;
                    if (right && !right->isBooleanArithmetic) return false;
                    return true;
                default:
                    return false; // Other operations are not boolean arithmetic
            }
        }

        bool calculateLinearMBA() {
            // Check if the operation is ADD or SUBTRACT
            if (type == OperationType::ADD || type == OperationType::SUBTRACT) {
                if (!left || !left->isLinearMBA) return false;
                if (!right || !right->isLinearMBA) return false;
                return true; // Both operands are valid linear MBAs
            }

            // If the operation is MULTIPLY, one or both sides must be a scalar
            if (type == OperationType::MULTIPLY) {
                // Check if either the left or right operand is a scalar
                if (left && std::holds_alternative<uint64_t>(left->value)) return right->isLinearMBA;
                if (right && std::holds_alternative<uint64_t>(right->value)) return left->isLinearMBA;
                return false; // Neither side is a scalar
            }

            return isBooleanArithmetic; // Boolean expressions are also linear MBAs
        }

        bool calculateConstant() {
            if ((left && !left->isConstant) || (right && !right->isConstant)) {
                return false;
            }
            return true;
        }
    };

    struct Node {
        std::variant<uint64_t, std::shared_ptr<Operation>, int> value; // Can hold a scalar, an operation, or a variable ID
        std::vector<std::shared_ptr<Node>> alternateRepresentations;   // Alternate representations
        bool isBooleanArithmetic; // Determines if the node is part of boolean arithmetic
        bool isLinearMBA; // Determines if the node is a Linear 
        bool isConstant; // Determines if the node is a constant
        std::optional<std::vector<uint64_t>> bitVector; // Optional bit vector

        // Constructor for scalar value
        Node(uint64_t scalarValue, Expression* expr = nullptr) : value(scalarValue), isConstant(true) {
            isBooleanArithmetic = calculateBooleanArithmetic();
            isLinearMBA = calculateLinearMBA();
            if (expr) {
                bitVector = std::vector<uint64_t>(1 << expr->n_variables, -scalarValue);
            }
        }

        // Constructor for variable
        Node(int variableID, Expression* expr = nullptr) : value(variableID) {
            isBooleanArithmetic = true; // Variables are considered boolean arithmetic
            isLinearMBA = true; // Variables are also Linear MBAs
            isConstant = false; // Variables are not constants
            if (expr) {
                auto& n = expr->n_variables;
                bitVector = std::vector<uint64_t>(1 << n, 0);
                for(int i = 0; i < (1 << n); i++) {
                    if (i >> variableID & 1) (*bitVector)[i] = 1;
                }
            }
        }

        // Constructor for Operation
        Node(std::shared_ptr<Operation> op, Expression* expr = nullptr)
            : value(op), 
            isBooleanArithmetic(op->isBooleanArithmetic), 
            isLinearMBA(op->isLinearMBA),
            isConstant(op->isConstant) {
            
            if (expr) {
                bitVector = calculateBitVector(op->type, op->left, op->right, expr);
            }
        }

        // Constructor for operation (directly rather than from an Operation)
        Node(OperationType type, std::shared_ptr<Node> left, std::shared_ptr<Node> right = nullptr, Expression* expr = nullptr) {
            auto op = std::make_shared<Operation>(type, left, right);
            value = op;
            isBooleanArithmetic = op->isBooleanArithmetic;
            isLinearMBA = op->isLinearMBA;
            isConstant = op->isConstant;

            if (expr) bitVector = calculateBitVector(type, left, right, expr);
        }

        // Add an alternate representation
        void addAlternateRepresentation(std::shared_ptr<Node> alternate) {
            alternateRepresentations.push_back(alternate);
        }

    private:
        // Method to determine if the node is part of boolean arithmetic
        bool calculateBooleanArithmetic() const {
            if (std::holds_alternative<uint64_t>(value)) {
                // Check if it's a boolean constant (0 or -1)
                uint64_t constant = std::get<uint64_t>(value);
                return (constant == 0 || constant == static_cast<uint64_t>(-1));
            } 
            else if (std::holds_alternative<int>(value)) {
                return true; // Variables are boolean arithmetic
            } 
            else {
                // It's an operation
                auto operation = std::get<std::shared_ptr<Operation>>(value);
                return operation->isBooleanArithmetic; // Use the property from the Operation
            }
            return false; // Fallback
        }

        // Method to determine if the node is a Linear MBA
        bool calculateLinearMBA() const {
            if (std::holds_alternative<uint64_t>(value)) {
                return true; // Constants are Linear MBAs
            } else if (std::holds_alternative<int>(value)) {
                return true; // Variables are Linear MBAs
            } else {
                // It's an operation
                auto operation = std::get<std::shared_ptr<Operation>>(value);
                return operation->isLinearMBA; // Use the property from the Operation
            }
            return false; // Fallback
        }

        // Known issue: Combinations of variables don't have bit vectors or are considered linear, when they should be
        std::vector<uint64_t> calculateBitVector(OperationType type, std::shared_ptr<Node> left, std::shared_ptr<Node> right, Expression* expr) {
            auto bitVector = std::vector<uint64_t>(1 << expr->n_variables, 0);
            
            // Handle NOT operators separately
            if (type == OperationType::NOT) {
                if (isBooleanArithmetic && left && left->bitVector) {
                    for (size_t i = 0; i < (1 << expr->n_variables); i++) {
                        bitVector[i] = -(~(-(*left->bitVector)[i])); // Invert the bits
                    }
                }
            }
            // Handle NEGATIVE operators separately
            else if (type == OperationType::NEGATIVE) {
                if (isLinearMBA && left && left->bitVector) {
                    for (size_t i = 0; i < (1 << expr->n_variables); i++) {
                        bitVector[i] = -(*left->bitVector)[i]; // Negate the bits
                    }
                }
            }
            // Handle boolean arithmetic operations
            else if (isBooleanArithmetic) {
                if (left && left->bitVector && right && right->bitVector) {
                    for (size_t i = 0; i < (1 << expr->n_variables); i++) {
                        uint64_t leftVal = (*left->bitVector)[i];
                        uint64_t rightVal = (*right->bitVector)[i];
                        if (isConstant) {
                            leftVal = -leftVal;
                            rightVal = -rightVal;
                        }
                        switch (type) {
                            case OperationType::AND: bitVector[i] = leftVal & rightVal; break;
                            case OperationType::OR: bitVector[i] = leftVal | rightVal; break;
                            case OperationType::XOR: bitVector[i] = leftVal ^ rightVal; break;
                            default: break; // Ignore other types for now
                        }
                        if (isConstant) bitVector[i] = -bitVector[i];
                    }
                }
            }
            // Handle linear MBA operations
            else if (isLinearMBA) {
                if (left && left->bitVector && right && right->bitVector) {
                    for (size_t i = 0; i < (1 << expr->n_variables); i++) {
                        uint64_t leftVal = (*left->bitVector)[i];
                        uint64_t rightVal = (*right->bitVector)[i];
                        if (left->isConstant) leftVal = -leftVal;
                        if (right->isConstant) rightVal = -rightVal;
                        switch (type) {
                            case OperationType::ADD: bitVector[i] = leftVal + rightVal; break;
                            case OperationType::SUBTRACT: bitVector[i] = leftVal - rightVal; break;
                            case OperationType::MULTIPLY: bitVector[i] = leftVal * rightVal; break;
                            default: break; // Ignore other types for now
                        }
                        if (isConstant) bitVector[i] = -bitVector[i];
                    }
                }
            }
            return bitVector;
        }
    };

    // Helper function to print the operation type as a string
    static std::string operationTypeToString(OperationType type) {
        switch (type) {
            case OperationType::ADD: return "+";
            case OperationType::SUBTRACT: return "-";
            case OperationType::MULTIPLY: return "*";
            case OperationType::NEGATIVE: return "-";
            case OperationType::AND: return "&";
            case OperationType::OR: return "|";
            case OperationType::XOR: return "^";
            case OperationType::NOT: return "~";
            default: return "?";
        }
    }

    // Prints a node in Polish or normal notation.
    static void printNode(const std::shared_ptr<Node>& node, Expression* expr = nullptr, bool isPolishNotation = false) {
        if (std::holds_alternative<uint64_t>(node->value)) {
            std::cout << std::get<uint64_t>(node->value);
        } else if (std::holds_alternative<int>(node->value)) {
            if (!expr) {
                std::cout << "x" << std::get<int>(node->value); // Print variables as x0, x1, x2, etc.
            } else {
                std::cout << expr->variables[std::get<int>(node->value)]; // Print variables with their names
            }
        } else {
            auto operation = std::get<std::shared_ptr<Operation>>(node->value);
            bool isPrefixOperator = isPolishNotation || operation->type == OperationType::NOT || operation->type == OperationType::NEGATIVE;
            if (isPrefixOperator) {
                std::cout << operationTypeToString(operation->type) << " ";
            }
            std::cout << "(";
            if (operation->left) {
                printNode(operation->left, expr);
            }
            if (!isPrefixOperator) {
                std::cout << " " << operationTypeToString(operation->type);
            }
            if (operation->right) {
                std::cout << " ";
                printNode(operation->right, expr);
            }
            std::cout << ")";
        }
    }
};

// For convenience when we work later
typedef Expression::Node Node;
typedef Expression::Operation Operation;
typedef Expression::OperationType OperationType;