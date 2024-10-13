#include <iostream>
#include <vector>
#include <memory>
#include <variant>
#include <unordered_map>
#include <string>
#include <optional>
#include <stdexcept>
#include <random>
#include <cassert>
#include <functional>

#include "matrix.hpp"
#include "linear_congruence.hpp"
#include "expression.hpp"

typedef uint64_t u64;
typedef __int128_t i128;

/*
Basis over the vector space Z_{2^64}^n
Algorithm from https://codeforces.com/blog/entry/98376
*/
struct Basis {
    static constexpr i128 m = (i128)1 << 64;
    size_t n;
    std::vector<Vector> basis;

    Basis(size_t n) : n(n), basis(n, Vector(n)) {}

    // Check if a vector can be represented in the basis
    bool insideBasis(Vector v) {
        if (v.size() != n) {
            throw std::invalid_argument("Vector size does not match basis dimension.");
        }

        for (auto i = 0; i < n; ++i) {
            const auto& b = basis[i];
            if (b[i] && v[i] % b[i] == 0) {
                v -= b * (v[i] / b[i]);
            }
            if (v[i] != 0) return false;
        }

        return true;
    }

    void insert(Vector v) {
        if (v.size() != n) {
            throw std::invalid_argument("Vector size does not match basis dimension.");
        }

        for (auto i = 0; i < n; ++i) {
            if (v[i] != 0) {
                auto& b = basis[i];

                if (b[i] == 0) {
                    auto [_ignore_1, s, _ignore_2] = extended_gcd(v[i], m);
                    b = v * s;
                    v = v * (m / b[i]);
                } else if (v[i] % b[i] == 0) {
                    v -= b * (v[i] / b[i]);
                } else {
                    auto w = b;
                    auto [_, s, t] = extended_gcd(v[i], w[i]);
                    b = v * s + w * t;
                    v -= b * (v[i] / b[i]);
                    w -= b * (w[i] / b[i]);
                    insert(w);
                    insert(b * (m / b[i]));
                }
            }
        }
    }

    const bool full() const {
        for(auto i = 0; i < n; ++i) {
            if (basis[i][i] != 1) {
                return false;
            }
        }
        return true;
    }

    // Estimate the proportion of possible vectors that can be represented with the basis.
    const long double fullness() const { 
        long double ans = 1;
        for(int i = 0; i < n; ++i) {
            if (basis[i][i] != 0) {
                ans /= basis[i][i];
            } else {
                ans /= m;
            }
        }
        return ans;
    }
};

template<typename IntType>
IntType getRand(IntType l, IntType r, std::mt19937_64& gen) {
    std::uniform_int_distribution<> dist(l, r);
    return dist(gen);
}

std::shared_ptr<Operation> generateLinearMBA(std::shared_ptr<Node> node, Expression* expr) {
    const auto n = expr->n_variables;
    const auto nMask = (size_t)1 << n;

    // Create Nodes corresponding to each variable in the Expression
    std::vector<std::shared_ptr<Node>> variableNodes;
    for (auto i = 0; i < expr->n_variables; ++i) {
        std::shared_ptr<Node> varNode = std::make_shared<Node>(i, expr);
        variableNodes.push_back(varNode);
    }
    
    // Set up random number generation for later
    std::random_device rd;
    std::mt19937_64 gen(rd());

    Basis basis(nMask);
    std::cerr << std::endl;
    std::vector<std::shared_ptr<Node>> newNodes;
    while (!basis.full()) {
        auto nodes = variableNodes;

        // Repeat until we have only one Node left
        while(nodes.size() > 1) {
            // Select random nodes from the nodes list
            auto leftNodeIdx = getRand<size_t>(0, nodes.size() - 1, gen);
            auto leftNode = nodes[leftNodeIdx];
            std::swap(nodes[leftNodeIdx], nodes.back());
            nodes.pop_back();
            auto rightNodeIdx = getRand<size_t>(0, nodes.size() - 1, gen);
            auto rightNode = nodes[rightNodeIdx];
            std::swap(nodes[rightNodeIdx], nodes.back());
            nodes.pop_back();

            // Check if left node is a NOT operation
            bool leftNodeIsNOT = false;
            if (std::holds_alternative<std::shared_ptr<Operation>>(leftNode->value)) {
                auto operation = std::get<std::shared_ptr<Operation>>(leftNode->value);
                if (operation->type == OperationType::NOT) {
                    leftNodeIsNOT = true;
                }
            }

            // Generate an operation type (cannot be NOT if the left node is NOT already)
            OperationType opType = static_cast<OperationType>(getRand(0, 3 - leftNodeIsNOT, gen));
            if (opType == OperationType::NOT) {
                auto notLeftNode = std::make_shared<Node>(opType, leftNode, nullptr, expr);
                nodes.push_back(notLeftNode);
                nodes.push_back(rightNode);
            } else {
                auto newNode = std::make_shared<Node>(opType, leftNode, rightNode, expr);
                nodes.push_back(newNode);
            }
        }

        auto finalNode = nodes[0];

        // Insert vector inside basis if need be
        Vector finalNodeBitVector(*(finalNode->bitVector));
        if (!basis.insideBasis(finalNodeBitVector)) {
            std::cerr << "new Node: " << std::endl;
            Expression::printNode(finalNode, expr);
            std::cerr << std::endl;
            std::cerr << "bit vector: " << std::endl;
            for (auto i = 0; i < nMask; i++) std::cerr << finalNodeBitVector[i] << " ";
            std::cerr << std::endl;
            basis.insert(finalNodeBitVector);
            newNodes.push_back(finalNode);
        }
    }

    Matrix A(nMask, newNodes.size());
    for (auto i = 0; i < nMask; ++i) {
        for (auto j = 0; j < newNodes.size(); ++j) {
            A[i][j] = (*(newNodes[j]->bitVector))[i];
        }
    }
    Vector b(*(node->bitVector));

    Vector x = solve_linear_system(A, b);

    std::cerr << "A:" << std::endl;
    for (auto i = 0; i < nMask; ++i) {
        for (auto j = 0; j < newNodes.size(); ++j) {
            std::cerr << A[i][j] << ' ';
        }
        std::cerr << std::endl;
    }
    std::cerr << "b:" << std::endl;
    for (auto i = 0; i < nMask; ++i) {
        std::cerr << b[i] << " ";
    }
    std::cerr << std::endl;
    for (auto i = 0; i < newNodes.size(); i++) {
        std::cerr << x[i] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Done calculations" << std::endl;

    auto createNewOperation = [&](size_t l, size_t r, const auto& self) -> std::shared_ptr<Operation> {
        std::cerr << "merging " << l << ' ' << r << std::endl;

        // Base case: single element
        if (l == r) {
            u64 coeff = x[l];
            if (coeff == 0) {
                return nullptr;
            } else if (coeff == 1) {
                return std::get<std::shared_ptr<Operation>>(newNodes[l]->value);
            } else if (coeff == (u64)-1) {
                return std::make_shared<Operation>(OperationType::NEGATIVE, newNodes[l], nullptr);
            } else {
                auto scalarNode = std::make_shared<Node>(coeff, expr);
                return std::make_shared<Operation>(OperationType::MULTIPLY, scalarNode, newNodes[l]);
            }
        }

        // Recursive case: split into two halves
        auto mid = (l + r) / 2;
        auto leftOp = self(l, mid, self);
        auto rightOp = self(mid + 1, r, self);

        if (!leftOp && !rightOp) {
            return nullptr;
        } else if (!leftOp) {
            return rightOp;
        } else if (!rightOp) {
            return leftOp;
        } else {
            // If both sides are not empty
            // Check if right operation is NEGATIVE, if so, turn ADD into SUBTRACT
            if (rightOp->type == OperationType::NEGATIVE) {
                const auto leftNode = std::make_shared<Node>(leftOp, expr);
                return std::make_shared<Operation>(OperationType::SUBTRACT, leftNode, rightOp->left);
            }

            // Combine with an ADD operation if no special case
            const auto leftNode = std::make_shared<Node>(leftOp, expr);
            const auto rightNode = std::make_shared<Node>(rightOp, expr);
            return std::make_shared<Operation>(OperationType::ADD, leftNode, rightNode);
        }
    };

    std::cerr << "Creating Node" << std::endl;
    auto newLinearMBA = createNewOperation(0, newNodes.size() - 1, createNewOperation);

    // Testing
    auto newLinearMBANode = std::make_shared<Node>(newLinearMBA, expr);
    Expression::printNode(newLinearMBANode, expr);

    return newLinearMBA;
}

/*
Generates a Semi-Linear MBA.
The ground truth value is the main value of the returned Node.
The obfuscated value is pushed as an alternate representation.
*/ 
std::shared_ptr<Node> generateSemiLinearMBA(std::shared_ptr<Node> originalNode, Expression* expr, int alternations = 2) {
    const auto n = expr->n_variables;
    const auto nMask = (size_t)1 << n;
    const auto nAlts = (size_t)1 << alternations;
    const auto basisSize = nMask * nAlts;

    // Create Nodes corresponding to each variable in the Expression
    std::vector<std::shared_ptr<Node>> variableNodes;
    for (auto i = 0; i < expr->n_variables; ++i) {
        std::shared_ptr<Node> varNode = std::make_shared<Node>(i, expr);
        variableNodes.push_back(varNode);
    }

    // Create Nodes for zero and (negative) one
    auto zeroNode = std::make_shared<Node>((u64)0, expr);
    auto oneNode = std::make_shared<Node>((u64)-1, expr);
    
    // Set up random number generation for later
    std::random_device rd;
    std::mt19937_64 gen(rd());

    Basis basis(basisSize);
    std::vector<std::pair<std::vector<std::pair<std::shared_ptr<Node>, uint32_t>>, Vector>> newNodeList;
    int cnt = 0;
    while (!basis.full()) {
        auto nodes = variableNodes;
        int alternationsRemaining = alternations;

        // Repeat until we have only one Node left
        while (nodes.size() > 1) {
            // Select random nodes from the nodes list
            auto leftNodeIdx = getRand<size_t>(0, nodes.size() - 1, gen);
            auto leftNode = nodes[leftNodeIdx];
            std::swap(nodes[leftNodeIdx], nodes.back());
            nodes.pop_back();
            auto rightNodeIdx = getRand<size_t>(0, nodes.size() - 1, gen);
            auto rightNode = nodes[rightNodeIdx];
            std::swap(nodes[rightNodeIdx], nodes.back());
            nodes.pop_back();

            // Check if left node is a NOT operation
            bool leftNodeIsNOT = false;
            if (std::holds_alternative<std::shared_ptr<Operation>>(leftNode->value)) {
                auto operation = std::get<std::shared_ptr<Operation>>(leftNode->value);
                if (operation->type == OperationType::NOT) {
                    leftNodeIsNOT = true;
                }
            }

            // Generate an operation type (cannot be NOT if the left node is NOT already)
            OperationType opType = static_cast<OperationType>(getRand(0, 3 - leftNodeIsNOT, gen));
            if (opType == OperationType::NOT) {
                auto notLeftNode = std::make_shared<Node>(opType, leftNode, nullptr, expr);
                nodes.push_back(notLeftNode);
                nodes.push_back(rightNode);
            } else {
                auto newNode = std::make_shared<Node>(opType, leftNode, rightNode, expr);
                nodes.push_back(newNode);
            }
        }

        /*
        Expression::printNode(nodes[0]);
        std::cerr << std::endl;
        */
       
        // Apply alternations recursively
        auto applyAlternations = [&](std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> nodes, int& alternationsLeft, const auto& self, bool isRoot = true) -> std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> {
            static int depth = -1;
            depth++;

            /*
            std::cerr << "applying alternations " << depth << std::endl;
            Expression::printNode(nodes[0].first);
            std::cerr << std::endl;
            */

            if (alternationsLeft == 0 || !std::holds_alternative<std::shared_ptr<Operation>>(nodes[0].first->value)) {
                depth--;
                return nodes;
            }

            // Gather left and right children for all nodes in the vector
            std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> leftChildren;
            std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> rightChildren;
            for (const auto& node : nodes) {
                auto operation = std::get<std::shared_ptr<Operation>>(node.first->value);
                if (operation->left) {
                    leftChildren.push_back({operation->left, node.second});
                }
                if (operation->right) {
                    rightChildren.push_back({operation->right, node.second});
                }
            }

            // Recursively descend down the left and right sides of the tree
            std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> newLeftChildren, newRightChildren;
            if (!leftChildren.empty()) newLeftChildren = self(leftChildren, alternationsLeft, self, false);
            if (!rightChildren.empty()) newRightChildren = self(rightChildren, alternationsLeft, self, false);

            // Create a new vector of nodes by combining new left and right children
            std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> combinedNodes;
            auto operation = std::get<std::shared_ptr<Operation>>(nodes[0].first->value);
            for (const auto& leftChild : newLeftChildren) {
                std::shared_ptr<Node> newNode;

                if (operation->type == OperationType::NOT) {
                    // Only use the left child for NOT operation
                    newNode = std::make_shared<Node>(OperationType::NOT, leftChild.first, nullptr, expr);
                    combinedNodes.push_back({newNode, leftChild.second});
                } else {
                    for (const auto& rightChild : newRightChildren) {
                        newNode = std::make_shared<Node>(operation->type, leftChild.first, rightChild.first, expr);
                        combinedNodes.push_back({newNode, leftChild.second | rightChild.second});
                    }
                }
            }

            // Apply an alternation with a 1/n chance
            const int tries = isRoot ? alternationsLeft : 2;
            for (int iter = 1; iter <= tries; ++iter) {
                if (alternationsLeft > 0 && (isRoot || getRand<size_t>(0, n - 1, gen) == 0)) {
                    /*
                    std::cerr << "Split ";
                    Expression::printNode(nodes[0].first);
                    std::cerr << std::endl;
                    std::cerr << "Previous combined nodes size: " << combinedNodes.size() << std::endl;
                    */

                    alternationsLeft--;
                    OperationType altType = static_cast<OperationType>(getRand(0, 2, gen)); // AND, OR, XOR
                    if (getRand(0, 4, gen)) altType = OperationType::XOR; // prefer XOR
                    std::vector<std::pair<std::shared_ptr<Node>, uint32_t>> alternatedNodes;

                    // Duplicate each node and apply the alternation
                    for (const auto& node: combinedNodes) {
                        std::shared_ptr<Node> alternationNode;
                        uint32_t alternationBitmask = node.second | (1u << (alternations - alternationsLeft - 1));

                        switch (altType) {
                            case OperationType::AND:
                                alternationNode = std::make_shared<Node>(OperationType::AND, node.first, zeroNode, expr);
                                break;
                            case OperationType::OR:
                                alternationNode = std::make_shared<Node>(OperationType::OR, node.first, oneNode, expr);
                                break;
                            case OperationType::XOR:
                                alternationNode = std::make_shared<Node>(OperationType::XOR, node.first, oneNode, expr);
                                break;
                            default:
                                break;
                        }
                        alternatedNodes.push_back({alternationNode, alternationBitmask});
                    }

                    // Insert the alternated nodes into the combined nodes
                    combinedNodes.insert(combinedNodes.end(), alternatedNodes.begin(), alternatedNodes.end());
                    // std::cerr << "combined nodes size: " << combinedNodes.size() << std::endl;
                }
            }

            depth--;
            return combinedNodes;
        };
        
        // Apply alternations to the final node
        auto newNodes = applyAlternations({std::make_pair(nodes[0], 0u)}, alternationsRemaining, applyAlternations);
        std::shuffle(newNodes.begin(), newNodes.end(), gen);

        /*
        for (auto& [node, mask]: newNodes) {
            Expression::printNode(node); 
            std::cerr << " " << mask << std::endl;
        }
        std::cerr << "nodes size: " << newNodes.size() << std::endl;
        */

        // Generate the vector of differences between alternations
        Vector vec(basisSize);
        for (auto i = 0; i < nMask; ++i) {
            vec[i] = (*(newNodes[0].first->bitVector))[i];
        }
        for (auto i = 1; i < nAlts; ++i) {
            const auto& firstBitVector = *(newNodes[0].first->bitVector);
            const auto& curBitVector = *(newNodes[i].first->bitVector);
            for (auto j = 0; j < nMask; ++j) {
                vec[i * nMask + j] = firstBitVector[j] - curBitVector[j] ;
            }
        }

        if (!basis.insideBasis(vec)) {
            basis.insert(vec);
            for (const auto& [node, mask]: newNodes) {
                Expression::printNode(node, expr);
                std::cerr << std::endl;
            }
            std::cerr << "Acquired vector: ";
            vec.print();
            std::cerr << " at iteration " << cnt + 1 << std::endl;
            newNodeList.push_back(std::make_pair(newNodes, vec));
        }

        cnt++;
        if (cnt > 10000) break;
    }
    std::cerr << "basis fullness: " << basis.fullness() << std::endl;
    for (auto i = 0; i < basis.basis.size(); ++i) {
        std::cerr << basis.basis[i][i] << " ";
    }
    std::cerr << std::endl;

    Matrix A(basisSize, newNodeList.size());
    for (auto i = 0; i < newNodeList.size(); ++i) {
        const auto& vec = newNodeList[i].second;
        for (auto j = 0; j < basisSize; ++j) {
            A[j][i] = vec[j];
        }
    }
    Vector b(basisSize);
    for (auto i = 0; i < nMask; ++i) {
        b[i] = (*(originalNode->bitVector))[i];
    }

    std::vector<u64> compensation(nAlts, 0);
    for (auto i = 1; i < nAlts; ++i) {
        compensation[i] = getRand<u64>(0, ULLONG_MAX, gen) * 2 + 1;
        for (auto j = 0; j < nMask; ++j) {
            b[i * nMask + j] = compensation[i];
        }
    }

    /* std::cerr << "A:" << std::endl;
    for (auto i = 0; i < A.size(); ++i) {
        for (auto j = 0; j < A[i].size(); ++j) {
            std::cerr << A[i][j] << ' ';
        }
        std::cerr << std::endl;
    } */ 
    std::cerr << "b:" << std::endl;
    for (auto i = 0; i < b.size(); ++i) {
        std::cerr << b[i] << " ";
    }
    std::cerr << std::endl;

    Vector x = solve_linear_system(A, b);

    std::cerr << "x:" << std::endl;
    for (auto i = 0; i < x.size(); i++) {
        std::cerr << x[i] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Done calculations" << std::endl;

    // Decide which variant to use
    int var[64];
    u64 totalCompensation = 0;
    for (auto i = 0; i < 64; ++i) {
        var[i] = getRand<size_t>(0, nAlts - 1, gen);
        totalCompensation += compensation[var[i]] * ((u64)1 << i);
    }

    // create final Node
    std::shared_ptr<Node> totalNode = std::make_shared<Node>(totalCompensation, expr);

    std::cerr << "Variants chosen: ";
    for (auto i = 0; i < 64; ++i) {
        std::cerr << var[i] << " ";
    }
    std::cerr << std::endl;

    // Generate new Nodes and add them up
    std::vector<std::shared_ptr<Node>> semiLinearNodes;
    for (auto i = 0; i < newNodeList.size(); ++i) {
        if (x[i] == 0) continue;

        const auto& [nodes, _] = newNodeList[i];

        u64 val[alternations] = {0};
        for (int bit = 0; bit < 64; ++bit) {
            int chosen = var[bit];
            auto msk = nodes[chosen].second;
            for (auto i = 0; i < alternations; ++i) {
                if ((msk >> i) & 1) {
                    val[i] |= ((u64)1 << bit);
                }
            }
        }

        int mainNodeIdx = 0;
        for (auto i = 0; i < nAlts; i++) {
            if (nodes[i].second == nAlts - 1) mainNodeIdx = i;
        }
        const std::shared_ptr<Node>& mainNode = nodes[mainNodeIdx].first;

        /*
        std::cerr << "mask\n";
        for (auto i = 0; i < nAlts; ++i) {
            std::cerr << nodes[i].second << ' ';
        }
        std::cerr << std::endl;
        std::cerr << "vals\n";
        for (auto i = 0; i < alternations; ++i) {
            std::cerr << val[i] << ' ';
        }
        std::cerr << std::endl;
        */

        int alt = 0;
        // I have no idea why this must be wrapped in an std::function, but so be it
        std::function<std::shared_ptr<Node>(std::shared_ptr<Node>)> assignSemiLinearMBA = 
        [&](std::shared_ptr<Node> node) -> std::shared_ptr<Node> {
            // If the node is a scalar value
            if (std::holds_alternative<u64>(node->value)) {
                u64 scalarValue = std::get<u64>(node->value);
                // std::cerr << "alt " << alt << std::endl;
                u64 MBAValue = scalarValue == 0 ? ~val[alt] : val[alt];
                auto scalarNode = std::make_shared<Node>(MBAValue, expr);
                alt++;
                return scalarNode;
            }

            // If the node is an operation
            else if (std::holds_alternative<std::shared_ptr<Operation>>(node->value)) {
                auto operation = std::get<std::shared_ptr<Operation>>(node->value);

                // Special handling for NOT operation
                if (operation->type == OperationType::NOT) {
                    std::shared_ptr<Node> newLeft = operation->left ? assignSemiLinearMBA(operation->left) : nullptr;
                    return std::make_shared<Node>(operation->type, newLeft, nullptr, expr);
                }

                // Recursively assign MBA for left and right children for other operations
                std::shared_ptr<Node> newLeft = operation->left ? assignSemiLinearMBA(operation->left) : nullptr;
                std::shared_ptr<Node> newRight = operation->right ? assignSemiLinearMBA(operation->right) : nullptr;

                // Reconstruct the node with the new left and right children
                return std::make_shared<Node>(operation->type, newLeft, newRight, expr);
            }

            // In case of unknown types, return the original node
            return node;
        };

        auto finalNode = std::make_shared<Node>(
            OperationType::MULTIPLY, 
            std::make_shared<Node>(x[i], expr), 
            assignSemiLinearMBA(mainNode), 
            expr
        );
        totalNode = std::make_shared<Node>(OperationType::ADD, totalNode, finalNode, expr);
        // Expression::printNode(finalNode);
        // std::cerr << std::endl;
    }

    Expression::printNode(totalNode, expr);
    return totalNode;
}

int main() {
    // Create a mock Expression
    Expression expr;
    expr.n_variables = 3;
    expr.variables = {"x", "y", "z"};
    expr.variableMap = { {"x", 0}, {"y", 1}, {"z", 2} };

    auto printVector = []<typename T>(const std::vector<T>& vec) -> void {
        for (const auto& elem: vec) std::cerr << elem << " ";
        std::cerr << std::endl;
    };

    // Create a root node (which will be updated by generateLinearMBA)
    std::shared_ptr<Node> xNode = std::make_shared<Node>((int)0, &expr);
    std::shared_ptr<Node> yNode = std::make_shared<Node>((int)1, &expr);
    std::shared_ptr<Node> zNode = std::make_shared<Node>((int)2, &expr);

    auto twoNode = std::make_shared<Node>((u64)2, &expr);
    auto xyNode = std::make_shared<Node>(OperationType::XOR, xNode, yNode, &expr);
    auto xyzNode = std::make_shared<Node>(OperationType::OR, zNode, xyNode, &expr);
    auto mulzNode = std::make_shared<Node>(OperationType::MULTIPLY, twoNode, zNode, &expr);
    auto rootNode = std::make_shared<Node>(OperationType::SUBTRACT, xyzNode, mulzNode, &expr);
    expr.root = rootNode;

    printVector(*(zNode->bitVector));
    printVector(*(twoNode->bitVector));
    printVector(*(xyzNode->bitVector));
    printVector(*(mulzNode->bitVector));
    printVector(*(rootNode->bitVector));

    // Generate a random boolean-arithmetic expression
    // generateLinearMBA(rootNode, &expr);

    // Output the generated root node (for basic testing, you'd want to implement an actual traversal or printing logic)
    // std::cout << "Linear MBA generation complete.\n";

    // generateSemiLinearMBA(rootNode, &expr, 1);
    generateSemiLinearMBA(rootNode, &expr, 2);
    
    // Additional assertions or checks can be placed here to validate properties of the generated expression.
}

