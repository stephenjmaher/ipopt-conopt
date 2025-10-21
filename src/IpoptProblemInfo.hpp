/**
 * @file IpoptProblemInfo.hpp
 * @brief Struct to store problem information retrieved from Ipopt TNLP interface.
 */

#ifndef IPOPT_PROBLEM_INFO_HPP
#define IPOPT_PROBLEM_INFO_HPP

#include <vector>
#include <string>
#include "IpoptTypes.hpp"

namespace Ipopt {

    /**
     * @brief Enumeration for variable/constraint linearity types
     */
    enum LinearityType {
        LINEAR = 0,
        NONLINEAR = 1
    };

    /**
     * @brief Enumeration for index style (C vs Fortran)
     */
    enum IndexStyleEnum {
        C_STYLE = 0,      // 0-based indexing
        FORTRAN_STYLE = 1 // 1-based indexing
    };

    /**
     * @brief Enumeration for the constraint types.
     *
     * NOTE: the RANGE type is not used by CONOPT. This is used to help process the constraints from Ipopt.
     */
    enum ConoptConstraintType {
        EQUAL = 0,
        GREATEREQ = 1,
        LESSEQ = 2,
        FREE = 3,
        RANGE = 4
    };

    /**
     * @brief Struct to hold all problem information retrieved from TNLP
     * This information is gathered before solving and used by CONOPT callbacks
     */
    struct IpoptProblemInfo {

        // === Problem Dimensions ===
        Index n;                    // Number of variables
        Index m;                    // Number of constraints
        Index nnz_jac_g;           // Number of non-zeros in constraint Jacobian
        Index nnz_h_lag;           // Number of non-zeros in Hessian of Lagrangian
        IndexStyleEnum index_style; // C-style (0-based) or Fortran-style (1-based)

        // === Variable Information ===
        std::vector<Number> x_l;           // Variable lower bounds (length n)
        std::vector<Number> x_u;           // Variable upper bounds (length n)
        std::vector<Number> x_init;        // Initial variable values (length n)
        std::vector<Number> z_L_init;      // Initial lower bound multipliers (length n)
        std::vector<Number> z_U_init;      // Initial upper bound multipliers (length n)
        std::vector<LinearityType> var_linearity; // Variable linearity types (length n)
        Index num_nonlin_vars;             // Number of nonlinear variables

        // === Constraint Information ===
        std::vector<Number> g_l;           // Constraint lower bounds (length m)
        std::vector<Number> g_u;           // Constraint upper bounds (length m)
        std::vector<Number> lambda_init;   // Initial constraint multipliers (length m)
        std::vector<LinearityType> const_linearity; // Constraint linearity types (length m)

        // === Split Constraint Information (for CONOPT) ===
        Index m_split;                     // Number of constraints after splitting
        std::vector<Number> g_rhs;         // RHS values for split constraints (length m_split)
        std::vector<Number> g_type;        // the constraint type for Conopt
        std::vector<LinearityType> const_linearity_split; // Linearity for split constraints (length m_split)
        std::vector<Index> original_constraint_map; // Maps split constraint index to original constraint index (length m_split)

        // === Jacobian Structure ===
        std::vector<Index> jac_g_iRow;     // Jacobian row indices (length nnz_jac_g)
        std::vector<Index> jac_g_jCol;     // Jacobian column indices (length nnz_jac_g)
        std::vector<Number> jac_g_values;  // Jacobian values (length nnz_jac_g)

        // === Split Jacobian Structure (for CONOPT) ===
        Index nnz_jac_g_split;             // Number of non-zeros in split Jacobian
        std::vector<Index> jac_g_iRow_split; // Split Jacobian row indices (length nnz_jac_g_split)
        std::vector<Index> jac_g_jCol_split; // Split Jacobian column indices (length nnz_jac_g_split)
        std::vector<Number> jac_g_values_split; // Split Jacobian values (length nnz_jac_g_split)

        // === Hessian Structure ===
        std::vector<Index> hess_iRow;      // Hessian row indices (length nnz_h_lag)
        std::vector<Index> hess_jCol;      // Hessian column indices (length nnz_h_lag)
        std::vector<Number> hess_values;   // Hessian values (length nnz_h_lag)

        // === Problem Metadata ===
        bool has_initial_x;                // Whether initial x values were provided
        bool has_initial_z;                // Whether initial z values were provided
        bool has_initial_lambda;           // Whether initial lambda values were provided
        bool has_variable_linearity;       // Whether variable linearity info was provided
        bool has_constraint_linearity;     // Whether constraint linearity info was provided
        bool has_nonlinear_vars;           // Whether nonlinear variable count was provided

        // === Constructor ===
        IpoptProblemInfo() :
            n(0), m(0), nnz_jac_g(0), nnz_h_lag(0), index_style(C_STYLE),
            m_split(0), nnz_jac_g_split(0),
            num_nonlin_vars(0),
            has_initial_x(false), has_initial_z(false), has_initial_lambda(false),
            has_variable_linearity(false), has_constraint_linearity(false),
            has_nonlinear_vars(false)
        {}

        // === Utility Methods ===

        /**
         * @brief Resize all vectors based on problem dimensions
         */
        void resize_vectors() {
            x_l.resize(n);
            x_u.resize(n);
            x_init.resize(n);
            z_L_init.resize(n);
            z_U_init.resize(n);
            var_linearity.resize(n);

            g_l.resize(m);
            g_u.resize(m);
            lambda_init.resize(m);
            const_linearity.resize(m);

            jac_g_iRow.resize(nnz_jac_g);
            jac_g_jCol.resize(nnz_jac_g);
            jac_g_values.resize(nnz_jac_g);

            hess_iRow.resize(nnz_h_lag);
            hess_jCol.resize(nnz_h_lag);
            hess_values.resize(nnz_h_lag);
        }

        /**
         * @brief utility method for getting the constraint type.
         */
        ConoptConstraintType constraint_type(int orig_row) {
            ConoptConstraintType type;
            bool has_lower = IsFiniteNumber(g_l[orig_row]);
            bool has_upper = IsFiniteNumber(g_u[orig_row]);

            if (has_lower && has_upper) {
                if (g_l[orig_row] == g_u[orig_row])
                    type = ConoptConstraintType::EQUAL;
                else
                    type = ConoptConstraintType::RANGE;
            }
            else if (has_upper)
                type = ConoptConstraintType::LESSEQ;
            else if (has_lower)
                type = ConoptConstraintType::GREATEREQ;
            else
                type = ConoptConstraintType::FREE;

            return type;
        }


        /**
         * @brief Split constraints that have both lower and upper bounds
         * This method processes the original constraints and creates split constraints
         * for CONOPT, where constraints with both bounds are duplicated and negated.
         */
        void split_constraints() {
            // Clear previous split data
            m_split = 0;
            nnz_jac_g_split = 0;
            g_rhs.clear();
            g_type.clear();
            const_linearity_split.clear();
            original_constraint_map.clear();
            jac_g_iRow_split.clear();
            jac_g_jCol_split.clear();
            jac_g_values_split.clear();

            // Count how many constraints we'll have after splitting
            for (Index i = 0; i < m; ++i) {
                ConoptConstraintType type = constraint_type(i);
                if (type == ConoptConstraintType::RANGE) {
                    // Both bounds: need two constraints
                    m_split += 2;
                } else {
                    // Single constraint: equality, lower bound only, upper bound only, or free
                    m_split += 1;
                }
            }

            // Resize split constraint vectors
            g_rhs.resize(m_split);
            g_type.resize(m_split);
            const_linearity_split.resize(m_split);
            original_constraint_map.resize(m_split);

            // Process each original constraint
            Index split_idx = 0;
            for (Index i = 0; i < m; ++i) {
                ConoptConstraintType type = constraint_type(i);

                if (type == ConoptConstraintType::EQUAL) {
                    g_rhs[split_idx] = g_l[i];
                    g_type[split_idx] = ConoptConstraintType::EQUAL;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;
                }
                else if (type == ConoptConstraintType::RANGE) {
                    // Both bounds: create two constraints
                    // First constraint: g(x) >= g_l
                    g_rhs[split_idx] = g_l[i];
                    g_type[split_idx] = ConoptConstraintType::GREATEREQ;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;

                    // Second constraint: g(x) <= g_u
                    g_rhs[split_idx] = g_u[i];
                    g_type[split_idx] = ConoptConstraintType::LESSEQ;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;
                }
                else if (type == ConoptConstraintType::GREATEREQ) {
                    // Only lower bound: g(x) >= g_l
                    g_rhs[split_idx] = g_l[i];
                    g_type[split_idx] = ConoptConstraintType::GREATEREQ;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;
                }
                else if (type == ConoptConstraintType::LESSEQ) {
                    // Only upper bound: g(x) <= g_u
                    g_rhs[split_idx] = g_u[i];
                    g_type[split_idx] = ConoptConstraintType::LESSEQ;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;
                }
                else if (type == ConoptConstraintType::FREE) {
                    // Free constraint: no bounds
                    g_rhs[split_idx] = 0.0; // RHS doesn't matter for free constraints
                    g_type[split_idx] = ConoptConstraintType::FREE;
                    const_linearity_split[split_idx] = const_linearity[i];
                    original_constraint_map[split_idx] = i;
                    split_idx++;
                }
            }

            // Now split the Jacobian structure
            split_jacobian_structure();
        }

        /**
         * @brief Split the Jacobian structure to match split constraints
         */
        void split_jacobian_structure() {
            nnz_jac_g_split = 0;
            jac_g_iRow_split.clear();
            jac_g_jCol_split.clear();
            jac_g_values_split.clear();

            // Count non-zeros for split Jacobian
            for (Index k = 0; k < nnz_jac_g; ++k) {
                Index orig_row = jac_g_iRow[k];
                ConoptConstraintType type = constraint_type(orig_row);

                if (type == ConoptConstraintType::RANGE) {
                    // Both bounds: this entry appears in both split constraints
                    nnz_jac_g_split += 2;
                } else if (type != ConoptConstraintType::FREE) {
                    // All other types: this entry appears once
                    nnz_jac_g_split += 1;
                }
            }

            // Resize split Jacobian vectors
            jac_g_iRow_split.resize(nnz_jac_g_split);
            jac_g_jCol_split.resize(nnz_jac_g_split);
            jac_g_values_split.resize(nnz_jac_g_split);

            // Create split Jacobian structure
            Index split_k = 0;
            for (Index k = 0; k < nnz_jac_g; ++k) {
                Index orig_row = jac_g_iRow[k];
                Index col = jac_g_jCol[k];
                ConoptConstraintType type = constraint_type(orig_row);

                if (type == ConoptConstraintType::RANGE) {
                    // Both bounds: create two entries
                    // Find the split constraint indices for this original constraint
                    Index split_row_lower = -1, split_row_upper = -1;
                    for (Index s = 0; s < m_split; ++s) {
                        if (original_constraint_map[s] == orig_row) {
                            if (g_type[s] == ConoptConstraintType::GREATEREQ) {
                                split_row_lower = s;
                            } else if (g_type[s] == ConoptConstraintType::LESSEQ) {
                                split_row_upper = s;
                            }
                        }
                    }

                    // Lower bound constraint
                    jac_g_iRow_split[split_k] = split_row_lower;
                    jac_g_jCol_split[split_k] = col;
                    jac_g_values_split[split_k] = 0.0; // Will be set during evaluation
                    split_k++;

                    // Upper bound constraint
                    jac_g_iRow_split[split_k] = split_row_upper;
                    jac_g_jCol_split[split_k] = col;
                    jac_g_values_split[split_k] = 0.0; // Will be set during evaluation
                    split_k++;

                } else if (type != ConoptConstraintType::FREE) {
                    // All other types: find the corresponding split constraint
                    Index split_row = -1;
                    for (Index s = 0; s < m_split; ++s) {
                        if (original_constraint_map[s] == orig_row) {
                            split_row = s;
                            break;
                        }
                    }
                    jac_g_iRow_split[split_k] = split_row;
                    jac_g_jCol_split[split_k] = col;
                    jac_g_values_split[split_k] = 0.0; // Will be set during evaluation
                    split_k++;
                }
            }
        }

        /**
         * @brief Check if all required information is available
         */
        bool is_complete() const {
            return (n > 0 && m >= 0 &&
                    x_l.size() == n && x_u.size() == n &&
                    g_l.size() == m && g_u.size() == m &&
                    jac_g_iRow.size() == nnz_jac_g &&
                    jac_g_jCol.size() == nnz_jac_g);
        }

        /**
         * @brief Clear all data
         */
        void clear() {
            n = m = nnz_jac_g = nnz_h_lag = 0;
            m_split = nnz_jac_g_split = 0;
            index_style = C_STYLE;
            num_nonlin_vars = 0;

            x_l.clear();
            x_u.clear();
            x_init.clear();
            z_L_init.clear();
            z_U_init.clear();
            var_linearity.clear();

            g_l.clear();
            g_u.clear();
            lambda_init.clear();
            const_linearity.clear();

            jac_g_iRow.clear();
            jac_g_jCol.clear();
            jac_g_values.clear();

            hess_iRow.clear();
            hess_jCol.clear();
            hess_values.clear();

            // Clear split constraint data
            g_rhs.clear();
            g_type.clear();
            const_linearity_split.clear();
            original_constraint_map.clear();
            jac_g_iRow_split.clear();
            jac_g_jCol_split.clear();
            jac_g_values_split.clear();

            has_initial_x = has_initial_z = has_initial_lambda = false;
            has_variable_linearity = has_constraint_linearity = false;
            has_nonlinear_vars = false;
        }

        /**
         * @brief Get a string representation for debugging
         */
        std::string to_string() const {
            std::string result = "IpoptProblemInfo:\n";
            result += "  Dimensions: n=" + std::to_string(n) +
                     ", m=" + std::to_string(m) +
                     ", nnz_jac_g=" + std::to_string(nnz_jac_g) +
                     ", nnz_h_lag=" + std::to_string(nnz_h_lag) + "\n";
            result += "  Split dimensions: m_split=" + std::to_string(m_split) +
                     ", nnz_jac_g_split=" + std::to_string(nnz_jac_g_split) + "\n";
            result += "  Index style: " + std::string(index_style == C_STYLE ? "C_STYLE" : "FORTRAN_STYLE") + "\n";
            result += "  Has initial x: " + std::string(has_initial_x ? "yes" : "no") + "\n";
            result += "  Has initial z: " + std::string(has_initial_z ? "yes" : "no") + "\n";
            result += "  Has initial lambda: " + std::string(has_initial_lambda ? "yes" : "no") + "\n";
            result += "  Complete: " + std::string(is_complete() ? "yes" : "no") + "\n";
            return result;
        }
    };

} // namespace Ipopt

#endif // IPOPT_PROBLEM_INFO_HPP
