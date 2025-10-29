/**
 * @file IpoptProblemInfo.hpp
 * @brief Struct to store problem information retrieved from Ipopt TNLP interface.
 */

#ifndef IPOPT_PROBLEM_INFO_HPP
#define IPOPT_PROBLEM_INFO_HPP

// CONOPT infinity value - used to represent unbounded constraints
#define CONOPT_INFINITY 1e12

#include <vector>
#include <string>
#include <algorithm>
#include "IpoptTypes.hpp"

namespace Ipopt {

/**
 * @brief Enumeration for variable/constraint linearity types
 */
enum LinearityType { LINEAR = 0, NONLINEAR = 1 };

/**
 * @brief Enumeration for index style (C vs Fortran)
 */
enum IndexStyleEnum {
   C_STYLE = 0,       // 0-based indexing
   FORTRAN_STYLE = 1  // 1-based indexing
};

/**
 * @brief Enumeration for the constraint types.
 *
 * NOTE: the RANGE type is not used by CONOPT. This is used to help process the constraints from
 * Ipopt.
 */
enum ConoptConstraintType { EQUAL = 0, GREATEREQ = 1, LESSEQ = 2, FREE = 3, RANGE = 4 };

/**
 * @brief Struct to hold all problem information retrieved from TNLP
 * This information is gathered before solving and used by CONOPT callbacks
 */
struct IpoptProblemInfo {

   // === Problem Dimensions ===
   Index n;                     // Number of variables
   Index m;                     // Number of constraints
   Index nnz_jac_g;             // Number of non-zeros in constraint Jacobian
   Index nnz_h_lag;             // Number of non-zeros in Hessian of Lagrangian
   IndexStyleEnum index_style;  // C-style (0-based) or Fortran-style (1-based)

   // === Variable Information ===
   std::vector<Number> x_l;                   // Variable lower bounds (length n)
   std::vector<Number> x_u;                   // Variable upper bounds (length n)
   std::vector<Number> x_init;                // Initial variable values (length n)
   std::vector<Number> z_L_init;              // Initial lower bound multipliers (length n)
   std::vector<Number> z_U_init;              // Initial upper bound multipliers (length n)
   std::vector<LinearityType> var_linearity;  // Variable linearity types (length n)
   Index num_nonlin_vars;                     // Number of nonlinear variables

   // === Constraint Information ===
   std::vector<Number> g_l;                     // Constraint lower bounds (length m)
   std::vector<Number> g_u;                     // Constraint upper bounds (length m)
   std::vector<Number> lambda_init;             // Initial constraint multipliers (length m)
   std::vector<LinearityType> const_linearity;  // Constraint linearity types (length m)

   // === Constraint Splitting Mapping (for CONOPT) ===
   Index m_split;  // Number of constraints after splitting (includes objective)
   std::vector<Index> original_constraint_map;  // Maps split constraint index to original
                                                // constraint index (length m_split)
   std::vector<Index> split_constraint_map;     // Maps original constraint index to first split
                                                // constraint index (length m)
   Index objective_row_index;  // Index of the objective row in the split constraint matrix (-1 if
                               // not added)

   // === Jacobian Structure ===
   std::vector<Index> jac_g_iRow;     // Jacobian row indices (length nnz_jac_g)
   std::vector<Index> jac_g_jCol;     // Jacobian column indices (length nnz_jac_g)
   std::vector<Number> jac_g_values;  // Jacobian values (length nnz_jac_g)

   // === Jacobian Splitting Mapping (for CONOPT) ===
   Index nnz_jac_g_split;                  // Number of non-zeros in split Jacobian
   std::vector<Index> jacobian_split_map;  // Maps split Jacobian index to original Jacobian index
                                           // (length nnz_jac_g_split)
   std::vector<Index>
         jacobian_split_rows;  // Direct mapping to split constraint rows (length nnz_jac_g_split)

   // === Hessian Structure ===
   std::vector<Index> hess_iRow;     // Hessian row indices (length nnz_h_lag)
   std::vector<Index> hess_jCol;     // Hessian column indices (length nnz_h_lag)
   std::vector<Number> hess_values;  // Hessian values (length nnz_h_lag)
   // Permutation mapping to enforce CONOPT ordering (sorted by column then row)
   // Maps sorted index -> original TNLP index
   std::vector<Index> hess_perm_sorted_to_orig;

   // === Problem Metadata ===
   bool init_x_req;                // Whether initial x values is required from the user
   bool init_z_req;                // Whether initial z values is required from the user
   bool init_lambda_req;           // Whether initial lambda values is required from the user
   bool has_variable_linearity;    // Whether variable linearity info was provided
   bool has_constraint_linearity;  // Whether constraint linearity info was provided
   bool has_nonlinear_vars;        // Whether nonlinear variable count was provided

   // === Constructor ===
   IpoptProblemInfo()
       : n(0), m(0), nnz_jac_g(0), nnz_h_lag(0), index_style(C_STYLE), m_split(0),
         objective_row_index(-1), nnz_jac_g_split(0), num_nonlin_vars(0), init_x_req(true),
         init_z_req(false), init_lambda_req(false), has_variable_linearity(false),
         has_constraint_linearity(false), has_nonlinear_vars(false) {}

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
    * @brief Compute and store permutation to enforce CONOPT-required Hessian ordering.
    * Sorts entries by column (non-decreasing) then row (increasing within same column).
    */
   void compute_hessian_permutation() {
      const Index nnz = nnz_h_lag;
      hess_perm_sorted_to_orig.clear();
      hess_perm_sorted_to_orig.reserve(nnz);
      for (Index k = 0; k < nnz; ++k)
         hess_perm_sorted_to_orig.push_back(k);
      std::stable_sort(hess_perm_sorted_to_orig.begin(), hess_perm_sorted_to_orig.end(),
            [&](Index a, Index b) {
               Index ca = hess_jCol[a];
               Index cb = hess_jCol[b];
               if (ca != cb)
                  return ca < cb;
               Index ra = hess_iRow[a];
               Index rb = hess_iRow[b];
               return ra < rb;
            });
   }

   /**
    * @brief utility method for getting the constraint type.
    */
   ConoptConstraintType constraint_type(int orig_row) const {
      ConoptConstraintType type;

      // Check if bounds represent infinity (IPOPT uses large finite numbers like 1e19, 2e19)
      bool has_lower = IsFiniteNumber(g_l[orig_row]) && g_l[orig_row] < CONOPT_INFINITY;
      bool has_upper = IsFiniteNumber(g_u[orig_row]) && g_u[orig_row] < CONOPT_INFINITY;

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
      objective_row_index = -1;
      original_constraint_map.clear();
      split_constraint_map.clear();
      jacobian_split_map.clear();

      // Initialize split constraint mapping for original constraints
      split_constraint_map.resize(m, -1);

      // Count how many constraints we'll have after splitting
      for (Index i = 0; i < m; ++i) {
         ConoptConstraintType type = constraint_type(i);
         if (type == ConoptConstraintType::RANGE) {
            // Both bounds: need two constraints
            m_split += 2;
         }
         else {
            // Single constraint: equality, lower bound only, upper bound only, or free
            m_split += 1;
         }
      }

      // Add one more constraint for the objective function
      m_split += 1;

      // Resize mapping vectors
      original_constraint_map.resize(m_split);

      // Process each original constraint and create mapping
      Index split_idx = 0;
      for (Index i = 0; i < m; ++i) {
         ConoptConstraintType type = constraint_type(i);
         split_constraint_map[i] = split_idx;  // Record first split index for this constraint

         if (type == ConoptConstraintType::RANGE) {
            // Both bounds: create two constraints
            original_constraint_map[split_idx] = i;
            split_idx++;
            original_constraint_map[split_idx] = i;
            split_idx++;
         }
         else {
            // for all other constraints we need one mapping.
            original_constraint_map[split_idx] = i;
            split_idx++;
         }
      }

      // Add the objective function as the last constraint (type FREE)
      objective_row_index = split_idx;
      original_constraint_map[split_idx] = -1;  // Special marker for objective (not a constraint)
      split_idx++;

      // Now split the Jacobian structure
      split_jacobian_structure();
   }

   /**
    * @brief Split the Jacobian structure to match split constraints
    */
   void split_jacobian_structure() {
      nnz_jac_g_split = 0;
      jacobian_split_map.clear();

      // Count non-zeros for split Jacobian
      for (Index k = 0; k < nnz_jac_g; ++k) {
         Index orig_row = jac_g_iRow[k];
         ConoptConstraintType type = constraint_type(orig_row);

         if (type == ConoptConstraintType::RANGE) {
            // Both bounds: this entry appears in both split constraints
            nnz_jac_g_split += 2;
         }
         else {
            // All other types: this entry appears once
            nnz_jac_g_split += 1;
         }
      }

      // Add dense objective function gradient (all n variables)
      nnz_jac_g_split += n;

      // Resize Jacobian mapping vectors
      jacobian_split_map.resize(nnz_jac_g_split);
      jacobian_split_rows.resize(nnz_jac_g_split);

      // Create split Jacobian mapping
      Index split_k = 0;
      for (Index k = 0; k < nnz_jac_g; ++k) {
         Index orig_row = jac_g_iRow[k];
         ConoptConstraintType type = constraint_type(orig_row);

         if (type == ConoptConstraintType::RANGE) {
            // Both bounds: create two entries
            Index first_split = split_constraint_map[orig_row];
            jacobian_split_map[split_k] = k;             // Map to original Jacobian entry
            jacobian_split_rows[split_k] = first_split;  // Lower bound constraint
            split_k++;
            jacobian_split_map[split_k] = k;                 // Map to original Jacobian entry
            jacobian_split_rows[split_k] = first_split + 1;  // Upper bound constraint
            split_k++;
         }
         else {
            // All other types: map to original entry
            Index split_row = split_constraint_map[orig_row];
            jacobian_split_map[split_k] = k;
            jacobian_split_rows[split_k] = split_row;
            split_k++;
         }
      }

      // Add objective function gradient entries (dense row)
      // These don't map to original Jacobian (they're new)
      for (Index j = 0; j < n; ++j) {
         jacobian_split_map[split_k] = -1;  // Special marker for objective gradient
         jacobian_split_rows[split_k] = objective_row_index;  // Objective constraint row
         split_k++;
      }
   }

   // === Helper methods for generating split data on-the-fly ===

   /**
    * @brief Get constraint type for a split constraint index
    */
   ConoptConstraintType get_split_constraint_type(Index split_idx) const {
      // if the index maps to the objective row, then we just return the
      // constraint type FREE
      if (split_idx == objective_row_index) {
         return ConoptConstraintType::FREE;
      }

      Index orig_idx = original_constraint_map[split_idx];
      ConoptConstraintType type = constraint_type(orig_idx);

      // if the constraint type is RANGE, then we need to get the split
      // inequalities. The first one is always >=, then second is <=.
      if (type == ConoptConstraintType::RANGE) {
         // Determine if this is the lower or upper bound constraint
         Index first_split = split_constraint_map[orig_idx];
         if (split_idx == first_split) {
            return ConoptConstraintType::GREATEREQ;
         }
         else {
            return ConoptConstraintType::LESSEQ;
         }
      }

      return type;
   }

   /**
    * @brief Get RHS value for a split constraint index
    */
   Number get_split_constraint_rhs(Index split_idx) const {
      // if the index maps to an objective row, then we just return 0.0
      if (split_idx == objective_row_index) {
         return 0.0;  // RHS doesn't matter for objective
      }

      Index orig_idx = original_constraint_map[split_idx];
      ConoptConstraintType type = constraint_type(orig_idx);

      // for RANGE constraints, the first is >=, then the second is <=
      if (type == ConoptConstraintType::RANGE) {
         // Determine if this is the lower or upper bound constraint
         Index first_split = split_constraint_map[orig_idx];
         if (split_idx == first_split) {
            return g_l[orig_idx];  // Lower bound (should be finite)
         }
         else {
            return g_u[orig_idx];  // Upper bound (should be finite)
         }
      }
      else if (type == ConoptConstraintType::EQUAL || type == ConoptConstraintType::GREATEREQ) {
         return g_l[orig_idx];  // Should be finite
      }
      else if (type == ConoptConstraintType::LESSEQ) {
         return g_u[orig_idx];  // Should be finite
      }
      else {
         return 0.0;  // FREE constraint
      }
   }

   /**
    * @brief Get split Jacobian row index for a split Jacobian entry
    */
   Index get_split_jacobian_row(Index split_k) const {
      return jacobian_split_rows[split_k];
   }

   /**
    * @brief Get split Jacobian column index for a split Jacobian entry
    */
   Index get_split_jacobian_col(Index split_k) const {
      if (jacobian_split_map[split_k] == -1) {
         // This is an objective gradient entry
         return split_k - (nnz_jac_g_split - n);  // Column index for objective
      }

      Index orig_k = jacobian_split_map[split_k];
      return jac_g_jCol[orig_k];
   }

   /**
    * @brief Check if all required information is available
    */
   bool is_complete() const {
      return (n > 0 && m >= 0 && x_l.size() == n && x_u.size() == n && g_l.size() == m &&
            g_u.size() == m && jac_g_iRow.size() == nnz_jac_g && jac_g_jCol.size() == nnz_jac_g);
   }

   /**
    * @brief Clear all data
    */
   void clear() {
      n = m = nnz_jac_g = nnz_h_lag = 0;
      m_split = nnz_jac_g_split = 0;
      objective_row_index = -1;
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
      hess_perm_sorted_to_orig.clear();

      // Clear split constraint mapping data
      original_constraint_map.clear();
      split_constraint_map.clear();
      jacobian_split_map.clear();
      jacobian_split_rows.clear();

      /* by default, as initial x is required. An initial z or lambda is not required. */
      init_x_req = true;
      init_z_req = init_lambda_req = false;
      has_variable_linearity = has_constraint_linearity = false;
      has_nonlinear_vars = false;
   }

   /**
    * @brief Get a string representation for debugging
    */
   std::string to_string() const {
      std::string result = "IpoptProblemInfo:\n";
      result += "  Dimensions: n=" + std::to_string(n) + ", m=" + std::to_string(m) +
            ", nnz_jac_g=" + std::to_string(nnz_jac_g) +
            ", nnz_h_lag=" + std::to_string(nnz_h_lag) + "\n";
      result += "  Split dimensions: m_split=" + std::to_string(m_split) +
            ", nnz_jac_g_split=" + std::to_string(nnz_jac_g_split) + "\n";
      result +=
            "  Index style: " + std::string(index_style == C_STYLE ? "C_STYLE" : "FORTRAN_STYLE") +
            "\n";
      result += "  Has initial x: " + std::string(init_x_req ? "yes" : "no") + "\n";
      result += "  Has initial z: " + std::string(init_z_req ? "yes" : "no") + "\n";
      result += "  Has initial lambda: " + std::string(init_lambda_req ? "yes" : "no") + "\n";
      result += "  Complete: " + std::string(is_complete() ? "yes" : "no") + "\n";
      return result;
   }
};

}  // namespace Ipopt

#endif  // IPOPT_PROBLEM_INFO_HPP
