/**
 * @file IpoptToConoptCallbacks.cpp
 * @brief Implements the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */

#include "IpoptToConoptCallbacks.hpp"
#include "Ipopt/TNLP.hpp"
#include "IpJournalist.hpp"
#include <vector>
#include <string>
#include <algorithm> // For std::max

// Helper to get the context struct from the cookie
static inline IpoptConoptContext* GetContext(void *USRMEM) {
   assert(USRMEM != nullptr);
   return static_cast<IpoptConoptContext*>(USRMEM);
}

// Helper to get TNLP from the context
static inline Ipopt::TNLP* GetTNLP(void *USRMEM) {
   return GetContext(USRMEM)->tnlp_;
}

// Helper to get Journalist
static inline Ipopt::Journalist* GetJournalist(void *USRMEM) {
    return GetContext(USRMEM)->journalist_;
}

// Helper to get Problem Info
static inline Ipopt::IpoptProblemInfo* GetProblemInfo(void *USRMEM) {
    return GetContext(USRMEM)->problem_info_;
}

// Cleanup function for context
void CleanupIpoptConoptContext(IpoptConoptContext* context) {
    if (context && context->fdeval_cache_) {
        delete context->fdeval_cache_;
        context->fdeval_cache_ = nullptr;
    }
}

// Helper function to get cached constraint value
bool GetCachedConstraintValue(IpoptConoptContext* context, int row_idx, double& value) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    if (row_idx < 0 || row_idx >= cache->num_constraints_) {
        return false;
    }

    if (cache->constraint_valid_[row_idx]) {
        value = cache->constraint_values_[row_idx];
        return true;
    }

    return false;
}

// Helper function to get cached objective value
bool GetCachedObjectiveValue(IpoptConoptContext* context, double& value) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    if (cache->objective_valid_) {
        value = cache->objective_value_;
        return true;
    }

    return false;
}

// Helper function to get cached jacobian value
bool GetCachedJacobianValue(IpoptConoptContext* context, int jacobian_idx, double& value) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    return cache->getCachedJacobianValue(jacobian_idx, value);
}

// Helper function to check if jacobian is cached
bool IsJacobianCached(IpoptConoptContext* context) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    return cache->isJacobianCached();
}

// Helper function to get cached objective gradient value
bool GetCachedObjectiveGradientValue(IpoptConoptContext* context, int var_idx, double& value) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    return cache->getCachedObjectiveGradientValue(var_idx, value);
}

// Helper function to check if objective gradient is cached
bool IsObjectiveGradientCached(IpoptConoptContext* context) {
    if (!context || !context->fdeval_cache_) {
        return false;
    }

    FDEvalCache* cache = context->fdeval_cache_;
    return cache->isObjectiveGradientCached();
}

// --- Implementation of the Trampolines ---

// Note: Conopt_ReadMatrix is tricky because Ipopt doesn't have a direct equivalent.
// Ipopt gets this info via get_nlp_info, get_bounds_info, and get_starting_point.
// This trampoline might need to cache info from previous calls or might not be needed
// if CONOPT can get the matrix structure differently. For now, let's assume it
// needs to be implemented based on cached data or default values.
int COI_CALLCONV Conopt_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPE[],
   double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[], int NUMVAR, int NUMCON,
   int NUMNZ, void *USRMEM)
{
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Problem info is NULL in ReadMatrix.\n");
      }
      return 1; // Critical error
   }

   // Verify dimensions match (use split dimensions for constraints)
   if (NUMVAR != problem_info->n || NUMCON != problem_info->m_split || NUMNZ != problem_info->nnz_jac_g_split) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Dimension mismatch in ReadMatrix. "
                      "Expected: n=%d, m_split=%d, nnz_split=%d. Got: n=%d, m=%d, nnz=%d\n",
                      problem_info->n, problem_info->m_split, problem_info->nnz_jac_g_split,
                      NUMVAR, NUMCON, NUMNZ);
      }
      return 1; // Critical error
   }

   try {
      // Populate variable bounds and current values
      for (Ipopt::Index i = 0; i < NUMVAR; ++i) {
         LOWER[i] = problem_info->x_l[i];
         UPPER[i] = problem_info->x_u[i];
         CURR[i] = problem_info->x_init[i];

         int varstatus = -1;
         if (CURR[i] == LOWER[i] || CURR[i] == UPPER[i])
            varstatus = 0;
         else
            varstatus = 1;
         VSTA[i] = varstatus;
      }

      // Populate constraint RHS values (generating split data on-the-fly)
      for (Ipopt::Index i = 0; i < NUMCON; ++i) {
         TYPE[i] = int(problem_info->get_split_constraint_type(i));
         RHS[i] = problem_info->get_split_constraint_rhs(i);

         // the equation status is not being given at this stage, since this requires the functions to be evaluated
         // first.
         // TODO: consider whether this is necessary.
         //ESTA[i] = 0;
      }

      // Populate Jacobian structure (using split Jacobian)
      // Convert from (iRow, jCol) pairs to CONOPT's sparse format
      // COLSTA[i] = starting index in ROWNO for column i
      // ROWNO[k] = row index for the k-th non-zero element

      // Count non-zeros per column
      std::vector<Ipopt::Index> col_counts(NUMVAR, 0);
      for (Ipopt::Index k = 0; k < NUMNZ; ++k) {
         Ipopt::Index col = problem_info->get_split_jacobian_col(k);
         if (col >= 0 && col < NUMVAR) {
            col_counts[col]++;
         }
      }

      // Set COLSTA as cumulative counts (starting indices)
      Ipopt::Index current_pos = 0;
      for (Ipopt::Index j = 0; j < NUMVAR; ++j) {
         COLSTA[j] = current_pos;
         current_pos += col_counts[j];
      }

      // Fill ROWNO array with row indices, maintaining column order
      std::vector<Ipopt::Index> col_positions(NUMVAR, 0);
      for (Ipopt::Index k = 0; k < NUMNZ; ++k) {
         Ipopt::Index row = problem_info->get_split_jacobian_row(k);
         Ipopt::Index col = problem_info->get_split_jacobian_col(k);

         if (col >= 0 && col < NUMVAR && row >= 0 && row < NUMCON) {
            Ipopt::Index pos = COLSTA[col] + col_positions[col];
            ROWNO[pos] = row;
            VALUE[pos] = 0.0; // Values will be computed by FDEval
            NLFLAG[pos] = 1; // the NLFLAG is set to 1 for all values because Ipopt evaluates all expressions, even
                             // linear expressions.
            col_positions[col]++;
         }
      }

      if (jnlst) {
         jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                      "CONOPT Shim: ReadMatrix populated successfully with split constraints.\n");
      }

   } catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Exception in ReadMatrix: %s\n", e.what());
      }
      return 1; // Critical error
   }

   return 0; // Success
}

int COI_CALLCONV Conopt_FDEval(const double X[], double *G, double JAC[], int ROWNO, const int JACNUM[], int MODE,
   int IGNERR, int *ERRCNT, int NUMVAR, int NUMJAC, int THREAD, void *USRMEM)
{
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!tnlp || !problem_info) {
      // Log error even if journalist is null? Maybe stderr.
      fprintf(stderr, "CONOPT Error: TNLP or ProblemInfo object is NULL in FDEval trampoline.\n");
      if (ERRCNT) (*ERRCNT)++;
      return 1; // Indicate critical error
   }

   // --- Adjust ROWNO based on Base ---
   // CONOPT doc mentions Base determines if ROWNO is 0/1 based. Ipopt is always C-style (0-based).
   // We assume CONOPT uses 0-based indexing (C_STYLE).
   // The objective function is now the last row in the constraint matrix.
   bool is_objective = (ROWNO == problem_info->objective_row_index); // CONOPT is 0-based
   Ipopt::Index conopt_constraint_idx = -1;
   Ipopt::Index ipopt_constraint_idx = -1;
   if (!is_objective) {
      conopt_constraint_idx = ROWNO; // CONOPT uses 0-based indexing
      if (conopt_constraint_idx < 0 || conopt_constraint_idx >= problem_info->m_split) {
         if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "CONOPT Shim Error: Invalid ROWNO %d received from CONOPT.\n", ROWNO);
         if (ERRCNT) (*ERRCNT)++;
         return 1; // Critical error
      }
      // Map CONOPT constraint index to Ipopt constraint index using the mapping
      ipopt_constraint_idx = problem_info->original_constraint_map[conopt_constraint_idx];
   }


   // --- Assume new_x=true for simplicity ---
   // Properly handling NEWPT/new_x requires more state management.
   bool new_x = true;

   bool success = true;
   int evaluation_errors = 0;

   try {
      // --- Evaluate Function Value (MODE 1 or 3) ---
      if (MODE & 1) {
         if (is_objective) {
            // Check for cached objective value first
            double cached_obj_value;
            if (GetCachedObjectiveValue(context, cached_obj_value)) {
               *G = cached_obj_value;
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
                               "CONOPT Shim: Using cached objective value.\n");
               }
            } else {
               // This is an error - FDEval should only be called for objective if it was in ROWLIST
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                               "CONOPT Shim Error: No cached value for objective. "
                               "Objective was not in the ROWLIST from FDEvalIni.\n");
               }
               evaluation_errors++;
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Return success but with error count
            }
         } else { // It's a constraint
            // Check for cached constraint value first
            double cached_constraint_value;
            if (GetCachedConstraintValue(context, conopt_constraint_idx, cached_constraint_value)) {
               *G = cached_constraint_value;
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
                               "CONOPT Shim: Using cached constraint value for row %d.\n", ROWNO);
               }
            } else {
               // This is an error - FDEval should only be called for rows that were in ROWLIST
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                               "CONOPT Shim Error: No cached value for constraint row %d. "
                               "This row was not in the ROWLIST from FDEvalIni.\n", ROWNO);
               }
               evaluation_errors++;
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Return success but with error count
            }
         }
      }

      // --- Evaluate Derivatives (MODE 2 or 3) ---
      if (MODE & 2) {
          // Initialize Jacobian row/gradient to zero, as Ipopt only provides non-zeros
          for(Ipopt::Index j=0; j<problem_info->n; ++j) {
              JAC[j] = 0.0;
          }

         if (is_objective) {
            // Check for cached objective gradient first
            if (IsObjectiveGradientCached(context)) {
               // Use cached objective gradient
               for (Ipopt::Index j = 0; j < problem_info->n; ++j) {
                  double cached_value;
                  if (GetCachedObjectiveGradientValue(context, j, cached_value)) {
                     JAC[j] = cached_value;
                  } else {
                     // This is an error - gradient should be fully cached
                     if (jnlst) {
                        jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                     "CONOPT Shim Error: No cached value for objective gradient at variable %d. "
                                     "Objective gradient was not properly cached in FDEvalIni.\n", j);
                     }
                     evaluation_errors++;
                     if (ERRCNT) (*ERRCNT)++;
                     return 0; // Return success but with error count
                  }
               }
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
                               "CONOPT Shim: Using cached objective gradient.\n");
               }
            } else {
               // This is an error - objective gradient should be cached
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                               "CONOPT Shim Error: No cached objective gradient. "
                               "Objective gradient was not cached in FDEvalIni.\n");
               }
               evaluation_errors++;
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Return success but with error count
            }
         } else { // It's a constraint Jacobian row
               // Handle constraint Jacobian row
               // Check if jacobian is cached first
               if (IsJacobianCached(context)) {
                  // Use cached jacobian values
                  for (Ipopt::Index k = 0; k < problem_info->nnz_jac_g_split; ++k) {
                     // Get the split Jacobian row and column for this entry
                     Ipopt::Index split_row = problem_info->get_split_jacobian_row(k);
                     Ipopt::Index split_col = problem_info->get_split_jacobian_col(k);

                     if (split_row == conopt_constraint_idx) {
                        // This entry belongs to the requested constraint row
                        if (split_col >= 0 && split_col < problem_info->n) {
                           // Get the value from the cached jacobian if it's not an objective entry
                           if (problem_info->jacobian_split_map[k] != -1) {
                              Ipopt::Index orig_k = problem_info->jacobian_split_map[k];
                              double cached_value;
                              if (GetCachedJacobianValue(context, orig_k, cached_value)) {
                                 JAC[split_col] = cached_value;
                              } else {
                                 // This is an error - jacobian should be fully cached
                                 if (jnlst) {
                                    jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                                 "CONOPT Shim Error: No cached value for jacobian entry %d. "
                                                 "Jacobian was not properly cached in FDEvalIni.\n", orig_k);
                                 }
                                 evaluation_errors++;
                                 if (ERRCNT) (*ERRCNT)++;
                                 return 0; // Return success but with error count
                              }
                           } else {
                              // This is an objective gradient entry - should not happen for constraints
                              JAC[split_col] = 0.0;
                           }
                        } else {
                           // Should not happen if structure is correct
                           if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "CONOPT Shim Error: Invalid column index %d from split Jacobian structure.\n", split_col);
                           evaluation_errors++;
                        }
                     }
                  }
               } else {
                  // This is an error - jacobian should be cached
                  if (jnlst) {
                     jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                  "CONOPT Shim Error: No cached jacobian. "
                                  "Jacobian was not cached in FDEvalIni.\n");
                  }
                  evaluation_errors++;
                  if (ERRCNT) (*ERRCNT)++;
                  return 0; // Return success but with error count
               }
            }
         }
      }
   } catch (const std::exception& e) {
      // Catch exceptions from user code
      if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                               "CONOPT Shim: Exception caught in user NLP evaluation: %s\n", e.what());
      evaluation_errors++;
   } catch (...) {
      if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                               "CONOPT Shim: Unknown exception caught in user NLP evaluation.\n");
      evaluation_errors++;
   }

   // --- Final Error Handling ---
   if (evaluation_errors > 0) {
      success = false;
      if (!IGNERR) { // Only count if IGNERR is 0
         if (ERRCNT) (*ERRCNT) += evaluation_errors;
      }
   }

   // CONOPT Docs: Return non-zero for serious/permanent error.
   // Evaluation errors are handled via ERRCNT.
   return 0; // Indicate success (or recoverable evaluation error)
}


/**
 * @brief CONOPT FDEvalIni callback - Function and Derivative Evaluator Initialization
 *
 * This optional callback is called each time the point of interest has changed,
 * and it defines the coming point and tells which constraints CONOPT will need
 * during the following calls to FDEval. This implementation caches constraint
 * and objective values for the rows in ROWLIST to optimize subsequent FDEval calls.
 *
 * @param X Vector with the point of evaluation for future FDEval calls
 * @param ROWLIST List of constraints that will be evaluated in future FDEval calls
 * @param MODE Mode of evaluation (1=function, 2=derivatives, 3=both)
 * @param LISTSIZE Number of elements in ROWLIST
 * @param NUMTHREAD Number of threads that will be used for following FDEval calls
 * @param IGNERR Whether CONOPT assumes the point to be safe (0) or potentially unsafe (1)
 * @param ERRCNT Function evaluation error counter
 * @param NUMVAR Number of variables
 * @param USRMEM User memory pointer containing IpoptConoptContext
 * @return 0 on success, 1 on critical error
 */
int COI_CALLCONV Conopt_FDEvalIni(const double X[], const int ROWLIST[], int MODE, int LISTSIZE, int NUMTHREAD,
   int IGNERR, int *ERRCNT, int NUMVAR, void *USRMEM)
{
   assert(LISTSIZE > 0);
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!tnlp || !problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: TNLP or ProblemInfo object is NULL in FDEvalIni.\n");
      }
      return 1; // Critical error
   }

   // Process based on MODE: 1=function, 2=derivatives, 3=both
   bool need_function_eval = (MODE & 1);
   bool need_derivative_eval = (MODE & 2);

   if (!need_function_eval && !need_derivative_eval) {
      // No evaluation needed
      return 0;
   }

   // Get the cache (should already be initialized in OptimizeTNLP)
   if (!context->fdeval_cache_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: FDEvalCache not initialized in FDEvalIni.\n");
      }
      return 1; // Critical error
   }

   FDEvalCache* cache = context->fdeval_cache_;

   // Reset validity flags for all constraints
   cache->invalidateAll();

   // first checking if the objective is the only row needing evaluation. It is determined later if the objective
   // needs to be evaluated.
   bool has_constraints = true;
   bool has_objective = false;
   if (LISTSIZE == 1 && ROWLIST[0] == problem_info->objective_row_index)
   {
      has_constraints = false;
      has_objective = true;
   }

   try {
      // Only evaluate functions if MODE requires it
      if (need_function_eval) {
         if (has_constraints) {
            // Evaluate all original constraints
            std::vector<Ipopt::Number> all_g(problem_info->m);
            if (!tnlp->eval_g(problem_info->n, X, true, problem_info->m, all_g.data())) {
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP,
                               "CONOPT Shim: eval_g failed in FDEvalIni.\n");
               }
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Non-critical error, continue
            }

            // Store constraint values for rows in ROWLIST
            for (int i = 0; i < LISTSIZE; ++i) {
               int row_idx = ROWLIST[i];
               if (row_idx == problem_info->objective_row_index) {
                  // Skip objective for now, will handle separately. Marking that the objective needs evaluation.
                  has_objective = true;
                  continue;
               } else if (row_idx >= 0 && row_idx < problem_info->m_split) {
                  // Map CONOPT constraint index to Ipopt constraint index
                  Ipopt::Index ipopt_constraint_idx = problem_info->original_constraint_map[row_idx];
                  if (ipopt_constraint_idx >= 0 && ipopt_constraint_idx < problem_info->m) {
                     cache->constraint_values_[row_idx] = all_g[ipopt_constraint_idx];
                     cache->constraint_valid_[row_idx] = true;
                  }
               }
            }
         }

         // Evaluate objective if it's in the ROWLIST
         if (has_objective) {
            Ipopt::Number obj_value;
            if (!tnlp->eval_f(problem_info->n, X, true, obj_value)) {
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP,
                               "CONOPT Shim: eval_f failed in FDEvalIni.\n");
               }
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Non-critical error, continue
            }
            cache->objective_value_ = obj_value;
            cache->objective_valid_ = true;
         }
      }

      // Evaluate derivatives if needed (MODE 2 or 3)
      if (need_derivative_eval) {
         // Check if objective gradient is needed
         for (int i = 0; i < LISTSIZE; ++i) {
            if (ROWLIST[i] == problem_info->objective_row_index) {
               has_objective = true;
               break;
            }
         }

         // Evaluate objective gradient if needed
         if (has_objective) {
            std::vector<Ipopt::Number> objective_gradient(problem_info->n);
            if (!tnlp->eval_grad_f(problem_info->n, X, true, objective_gradient.data())) {
               if (jnlst) {
                  jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP,
                               "CONOPT Shim: eval_grad_f failed in FDEvalIni.\n");
               }
               if (ERRCNT) (*ERRCNT)++;
               return 0; // Non-critical error, continue
            }

            // Cache the objective gradient
            cache->cacheObjectiveGradient(objective_gradient);
         }

         // Evaluate all jacobian values
         std::vector<Ipopt::Number> jacobian_values(problem_info->nnz_jac_g);
         if (!tnlp->eval_jac_g(problem_info->n, X, true, problem_info->m, problem_info->nnz_jac_g,
                               problem_info->jac_g_iRow.data(), problem_info->jac_g_jCol.data(),
                               jacobian_values.data())) {
            if (jnlst) {
               jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP,
                            "CONOPT Shim: eval_jac_g failed in FDEvalIni.\n");
            }
            if (ERRCNT) (*ERRCNT)++;
            return 0; // Non-critical error, continue
         }

         // Cache the jacobian values
         cache->cacheJacobian(jacobian_values);

         if (jnlst) {
            jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_NLP,
                         "CONOPT Shim: FDEvalIni cached jacobian values%s.\n",
                         need_objective_gradient ? " and objective gradient" : "");
         }
      }

      if (jnlst) {
         if (need_function_eval) {
            // Count how many constraint values were actually cached
            int cached_constraints = 0;
            bool cached_objective = false;
            for (int i = 0; i < LISTSIZE; ++i) {
               int row_idx = ROWLIST[i];
               if (row_idx == problem_info->objective_row_index) {
                  cached_objective = true;
               } else if (row_idx >= 0 && row_idx < problem_info->m_split) {
                  cached_constraints++;
               }
            }
            jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_NLP,
                         "CONOPT Shim: FDEvalIni cached %d constraint values%s.\n",
                         cached_constraints, cached_objective ? " and objective" : "");
         }
         if (need_derivative_eval) {
            jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_NLP,
                         "CONOPT Shim: FDEvalIni cached jacobian values.\n");
         }
      }

   } catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                      "CONOPT Shim: Exception in FDEvalIni: %s\n", e.what());
      }
      if (ERRCNT) (*ERRCNT)++;
      return 0; // Non-critical error, continue
   } catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                      "CONOPT Shim: Unknown exception in FDEvalIni.\n");
      }
      if (ERRCNT) (*ERRCNT)++;
      return 0; // Non-critical error, continue
   }

   return 0; // Success
}


/**
 * @brief CONOPT FDEvalEnd callback - Function and Derivative Evaluator Termination
 *
 * This optional callback is called at the end of the function evaluation stage.
 * This implementation cleans up any cached data generated in FDEvalIni that was
 * used to improve the efficiency of function and derivative evaluation.
 *
 * @param IGNERR Whether CONOPT assumes the point to be safe (0) or potentially unsafe (1)
 * @param ERRCNT Function evaluation error counter
 * @param USRMEM User memory pointer containing IpoptConoptContext
 * @return 0 on success, 1 on critical error
 */
int COI_CALLCONV Conopt_FDEvalEnd(int IGNERR, int *ERRCNT, void *USRMEM)
{
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!context) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Context is NULL in FDEvalEnd.\n");
      }
      return 1; // Critical error
   }

   // Clean up the FDEval cache if it exists
   if (context->fdeval_cache_) {
      // Mark all constraint values as invalid
      context->fdeval_cache_->invalidateAll();

      if (jnlst) {
         jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_NLP,
                      "CONOPT Shim: FDEvalEnd invalidated cached values.\n");
      }
   }

   return 0; // Success
}

int COI_CALLCONV Conopt_Status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void *USRMEM)
{
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);

   // This maps directly to Ipopt's intermediate_callback.
   // We need to translate CONOPT's MODSTA/SOLSTA to Ipopt::AlgorithmMode.
   // We also need to fetch other values like inf_pr, inf_du etc.
   // CONOPT might provide these through other means or not at all.

   Ipopt::AlgorithmMode mode = Ipopt::RegularMode; // Placeholder
   Ipopt::Number inf_pr = 0.0; // Placeholder
   Ipopt::Number inf_du = 0.0; // Placeholder
   // ... other placeholders ...

   // NOTE: Ipopt's intermediate_callback expects non-const pointers for ip_data and ip_cq.
   // CONOPT doesn't seem to pass these equivalents. We might need to pass NULL
   // or create dummy objects if the user's callback expects them. Passing NULL for now.
   bool continue_solve = tnlp->intermediate_callback(
      mode, ITER, OBJVAL, inf_pr, inf_du,
      0.0 /*mu*/, 0.0 /*d_norm*/, 0.0 /*regu_size*/, 0.0 /*alpha_du*/, 0.0 /*alpha_pr*/,
      0 /*ls_trials*/, nullptr /*ip_data*/, nullptr /*ip_cq*/
   );

   return continue_solve ? 0 : 1; // Return 1 to signal termination request
}


int COI_CALLCONV Conopt_Solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[],
   const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void *USRMEM)
{
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);

   // This maps directly to Ipopt's finalize_solution.
   // We need to translate CONOPT's status (presumably stored elsewhere)
   // into Ipopt::SolverReturn status.

   Ipopt::SolverReturn status = Ipopt::Solve_Succeeded; // Placeholder
   Ipopt::Number obj_value = 0.0; // Placeholder - CONOPT might pass this separately

   // Note: Ipopt needs g(x) values, CONOPT provides YVAL (constraint levels?).
   //       Ipopt needs z_L/z_U, CONOPT provides XMAR (variable marginals?).
   //       Ipopt needs lambda, CONOPT provides YMAR (constraint marginals?).
   //       Careful mapping is needed.
   //       Also, Ipopt expects non-const pointers for ip_data and ip_cq. Passing NULL.

   tnlp->finalize_solution(
      status, NUMVAR, XVAL, XMAR /*assuming maps to z_L*/, XMAR /*assuming maps to z_U*/,
      NUMCON, YVAL /*assuming maps to g*/, YMAR /*assuming maps to lambda*/,
      obj_value, nullptr /*ip_data*/, nullptr /*ip_cq*/
   );

   return 0; // Indicate success
}

int COI_CALLCONV Conopt_Message(int SMSG, int DMSG, int NMSG, char *MSGV[], void *USRMEM)
{
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   if (!jnlst) return 0; // Or handle error

   // Determine the Ipopt print level - requires mapping CONOPT levels (SMSG, DMSG, NMSG?)
   // to Ipopt::EJournalLevel. Let's make a simple guess for now.
   Ipopt::EJournalLevel level = Ipopt::J_ITERSUMMARY; // Default level

   // Determine the category
   Ipopt::EJournalCategory category = Ipopt::J_MAIN; // Default category

   int max_msg = std::max({SMSG, DMSG, NMSG}); // Assuming these are counts? Or levels?
   for (int i = 0; i < max_msg && MSGV[i] != nullptr; ++i) {
      // Use Journalist's Printf. It handles checking the print level internally.
      // Need to strip trailing newlines if CONOPT adds them, as Printf might add its own.
      std::string msg = MSGV[i];
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
          msg.pop_back();
      }
      jnlst->Printf(level, category, "CONOPT: %s\n", msg.c_str());
   }
   // Flush the buffer if needed (might depend on CONOPT's message frequency)
   // jnlst->FlushBuffer();

   return 0; // Indicate success
}


int COI_CALLCONV Conopt_ErrMsg(int ROWNO, int COLNO, int POSNO, const char *MSG, void *USRMEM)
{
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   if (!jnlst) return 0;

   // Errors should likely go to a high print level
   Ipopt::EJournalLevel level = Ipopt::J_ERROR;
   Ipopt::EJournalCategory category = Ipopt::J_MAIN; // Or maybe J_NLP

   std::string msg = MSG;
   while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
       msg.pop_back();
   }

   jnlst->Printf(level, category, "CONOPT Error (Row %d, Col %d, Pos %d): %s\n",
                 ROWNO, COLNO, POSNO, msg.c_str());
   // jnlst->FlushBuffer();

   return 0; // Indicate success
}

// ... Implementations for other trampolines (Progress, SDDir, etc.) ...
// These will follow the same pattern: cast USRMEM, call corresponding TNLP method
// (if one exists), translate arguments/results. Many might be simple stubs
// if Ipopt doesn't have a direct equivalent.
