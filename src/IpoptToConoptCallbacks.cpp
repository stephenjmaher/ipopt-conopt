/**
 * @file IpoptToConoptCallbacks.cpp
 * @brief Implements the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */

#include "IpoptToConoptCallbacks.hpp"
#include "IpTNLP.hpp"
#include "IpJournalist.hpp"
#include "Ipopt/IpSolveStatistics.hpp"
#include "IpoptProblemInfo.hpp"
#include <cassert>
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

bool CallFinalizeSolutionWithCachedData(IpoptConoptContext* context) {
   if (!context || !context->tnlp_ || !context->status_solution_) {
      return false;
   }

   Ipopt::Journalist* jnlst = context->journalist_;
   Ipopt::IpoptProblemInfo* problem_info = context->problem_info_;
   ConoptStatusSolution* status_sol = context->status_solution_;

   if (!problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: ProblemInfo is NULL in CallFinalizeSolutionWithCachedData.\n");
      }
      return false;
   }

   // Check if we have cached status and solution data
   if (!status_sol->status_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
                      "CONOPT Shim Warning: No status data cached for finalize_solution.\n");
      }
      return false;
   }

   if (!status_sol->solution_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
                      "CONOPT Shim Warning: No solution data cached for finalize_solution.\n");
      }
      return false;
   }

   // Translate CONOPT status codes to Ipopt ApplicationReturnStatus
   Ipopt::ApplicationReturnStatus status;

   // Map CONOPT MODSTA (Model Status) to Ipopt ApplicationReturnStatus
   switch (status_sol->conopt_modsta_) {
      case 1:  // Optimal
      case 15: // Solved Unique
      case 16: // Solved
      case 17: // Solved Singular
         status = Ipopt::Solve_Succeeded;
         break;
      case 2:  // Locally optimal
         status = Ipopt::Solved_To_Acceptable_Level;
         break;
      case 3:  // Unbounded
         status = Ipopt::Diverging_Iterates;
         break;
      case 4:  // Infeasible
      case 5:  // Locally infeasible
      case 6:  // Intermediate infeasible
         status = Ipopt::Infeasible_Problem_Detected;
         break;
      case 7:  // Intermediate non-optimal
         status = Ipopt::Solved_To_Acceptable_Level;
         break;
      case 12: // Unknown type of error
         status = Ipopt::Internal_Error;
         break;
      case 13: // Error no solution
         status = Ipopt::Error_In_Step_Computation;
         break;
      default:
         // For unknown status codes, check SOLSTA (Solver Status)
         switch (status_sol->conopt_solsta_) {
            case 1:  // Normal completion
               status = Ipopt::Solve_Succeeded;
               break;
            case 2:  // Iteration interrupt
               status = Ipopt::Maximum_Iterations_Exceeded;
               break;
            case 3:  // Resource interrupt
               status = Ipopt::Maximum_CpuTime_Exceeded;
               break;
            case 4:  // Terminated by solver
               status = Ipopt::Search_Direction_Becomes_Too_Small;
               break;
            case 5:  // Evaluation error limit
               status = Ipopt::Invalid_Number_Detected;
               break;
            case 8:  // User Interrupt
               status = Ipopt::User_Requested_Stop;
               break;
            case 9:  // Error: Setup failure
            case 10: // Error: Solver failure
            case 11: // Error: Internal solver error
               status = Ipopt::Internal_Error;
               break;
            default:
               status = Ipopt::Internal_Error;
               break;
         }
         break;
   }

   try {
      // Convert ApplicationReturnStatus to SolverReturn
      Ipopt::SolverReturn solver_status;
      switch (status) {
         case Ipopt::Solve_Succeeded:
            solver_status = Ipopt::SUCCESS;
            break;
         case Ipopt::Solved_To_Acceptable_Level:
            solver_status = Ipopt::STOP_AT_ACCEPTABLE_POINT;
            break;
         case Ipopt::Infeasible_Problem_Detected:
            solver_status = Ipopt::LOCAL_INFEASIBILITY;
            break;
         case Ipopt::Search_Direction_Becomes_Too_Small:
            solver_status = Ipopt::STOP_AT_TINY_STEP;
            break;
         case Ipopt::Diverging_Iterates:
            solver_status = Ipopt::DIVERGING_ITERATES;
            break;
         case Ipopt::User_Requested_Stop:
            solver_status = Ipopt::USER_REQUESTED_STOP;
            break;
         case Ipopt::Feasible_Point_Found:
            solver_status = Ipopt::FEASIBLE_POINT_FOUND;
            break;
         case Ipopt::Maximum_Iterations_Exceeded:
            solver_status = Ipopt::MAXITER_EXCEEDED;
            break;
         case Ipopt::Maximum_CpuTime_Exceeded:
            solver_status = Ipopt::CPUTIME_EXCEEDED;
            break;
         case Ipopt::Error_In_Step_Computation:
            solver_status = Ipopt::ERROR_IN_STEP_COMPUTATION;
            break;
         case Ipopt::Invalid_Number_Detected:
            solver_status = Ipopt::INVALID_NUMBER_DETECTED;
            break;
         case Ipopt::Internal_Error:
            solver_status = Ipopt::INTERNAL_ERROR;
            break;
         default:
            solver_status = Ipopt::INTERNAL_ERROR;
            break;
      }

      // Call finalize_solution with the cached data
      context->tnlp_->finalize_solution(
         solver_status,
         problem_info->n,
         status_sol->x_solution_.data(),
         status_sol->x_marginals_.data(),  // z_L (lower bound multipliers)
         status_sol->x_marginals_.data(),  // z_U (upper bound multipliers) - CONOPT provides combined marginals
         problem_info->m,
         status_sol->y_solution_.data(),   // g (constraint values)
         status_sol->y_marginals_.data(),  // lambda (constraint multipliers)
         status_sol->conopt_objval_,
         nullptr, // ip_data - not available from CONOPT
         nullptr  // ip_cq - not available from CONOPT
      );

      if (jnlst) {
         jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                      "CONOPT Shim: finalize_solution called with status %d (MODSTA=%d, SOLSTA=%d).\n",
                      static_cast<int>(status), status_sol->conopt_modsta_, status_sol->conopt_solsta_);
      }

      return true;

   } catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                      "CONOPT Shim: Exception in finalize_solution: %s\n", e.what());
      }
      return false;
   } catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                      "CONOPT Shim: Unknown exception in finalize_solution.\n");
      }
      return false;
   }
}

bool PopulateSolveStatistics(IpoptConoptContext* context) {
   if (!context || !context->stats_ || !context->status_solution_) {
      return false;
   }

   Ipopt::Journalist* jnlst = context->journalist_;
   ConoptStatusSolution* status_sol = context->status_solution_;

   // Check if we have cached data
   if (!status_sol->status_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
                      "CONOPT Shim Warning: No status data available for SolveStatistics.\n");
      }
      return false;
   }

   try {
      // Cast to our SolveStatistics implementation
      Ipopt::SolveStatistics* stats = dynamic_cast<Ipopt::SolveStatistics*>(context->stats_);
      if (!stats) {
         if (jnlst) {
            jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                         "CONOPT Shim Error: SolveStatistics is not a SolveStatistics instance.\n");
         }
         return false;
      }

      // Set basic solve information
      stats->SetIterationCount(status_sol->conopt_iter_);
      stats->SetFinalObjective(status_sol->conopt_objval_);
      stats->SetFinalScaledObjective(status_sol->conopt_objval_); // Assume no scaling for now

      // Translate CONOPT status to Ipopt solve status
      Ipopt::ApplicationReturnStatus solve_status;
      switch (status_sol->conopt_modsta_) {
         case 1:  // Optimal
         case 15: // Solved Unique
         case 16: // Solved
         case 17: // Solved Singular
            solve_status = Ipopt::Solve_Succeeded;
            break;
         case 2:  // Locally optimal
            solve_status = Ipopt::Solved_To_Acceptable_Level;
            break;
         case 3:  // Unbounded
            solve_status = Ipopt::Diverging_Iterates;
            break;
         case 4:  // Infeasible
         case 5:  // Locally infeasible
         case 6:  // Intermediate infeasible
            solve_status = Ipopt::Infeasible_Problem_Detected;
            break;
         case 7:  // Intermediate non-optimal
            solve_status = Ipopt::Solved_To_Acceptable_Level;
            break;
         case 12: // Unknown type of error
            solve_status = Ipopt::Internal_Error;
            break;
         case 13: // Error no solution
            solve_status = Ipopt::Error_In_Step_Computation;
            break;
         default:
            // Check SOLSTA for additional status information
            switch (status_sol->conopt_solsta_) {
               case 1:  // Normal completion
                  solve_status = Ipopt::Solve_Succeeded;
                  break;
               case 2:  // Iteration interrupt
                  solve_status = Ipopt::Maximum_Iterations_Exceeded;
                  break;
               case 3:  // Resource interrupt
                  solve_status = Ipopt::Maximum_CpuTime_Exceeded;
                  break;
               case 4:  // Terminated by solver
                  solve_status = Ipopt::Search_Direction_Becomes_Too_Small;
                  break;
               case 5:  // Evaluation error limit
                  solve_status = Ipopt::Invalid_Number_Detected;
                  break;
               case 8:  // User Interrupt
                  solve_status = Ipopt::User_Requested_Stop;
                  break;
               case 9:  // Error: Setup failure
               case 10: // Error: Solver failure
               case 11: // Error: Internal solver error
                  solve_status = Ipopt::Internal_Error;
                  break;
               default:
                  solve_status = Ipopt::Internal_Error;
                  break;
            }
            break;
      }
      stats->SetSolveStatus(solve_status);

      // Set convergence information
      double dual_inf = 0.0;
      double constr_viol = 0.0;
      double bound_viol = 0.0;
      double complementarity = 0.0;
      double kkt_error = 0.0;

      if (status_sol->solution_cached_) {
         // Estimate primal infeasibility from constraint violations
         double max_violation = 0.0;
         for (size_t i = 0; i < status_sol->y_solution_.size(); ++i) {
            // This is a simplified estimation - in practice you'd need to check against bounds
            double violation = std::abs(status_sol->y_solution_[i]);
            max_violation = std::max(max_violation, violation);
         }
         constr_viol = max_violation;

         // Estimate dual infeasibility from marginal values
         double max_dual_inf = 0.0;
         for (size_t i = 0; i < status_sol->x_marginals_.size(); ++i) {
            double dual_inf_val = std::abs(status_sol->x_marginals_[i]);
            max_dual_inf = std::max(max_dual_inf, dual_inf_val);
         }
         dual_inf = max_dual_inf;

         // Estimate bound violations (simplified)
         bound_viol = 0.0; // Would need to check against variable bounds

         // Estimate complementarity (simplified)
          complementarity = 0.0; // Would need to compute from solution and bounds

         // Estimate KKT error (simplified)
         kkt_error = std::max(dual_inf, constr_viol);
      }

      stats->SetInfeasibilities(dual_inf, constr_viol, bound_viol, complementarity, kkt_error);
      stats->SetScaledInfeasibilities(dual_inf, constr_viol, bound_viol, complementarity, kkt_error); // Assume no scaling

      // Function evaluation counts are now automatically tracked in FDEvalIni
      // Timing is now automatically tracked in OptimizeTNLP
      // No need to set these manually anymore

      if (jnlst) {
         jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                      "CONOPT Shim: SolveStatistics populated (iterations=%d, obj=%g, status=%d).\n",
                      status_sol->conopt_iter_, status_sol->conopt_objval_, static_cast<int>(solve_status));
      }

      return true;

   } catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim: Exception in PopulateSolveStatistics: %s\n", e.what());
      }
      return false;
   } catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim: Unknown exception in PopulateSolveStatistics.\n");
      }
      return false;
   }
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
               return 1; // there is an error in the interface.
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
               return 1; // there is an error in the interface.
            }
         }
      }

      // --- Evaluate Derivatives (MODE 2 or 3) ---
      if (MODE & 2) {
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
                     return 1; // There is an issue with the interface
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
               return 1; // There is an issue with the interface
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
                              return 1; // There is an issue with the interface
                           }
                        } else {
                           if (jnlst) {
                              jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                           "CONOPT Shim Error: The constraint mapping is not consistent.\n");
                           }
                           return 1; // there is an issue with the interface.
                        }
                     } else {
                        // Should not happen if structure is correct
                        if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "CONOPT Shim Error: Invalid column index %d from split Jacobian structure.\n", split_col);
                        return 1; // there is an issue with the interface
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
               return 1; // there is an issue with the interface
            }
         }
      }
   } catch (const std::exception& e) {
      // Catch exceptions from user code
      if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                               "CONOPT Shim: Exception caught in user NLP evaluation: %s\n", e.what());
      return 1;
   } catch (...) {
      if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                               "CONOPT Shim: Unknown exception caught in user NLP evaluation.\n");
      return 1;
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
            } else {
               // Increment constraint evaluation count
               if (context->stats_) {
                  context->stats_->IncrementConstraintEvaluations();
               }
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
            } else {
               // Increment objective evaluation count
               if (context->stats_) {
                  context->stats_->IncrementObjectiveEvaluations();
               }
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
            } else {
               // Increment objective gradient evaluation count
               if (context->stats_) {
                  context->stats_->IncrementObjectiveGradientEvaluations();
               }
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
         } else {
            // Increment constraint jacobian evaluation count
            if (context->stats_) {
               context->stats_->IncrementConstraintJacobianEvaluations();
            }
         }

         // Cache the jacobian values
         cache->cacheJacobian(jacobian_values);

         if (jnlst) {
            jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_NLP,
                         "CONOPT Shim: FDEvalIni cached jacobian values%s.\n",
                         has_objective ? " and objective gradient" : "");
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
   } catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
                      "CONOPT Shim: Unknown exception in FDEvalIni.\n");
      }
      if (ERRCNT) (*ERRCNT)++;
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
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!context || !context->status_solution_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Context or StatusSolution is NULL in Status.\n");
      }
      return 1; // Critical error
   }

   // Cache the status information for later use in finalize_solution
   context->status_solution_->conopt_modsta_ = MODSTA;
   context->status_solution_->conopt_solsta_ = SOLSTA;
   context->status_solution_->conopt_iter_ = ITER;
   context->status_solution_->conopt_objval_ = OBJVAL;
   context->status_solution_->status_cached_ = true;

   if (jnlst) {
      jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                   "CONOPT Shim: Status cached (MODSTA=%d, SOLSTA=%d, ITER=%d, OBJVAL=%g).\n",
                   MODSTA, SOLSTA, ITER, OBJVAL);
   }

   return 0; // Success
}


int COI_CALLCONV Conopt_Solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[],
   const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void *USRMEM)
{
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!context || !problem_info || !context->status_solution_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Context, ProblemInfo, or StatusSolution is NULL in Solution.\n");
      }
      return 1; // Critical error
   }

   // Verify dimensions match
   if (NUMVAR != problem_info->n || NUMCON != problem_info->m_split) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                      "CONOPT Shim Error: Dimension mismatch in Solution. "
                      "Expected: n=%d, m_split=%d. Got: n=%d, m=%d\n",
                      problem_info->n, problem_info->m_split, NUMVAR, NUMCON);
      }
      return 1; // Critical error
   }

   ConoptStatusSolution* status_sol = context->status_solution_;

   // Cache the solution data for later use in finalize_solution
   status_sol->x_solution_.assign(XVAL, XVAL + NUMVAR);
   status_sol->x_marginals_.assign(XMAR, XMAR + NUMVAR);
   status_sol->x_basis_.assign(XBAS, XBAS + NUMVAR);
   status_sol->x_status_.assign(XSTA, XSTA + NUMVAR);

   // Map CONOPT constraint data to original constraint indices
   status_sol->y_solution_.resize(problem_info->m, 0.0);
   status_sol->y_marginals_.resize(problem_info->m, 0.0);
   status_sol->y_basis_.resize(problem_info->m, 0);
   status_sol->y_status_.resize(problem_info->m, 0);

   for (int i = 0; i < NUMCON; ++i) {
      int orig_idx = problem_info->original_constraint_map[i];
      if (orig_idx >= 0 && orig_idx < problem_info->m) {
         status_sol->y_solution_[orig_idx] = YVAL[i];
         status_sol->y_marginals_[orig_idx] = YMAR[i];
         status_sol->y_basis_[orig_idx] = YBAS[i];
         status_sol->y_status_[orig_idx] = YSTA[i];
      }
   }

   status_sol->solution_cached_ = true;

   if (jnlst) {
      jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                   "CONOPT Shim: Solution cached (n=%d, m=%d).\n", NUMVAR, NUMCON);
   }

   return 0; // Success
}

int COI_CALLCONV Conopt_Message(int SMSG, int DMSG, int NMSG, char *MSGV[], void *USRMEM)
{
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   if (!jnlst) return 0;

   // Process messages for each stream according to CONOPT documentation:
   // SMSG lines -> Screen file (immediate display, progress updates)
   // NMSG lines -> Status file (summary for solution status)
   // DMSG lines -> Documentation file (detailed iteration log and debugging)

   // Process Screen file messages (SMSG lines) - immediate display
   for (int i = 0; i < SMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      // Strip trailing newlines
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
          msg.pop_back();
      }
      // Screen messages are for immediate display - use SUMMARY level
      jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "CONOPT: %s\n", msg.c_str());
   }

   // Process Status file messages (NMSG lines) - solution summary
   for (int i = 0; i < NMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
          msg.pop_back();
      }
      // Status messages are summaries - use SUMMARY level with SOLUTION category
      jnlst->Printf(Ipopt::J_SUMMARY, Ipopt::J_SOLUTION, "CONOPT Status: %s\n", msg.c_str());
   }

   // Process Documentation file messages (DMSG lines) - detailed debugging
   for (int i = 0; i < DMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
          msg.pop_back();
      }
      // Documentation messages are detailed - use DETAILED level with NLP category
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP, "CONOPT Debug: %s\n", msg.c_str());
   }

   // Flush buffer to ensure messages are displayed immediately
   jnlst->FlushBuffer();

   return 0; // Success
}


int COI_CALLCONV Conopt_ErrMsg(int ROWNO, int COLNO, int POSNO, const char *MSG, void *USRMEM)
{
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   if (!jnlst) return 0;

   // Clean up the message
   std::string msg = MSG;
   while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
       msg.pop_back();
   }

   // Determine error category and level based on the type of error
   Ipopt::EJournalLevel level = Ipopt::J_ERROR;
   Ipopt::EJournalCategory category = Ipopt::J_NLP; // Default to NLP category for model-specific errors

   // Handle special cases for Base=0 (C conventions) as per CONOPT documentation:
   // COLNO = -1: message about a row (ROWNO between 0 and NumCon-1)
   // ROWNO = -1: message about a column (COLNO between 0 and NumVar-1)
   // ROWNO and COLNO both non-negative: message about Jacobian element or (row,column)-pair

   if (COLNO == -1) {
      // Message about a row
      if (ROWNO >= 0) {
         jnlst->Printf(level, category, "CONOPT Row Error (Row %d): %s\n", ROWNO, msg.c_str());
      } else {
         jnlst->Printf(level, category, "CONOPT Row Error: %s\n", msg.c_str());
      }
   } else if (ROWNO == -1) {
      // Message about a column
      if (COLNO >= 0) {
         jnlst->Printf(level, category, "CONOPT Column Error (Col %d): %s\n", COLNO, msg.c_str());
      } else {
         jnlst->Printf(level, category, "CONOPT Column Error: %s\n", msg.c_str());
      }
   } else if (ROWNO >= 0 && COLNO >= 0) {
      // Message about Jacobian element or (row,column)-pair
      if (POSNO >= 0) {
         // Specific Jacobian element
         jnlst->Printf(level, category, "CONOPT Jacobian Error (Row %d, Col %d, Pos %d): %s\n",
                       ROWNO, COLNO, POSNO, msg.c_str());
      } else {
         // (row,column)-pair
         jnlst->Printf(level, category, "CONOPT Matrix Error (Row %d, Col %d): %s\n",
                       ROWNO, COLNO, msg.c_str());
      }
   } else {
      // General error message
      jnlst->Printf(level, category, "CONOPT Error: %s\n", msg.c_str());
   }

   // Flush buffer to ensure error messages are displayed immediately
   jnlst->FlushBuffer();

   return 0; // Success
}


// ... Implementations for other trampolines (Progress, SDDir, etc.) ...
// These will follow the same pattern: cast USRMEM, call corresponding TNLP method
// (if one exists), translate arguments/results. Many might be simple stubs
// if Ipopt doesn't have a direct equivalent.
