/**
 * @file IpoptToConoptCallbacks.cpp
 * @brief Implements the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */

#include "IpoptToConoptCallbacks.hpp"
#include "IpTNLP.hpp"
#include "IpJournalist.hpp"
#include "Ipopt/IpSolveStatistics.hpp"
#include "Ipopt/IpOptionsList.hpp"
#include "Ipopt/IpIpoptData.hpp"
#include "IpAlgTypes.hpp"
#include "IpoptProblemInfo.hpp"
#include <cassert>
#include <vector>
#include <string>
#include <algorithm> /*  For std::max, std::copy, std::transform */
#include <cmath>     /*  For std::fabs */
#include <cstring>
#include <cstddef>
#include <memory>

/*  Helper functions for status code conversions */
namespace {

/**
 * @brief Convert CONOPT status codes to Ipopt ApplicationReturnStatus
 * @param modsta CONOPT model status
 * @param solsta CONOPT solver status
 * @return Corresponding Ipopt ApplicationReturnStatus
 */
Ipopt::ApplicationReturnStatus ConvertConoptToIpoptStatus(int modsta, int solsta) {
   /*  Handle explicit user stop by MODSTA regardless of SOLSTA */
   if (modsta == 10) {
      return Ipopt::User_Requested_Stop;
   }

   switch (solsta) {
   case 1: { /*  Normal completion */
      switch (modsta) {
      case 1: /*  Optimal */
      case 2: /*  Locally Optimal */
      case 7: /*  Feasible Solution */
         return Ipopt::Solve_Succeeded;
      case 4: /*  Locally Infeasible */
      case 5: /*  Infeasible */
         return Ipopt::Infeasible_Problem_Detected;
      default:
         return Ipopt::Solved_To_Acceptable_Level; /*  conservative success fallback */
      }
   }
   case 2: /*  Iteration Interrupt */
      return Ipopt::Maximum_Iterations_Exceeded;
   case 3: /*  Resource/CPU Time */
      return Ipopt::Maximum_CpuTime_Exceeded;
   case 4: { /*  Terminated by Solver */
      switch (modsta) {
      case 4: /*  Locally Infeasible (termination) */
         return Ipopt::Restoration_Failed;
      case 6: /*  Intermediate Infeasible (termination) */
         return Ipopt::Search_Direction_Becomes_Too_Small;
      case 3: /*  Unbounded (termination) */
         return Ipopt::Diverging_Iterates;
      case 5: /*  Infeasible (termination) */
         return Ipopt::Infeasible_Problem_Detected;
      case 13: /*  Error No Solution (termination) */
         return Ipopt::Error_In_Step_Computation;
      default:
         return Ipopt::Internal_Error;
      }
   }
   case 5: /*  Evaluation Error Limit */
      return Ipopt::Invalid_Number_Detected;
   case 6: /*  Error */
      return Ipopt::Internal_Error;
   case 8: /*  User Interrupt */
      return Ipopt::User_Requested_Stop;
   case 9: /*  Out of Memory */
      return Ipopt::Insufficient_Memory;
   case 10: { /*  System Error / Invalid setup/options */
      switch (modsta) {
      case 13:
         return Ipopt::Invalid_Problem_Definition; /*  or Invalid_Option */
      default:
         return Ipopt::Internal_Error;
      }
   }
   default:
      return Ipopt::Internal_Error;
   }
}

/**
 * @brief Convert Ipopt ApplicationReturnStatus to SolverReturn
 * @param status Ipopt ApplicationReturnStatus
 * @return Corresponding Ipopt SolverReturn
 */
Ipopt::SolverReturn ConvertApplicationStatusToSolverReturn(Ipopt::ApplicationReturnStatus status) {
   switch (status) {
   case Ipopt::Solve_Succeeded:
      return Ipopt::SUCCESS;
   case Ipopt::Solved_To_Acceptable_Level:
      return Ipopt::STOP_AT_ACCEPTABLE_POINT;
   case Ipopt::Infeasible_Problem_Detected:
      return Ipopt::LOCAL_INFEASIBILITY;
   case Ipopt::Search_Direction_Becomes_Too_Small:
      return Ipopt::STOP_AT_TINY_STEP;
   case Ipopt::Diverging_Iterates:
      return Ipopt::DIVERGING_ITERATES;
   case Ipopt::User_Requested_Stop:
      return Ipopt::USER_REQUESTED_STOP;
   case Ipopt::Feasible_Point_Found:
      return Ipopt::FEASIBLE_POINT_FOUND;
   case Ipopt::Maximum_Iterations_Exceeded:
      return Ipopt::MAXITER_EXCEEDED;
   case Ipopt::Maximum_CpuTime_Exceeded:
      return Ipopt::CPUTIME_EXCEEDED;
   case Ipopt::Error_In_Step_Computation:
      return Ipopt::ERROR_IN_STEP_COMPUTATION;
   case Ipopt::Invalid_Number_Detected:
      return Ipopt::INVALID_NUMBER_DETECTED;
   case Ipopt::Internal_Error:
      return Ipopt::INTERNAL_ERROR;
   default:
      return Ipopt::INTERNAL_ERROR;
   }
}

} /*  anonymous namespace */

/*  Helper to get the context struct from the cookie */
static inline IpoptConoptContext* GetContext(void* USRMEM) {
   assert(USRMEM != nullptr);
   return static_cast<IpoptConoptContext*>(USRMEM);
}

/*  Helper to get TNLP from the context */
static inline Ipopt::TNLP* GetTNLP(void* USRMEM) {
   return GetContext(USRMEM)->tnlp_;
}

/*  Helper to get Journalist */
static inline Ipopt::Journalist* GetJournalist(void* USRMEM) {
   return GetContext(USRMEM)->journalist_;
}

/*  Helper to get Problem Info */
static inline Ipopt::IpoptProblemInfo* GetProblemInfo(void* USRMEM) {
   return GetContext(USRMEM)->problem_info_;
}

/*  Cleanup function for context */
void CleanupIpoptConoptContext(IpoptConoptContext* context) {
   if (context && context->fdeval_cache_) {
      delete context->fdeval_cache_;
      context->fdeval_cache_ = nullptr;
   }
}

/*  Helper function to check if jacobian is cached */
bool IsJacobianCached(IpoptConoptContext* context) {
   if (!context || !context->fdeval_cache_) {
      return false;
   }

   FDEvalCache* cache = context->fdeval_cache_;
   return cache->isJacobianCached();
}

/**
 * @brief Populate IpoptData from CONOPT solution data
 * @param context The context containing CONOPT solution data and IpoptData instance
 * @param problem_info Problem information
 * @param status_sol Status and solution data from CONOPT
 */
static void PopulateIpoptDataFromConoptSolution(IpoptConoptContext* context,
      Ipopt::IpoptProblemInfo* problem_info, ConoptStatusSolution* status_sol) {
   if (!context || !context->ip_data_ || !problem_info || !status_sol) {
      return;
   }

   Ipopt::IpoptData* ip_data = context->ip_data_;

   // Set iteration count
   ip_data->Set_iter_count(status_sol->conopt_iter_);

   // Set tolerance (use default if not available from options)
   if (context->options_list_) {
      Ipopt::Number tol = 1e-8;
      context->options_list_->GetNumericValue("tol", tol, "");
      ip_data->Set_tol(tol);
   }

   // Set barrier parameter (mu) - CONOPT doesn't use barrier method, so set to 0
   ip_data->Set_mu(0.0);

   // Set fraction to boundary parameter (tau) - CONOPT doesn't use this, so set to 1.0
   ip_data->Set_tau(1.0);

   // Set info fields from CONOPT data
   ip_data->Set_info_alpha_primal(1.0);  // Assume full step for final solution
   ip_data->Set_info_alpha_dual(1.0);
   ip_data->Set_info_ls_count(0);

   // Note: We don't populate IteratesVector (curr_, trial_, etc.) because:
   // 1. CONOPT manages its own iterates
   // 2. The full IteratesVector structure requires complex initialization
   // 3. User callbacks typically only need iteration count, mu, tau, and tolerance
}

/**
 * @brief Populate IpoptData from CONOPT progress data
 * @param context The context containing CONOPT progress data and IpoptData instance
 * @param problem_info Problem information
 * @param iter Iteration number
 * @param objval Objective value
 * @param step Step length
 */
static void PopulateIpoptDataFromConoptProgress(IpoptConoptContext* context,
      Ipopt::IpoptProblemInfo* problem_info, int iter, double objval, double step) {
   if (!context || !context->ip_data_ || !problem_info) {
      return;
   }

   Ipopt::IpoptData* ip_data = context->ip_data_;

   // Set iteration count
   ip_data->Set_iter_count(static_cast<Ipopt::Index>(iter));

   // Set tolerance (use default if not available from options)
   if (context->options_list_) {
      Ipopt::Number tol = 1e-8;
      context->options_list_->GetNumericValue("tol", tol, "");
      ip_data->Set_tol(tol);
   }

   // Set barrier parameter (mu) - CONOPT doesn't use barrier method, so set to 0
   ip_data->Set_mu(0.0);

   // Set fraction to boundary parameter (tau) - CONOPT doesn't use this, so set to 1.0
   ip_data->Set_tau(1.0);

   // Set info fields from CONOPT progress data
   ip_data->Set_info_alpha_primal(step);  // Use step length as primal step size
   ip_data->Set_info_alpha_dual(step);    // Use step length as dual step size
   ip_data->Set_info_ls_count(0);
}

bool CallFinalizeSolutionWithCachedData(IpoptConoptContext* context) {
   bool success = false;

   if (!context || !context->tnlp_ || !context->status_solution_) {
      return false;
   }

   Ipopt::Journalist* jnlst = context->journalist_;
   Ipopt::IpoptProblemInfo* problem_info = context->problem_info_;
   ConoptStatusSolution* status_sol = context->status_solution_;

   if (!problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: ProblemInfo is NULL in CallFinalizeSolutionWithCachedData.\n");
      }
      return false;
   }

   /*  Check if we have cached status and solution data */
   if (!status_sol->status_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
               "CONOPT Bridge Warning: No status data cached for finalize_solution.\n");
      }
      return false;
   }

   if (!status_sol->solution_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
               "CONOPT Bridge Warning: No solution data cached for finalize_solution.\n");
      }
      return false;
   }

   /*  Translate CONOPT status codes to Ipopt ApplicationReturnStatus */
   Ipopt::ApplicationReturnStatus status =
         ConvertConoptToIpoptStatus(status_sol->conopt_modsta_, status_sol->conopt_solsta_);

   try {
      /*  Convert ApplicationReturnStatus to SolverReturn */
      Ipopt::SolverReturn solver_status = ConvertApplicationStatusToSolverReturn(status);

      /*  Populate IpoptData from CONOPT solution data */
      PopulateIpoptDataFromConoptSolution(context, problem_info, status_sol);

      /*  Call finalize_solution with the cached data and populated IpoptData */
      context->tnlp_->finalize_solution(
            solver_status, problem_info->n, status_sol->x_solution_.data(),
            status_sol->x_marginals_.data(), /*  z_L (lower bound multipliers) */
            status_sol->x_marginals_
                  .data(), /*  z_U (upper bound multipliers) - CONOPT provides combined marginals */
            problem_info->m, status_sol->y_solution_.data(), /*  g (constraint values) */
            status_sol->y_marginals_.data(),                 /*  lambda (constraint multipliers) */
            status_sol->conopt_objval_, context->ip_data_,   /*  ip_data - populated from CONOPT */
            nullptr /*  ip_cq - not available from CONOPT */
      );

      if (jnlst) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Bridge: finalize_solution called with status %d (MODSTA=%d, SOLSTA=%d).\n",
               static_cast<int>(status), status_sol->conopt_modsta_, status_sol->conopt_solsta_);
      }

      success = true;
   }
   catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Exception in finalize_solution: %s\n", e.what());
      }
      success = false;
   }
   catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Unknown exception in finalize_solution.\n");
      }
      success = false;
   }

   return success;
}

bool PopulateSolveStatistics(IpoptConoptContext* context) {
   bool success = false;

   if (!context || !context->stats_ || !context->status_solution_) {
      return false;
   }

   Ipopt::Journalist* jnlst = context->journalist_;
   ConoptStatusSolution* status_sol = context->status_solution_;

   /*  Check if we have cached data */
   if (!status_sol->status_cached_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
               "CONOPT Bridge Warning: No status data available for SolveStatistics.\n");
      }
      return false;
   }

   try {
      /*  Cast to our SolveStatistics implementation */
      Ipopt::SolveStatistics* stats = dynamic_cast<Ipopt::SolveStatistics*>(context->stats_);
      if (!stats) {
         if (jnlst) {
            jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Bridge Error: SolveStatistics is not a SolveStatistics instance.\n");
         }
         return false;
      }

      /*  Set basic solve information */
      stats->SetIterationCount(status_sol->conopt_iter_);
      stats->SetFinalObjective(status_sol->conopt_objval_);              /*  Original units */
      stats->SetFinalScaledObjective(status_sol->conopt_objval_scaled_); /*  As sent to CONOPT */

      /*  Translate CONOPT status to Ipopt solve status */
      Ipopt::ApplicationReturnStatus solve_status =
            ConvertConoptToIpoptStatus(status_sol->conopt_modsta_, status_sol->conopt_solsta_);
      stats->SetSolveStatus(solve_status);

      /*  Set convergence information */
      double dual_inf = 0.0;
      double constr_viol = 0.0;
      double bound_viol = 0.0;
      double complementarity = 0.0;
      double kkt_error = 0.0;

      if (status_sol->solution_cached_) {
         /*  Estimate primal infeasibility from constraint violations */
         double max_violation = 0.0;
         for (size_t i = 0; i < status_sol->y_solution_.size(); ++i) {
            /*  This is a simplified estimation - in practice you'd need to check against bounds */
            double violation = std::abs(status_sol->y_solution_[i]);
            max_violation = std::max(max_violation, violation);
         }
         constr_viol = max_violation;

         /*  Estimate dual infeasibility from marginal values */
         double max_dual_inf = 0.0;
         for (size_t i = 0; i < status_sol->x_marginals_.size(); ++i) {
            double dual_inf_val = std::abs(status_sol->x_marginals_[i]);
            max_dual_inf = std::max(max_dual_inf, dual_inf_val);
         }
         dual_inf = max_dual_inf;

         /*  Estimate bound violations (simplified) */
         bound_viol = 0.0; /*  Would need to check against variable bounds */

         /*  Estimate complementarity (simplified) */
         complementarity = 0.0; /*  Would need to compute from solution and bounds */

         /*  Estimate KKT error (simplified) */
         kkt_error = std::max(dual_inf, constr_viol);
      }

      stats->SetInfeasibilities(dual_inf, constr_viol, bound_viol, complementarity, kkt_error);
      stats->SetScaledInfeasibilities(
            dual_inf, constr_viol, bound_viol, complementarity, kkt_error); /*  Assume no scaling */

      /*
       * Function evaluation counts are now automatically tracked in FDEvalIni
       * Timing is now automatically tracked in OptimizeTNLP
       * No need to set these manually anymore
       */

      if (jnlst) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Bridge: SolveStatistics populated (iterations=%d, obj=%g, status=%d).\n",
               status_sol->conopt_iter_, status_sol->conopt_objval_,
               static_cast<int>(solve_status));
      }

      success = true;
   }
   catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge: Exception in PopulateSolveStatistics: %s\n", e.what());
      }
      success = false;
   }
   catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge: Unknown exception in PopulateSolveStatistics.\n");
      }
      success = false;
   }

   return success;
}

/*  --- Helper Functions for FDEval --- */

/**
 * @brief Get the objective scaling magnitude to apply to the objective value, gradient,
 * and Hessian contribution sent to CONOPT.
 *
 * Direction (the sign of obj_scaling) is handled separately via COIDEF_OptDir in
 * SetupConoptProblem; this returns only the magnitude, matching how TNLP::get_scaling_parameters()
 * documents obj_scaling: "if this number is chosen to be 10, then Ipopt solves internally an
 * optimization problem that has 10 times the value of the original objective function." Falls
 * back to 1.0 (no scaling) if obj_scaling is exactly zero, since scaling by zero would zero out
 * the objective entirely, which is never the intent.
 */
static double GetObjectiveScaleMagnitude(Ipopt::IpoptProblemInfo* problem_info) {
   if (!problem_info || problem_info->obj_scaling == 0.0) {
      return 1.0;
   }
   return std::fabs(problem_info->obj_scaling);
}

/**
 * @brief Evaluate a constraint's function value from FDEvalIni's cache
 *
 * The cache is keyed by ORIGINAL constraint index (see EvaluateAndCacheConstraints), so the
 * split-row-to-original translation happens here, on demand, for just this one row.
 *
 * @param context The context containing the cache
 * @param problem_info Problem information (for original_constraint_map)
 * @param jnlst Journalist for logging
 * @param conopt_constraint_idx CONOPT (split-row) constraint index
 * @param ROWNO Row number for error messages
 * @param G Output parameter for function value
 */
static void EvaluateConstraintValue(IpoptConoptContext* context,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst,
      Ipopt::Index conopt_constraint_idx, int ROWNO, double* G) {
   Ipopt::Index orig_idx = problem_info->original_constraint_map[conopt_constraint_idx];
   *G = GetCachedConstraintValue(context, orig_idx);

   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
            "CONOPT Bridge: Using cached constraint value for row %d.\n", ROWNO);
   }
}

/**
 * @brief Evaluate the objective function value directly from the TNLP.
 *
 * Unlike constraints - whose eval_g() always computes every row at once, so batching
 * them ahead of time in FDEvalIni is worthwhile - the objective is a single independent
 * eval_f() call. There is nothing to gain from caching it in FDEvalIni: it is simply
 * computed here, on demand, exactly when CONOPT's ROWNO says it is needed.
 *
 * @param context The context (for evaluation-count stats)
 * @param tnlp TNLP object for evaluation
 * @param problem_info Problem information
 * @param jnlst Journalist for logging
 * @param X Point of evaluation
 * @param ERRCNT Error counter (incremented on a recoverable evaluation failure)
 * @param G Output parameter for the objective value
 */
static void EvaluateObjectiveValue(IpoptConoptContext* context, Ipopt::TNLP* tnlp,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst, const double X[],
      int* ERRCNT, double* G) {
   Ipopt::Number obj_value;
   if (!tnlp->eval_f(problem_info->n, X, true, obj_value)) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Bridge: eval_f failed in FDEval.\n");
      }
      (*ERRCNT)++;
   }

   if (context->stats_) {
      context->stats_->IncrementObjectiveEvaluations();
   }

   *G = obj_value * GetObjectiveScaleMagnitude(problem_info);
}

/**
 * @brief Evaluate the objective gradient directly from the TNLP.
 *
 * @param context The context (for evaluation-count stats)
 * @param tnlp TNLP object for evaluation
 * @param problem_info Problem information
 * @param jnlst Journalist for logging
 * @param X Point of evaluation
 * @param ERRCNT Error counter (incremented on a recoverable evaluation failure)
 * @param JAC Output array for the objective gradient
 */
static void EvaluateObjectiveGradient(IpoptConoptContext* context, Ipopt::TNLP* tnlp,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst, const double X[],
      int* ERRCNT, double JAC[]) {
   std::vector<Ipopt::Number> gradient(problem_info->n);
   if (!tnlp->eval_grad_f(problem_info->n, X, true, gradient.data())) {
      if (jnlst) {
         jnlst->Printf(
               Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Bridge: eval_grad_f failed in FDEval.\n");
      }
      (*ERRCNT)++;
   }

   if (context->stats_) {
      context->stats_->IncrementObjectiveGradientEvaluations();
   }

   const double scale = GetObjectiveScaleMagnitude(problem_info);
   if (scale == 1.0) {
      std::copy(gradient.begin(), gradient.end(), JAC);
   }
   else {
      std::transform(gradient.begin(), gradient.end(), JAC, [scale](double g) { return g * scale; });
   }
}

/**
 * @brief Evaluate constraint Jacobian row from cache
 *
 * @param context The context containing the cache
 * @param problem_info Problem information
 * @param jnlst Journalist for logging
 * @param conopt_constraint_idx CONOPT constraint index
 * @param JAC Output array for Jacobian row
 * @return 0 on success, 1 on error
 */
static int EvaluateConstraintJacobianRow(IpoptConoptContext* context,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst,
      Ipopt::Index conopt_constraint_idx, double JAC[]) {
   int result = 0; /*  Default to success */

   /*
    * Handle constraint Jacobian row
    * Check if jacobian is cached first
    */
   if (IsJacobianCached(context)) {
      /*  Use the precomputed CSR-style row index to visit only the entries that belong to
       *  this row, instead of scanning the entire split Jacobian for matches.
       */
      const Ipopt::Index row_start = problem_info->jac_row_start[conopt_constraint_idx];
      const Ipopt::Index row_end = problem_info->jac_row_start[conopt_constraint_idx + 1];
      for (Ipopt::Index idx = row_start; idx < row_end && result == 0; ++idx) {
         const Ipopt::JacRowEntry& entry = problem_info->jac_row_entries[idx];
         Ipopt::Index split_col = entry.col;

         if (split_col >= 0 && split_col < problem_info->n) {
            /*  Trusts that this row's bucket never contains an objective-gradient marker
             *  entry: Conopt_FDEval only reaches this function in its !is_objective branch,
             *  and build_jac_row_index() groups CSR entries strictly by
             *  jacobian_split_rows[k], with every objective-marker entry's row forced to
             *  objective_row_index by split_jacobian_structure() - so a non-objective row's
             *  bucket can never contain one. Checked only via assert() (compiles away in
             *  release/NDEBUG builds), matching GetCachedJacobianValue's own trust model. */
            assert(entry.orig_k != -1 &&
                  "CONOPT Bridge: constraint mapping is not consistent (unexpected "
                  "objective entry found while evaluating a non-objective row)");
            JAC[split_col] = GetCachedJacobianValue(context, entry.orig_k);
         }
         else {
            /*  Should not happen if structure is correct - kept as a genuine runtime check
             *  (not an assert): split_col ultimately derives from the user's TNLP
             *  eval_jac_g structure callback, not just this bridge's own internal
             *  bookkeeping, and an out-of-range write into JAC here would be memory
             *  corruption, not just a bad read. */
            if (jnlst)
               jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                     "CONOPT Bridge Error: Invalid column index %d from split Jacobian "
                     "structure.\n",
                     split_col);
            result = 1; /*  there is an issue with the interface */
         }
      }
   }
   else {
      /*  This is an error - jacobian should be cached */
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: No cached jacobian. "
               "Jacobian was not cached in FDEvalIni.\n");
      }
      result = 1; /*  there is an issue with the interface */
   }

   return result;
}

/*  --- Helper Functions for FDEvalIni --- */

/**
 * @brief Evaluate and cache constraint function values
 *
 * eval_g() writes directly into the cache, keyed by ORIGINAL constraint index (not the
 * split row index CONOPT uses). Translating a split row to its original constraint - a
 * single array lookup via original_constraint_map - happens on demand in FDEval
 * (EvaluateConstraintValue) for whichever row is actually asked for, rather than eagerly
 * remapping every split row here regardless of whether FDEval will ever ask for it.
 *
 * @param context The context containing the cache
 * @param tnlp TNLP object for evaluation
 * @param problem_info Problem information
 * @param jnlst Journalist for logging
 * @param X Point of evaluation
 * @param ERRCNT Error counter (will be incremented on error)
 * @return true on success, false on error
 */
static bool EvaluateAndCacheConstraints(IpoptConoptContext* context, Ipopt::TNLP* tnlp,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst, const double X[],
      int* ERRCNT) {
   bool success = true;

   FDEvalCache* cache = context->fdeval_cache_;
   if (!tnlp->eval_g(
             problem_info->n, X, true, problem_info->m, cache->constraint_values_.data())) {
      if (jnlst) {
         jnlst->Printf(
               Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Bridge: eval_g failed in FDEvalIni.\n");
      }
      (*ERRCNT)++;
      success = false;
   }
   else {
      if (context->stats_) {
         context->stats_->IncrementConstraintEvaluations();
      }
   }

   return success;
}

/**
 * @brief Evaluate and cache constraint Jacobian
 * @param context The context containing the cache
 * @param tnlp TNLP object for evaluation
 * @param problem_info Problem information
 * @param jnlst Journalist for logging
 * @param X Point of evaluation
 * @param ERRCNT Error counter (will be incremented on error)
 * @return true on success, false on error
 */
static bool EvaluateAndCacheJacobian(IpoptConoptContext* context, Ipopt::TNLP* tnlp,
      Ipopt::IpoptProblemInfo* problem_info, Ipopt::Journalist* jnlst, const double X[],
      int* ERRCNT) {
   bool success = true;

   /*  Evaluate all jacobian values */
   /*  Note: According to Ipopt documentation, structure arrays (iRow, jCol) should only
    *  be provided on the first call. Subsequent calls for values should pass nullptr.
    *  The structure was already retrieved in RetrieveProblemInfo.
    */
   std::vector<Ipopt::Number> jacobian_values(problem_info->nnz_jac_g);

   if (!tnlp->eval_jac_g(problem_info->n, X, true, problem_info->m, problem_info->nnz_jac_g,
             nullptr, /*  Structure arrays not needed for value evaluation */
             nullptr, /*  Structure arrays not needed for value evaluation */
             jacobian_values.data())) {
      if (jnlst) {
         jnlst->Printf(
               Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Bridge: eval_jac_g failed in FDEvalIni.\n");
      }
      (*ERRCNT)++;
      success = false;
   }
   else {
      /*  Increment constraint jacobian evaluation count */
      if (context->stats_) {
         context->stats_->IncrementConstraintJacobianEvaluations();
      }
   }

   if (success) {
      /*  Cache the jacobian values */
      FDEvalCache* cache = context->fdeval_cache_;
      cache->cacheJacobian(jacobian_values);
   }

   return success;
}

/*  --- Implementation of the Trampolines --- */

/*
 * Note: Conopt_ReadMatrix is tricky because Ipopt doesn't have a direct equivalent.
 * Ipopt gets this info via get_nlp_info, get_bounds_info, and get_starting_point.
 * This trampoline might need to cache info from previous calls or might not be needed
 * if CONOPT can get the matrix structure differently. For now, let's assume it
 * needs to be implemented based on cached data or default values.
 */
int COI_CALLCONV Conopt_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[],
      int TYPE[], double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[],
      int NUMVAR, int NUMCON, int NUMNZ, void* USRMEM) {
   int result = 0; /*  Default to success */

   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Problem info is NULL in ReadMatrix.\n");
      }
      return 1; /*  Guard return */
   }

   /*  Verify dimensions match (use split dimensions for constraints) */
   if (NUMVAR != problem_info->n || NUMCON != problem_info->m_split ||
         NUMNZ != problem_info->nnz_jac_g_split) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Dimension mismatch in ReadMatrix. "
               "Expected: n=%d, m_split=%d, nnz_split=%d. Got: n=%d, m=%d, nnz=%d\n",
               problem_info->n, problem_info->m_split, problem_info->nnz_jac_g_split, NUMVAR,
               NUMCON, NUMNZ);
      }
      return 1; /*  Critical error */
   }

   try {
      /*  Populate variable bounds and current values */
      for (Ipopt::Index i = 0; i < NUMVAR; ++i) {
         LOWER[i] = problem_info->x_l[i];
         UPPER[i] = problem_info->x_u[i];
         CURR[i] = problem_info->x_init[i];
         if (CURR[i] > UPPER[i])
            CURR[i] = UPPER[i];
         else if (CURR[i] < LOWER[i])
            CURR[i] = LOWER[i];

         int varstatus = -1;
         if (CURR[i] == LOWER[i] || CURR[i] == UPPER[i])
            varstatus = 0;
         else
            varstatus = 1;
         VSTA[i] = varstatus;
      }

      /*  Populate constraint RHS values (generating split data on-the-fly) */
      for (Ipopt::Index i = 0; i < NUMCON; ++i) {
         TYPE[i] = int(problem_info->get_split_constraint_type(i));
         RHS[i] = problem_info->get_split_constraint_rhs(i);

         /*  the equation status is not being given at this stage, since this requires the functions
          *  to be evaluated first.
          *
          *  TODO: consider whether this is necessary.
          *  ESTA[i] = 0;
          */
      }

      /*
       * Populate Jacobian structure (using split Jacobian)
       * Convert from (iRow, jCol) pairs to CONOPT's sparse format
       * COLSTA[i] = starting index in ROWNO for column i
       * ROWNO[k] = row index for the k-th non-zero element
       */

      /*  Count non-zeros per column */
      std::vector<Ipopt::Index> col_counts(NUMVAR, 0);
      for (Ipopt::Index k = 0; k < NUMNZ; ++k) {
         Ipopt::Index col = problem_info->get_split_jacobian_col(k);
         if (col >= 0 && col < NUMVAR) {
            col_counts[col]++;
         }
      }

      /*  Set COLSTA as cumulative counts (starting indices) */
      Ipopt::Index current_pos = 0;
      for (Ipopt::Index j = 0; j < NUMVAR; ++j) {
         COLSTA[j] = current_pos;
         current_pos += col_counts[j];
      }
      /*  Set the last element to the total number of non-zeros */
      COLSTA[NUMVAR] = NUMNZ;

      /*  Fill ROWNO array with row indices, maintaining column order */
      std::vector<Ipopt::Index> col_positions(NUMVAR, 0);
      for (Ipopt::Index k = 0; k < NUMNZ; ++k) {
         Ipopt::Index row = problem_info->get_split_jacobian_row(k);
         Ipopt::Index col = problem_info->get_split_jacobian_col(k);

         if (col >= 0 && col < NUMVAR && row >= 0 && row < NUMCON) {
            Ipopt::Index pos = COLSTA[col] + col_positions[col];
            ROWNO[pos] = row;
            VALUE[pos] = 0.0; /*  Values will be computed by FDEval */
            /*  Set NLFLAG based on whether this entry is nonlinear */
            NLFLAG[pos] = problem_info->is_jacobian_entry_nonlinear(k) ? 1 : 0;
            col_positions[col]++;
         }
      }

      if (jnlst) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Bridge: ReadMatrix populated successfully with split constraints.\n");
      }
   }
   catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Exception in ReadMatrix: %s\n", e.what());
      }
      return 1; /*  Critical error */
   }

   return result;
}

int COI_CALLCONV Conopt_FDEval(const double X[], double* G, double JAC[], int ROWNO,
      const int JACNUM[], int MODE, int IGNERR, int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD,
      void* USRMEM) {
   int result = 0; /*  Default to success */

   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!tnlp || !problem_info) {
      /*  Log error even if journalist is null? Maybe stderr. */
      fprintf(stderr, "CONOPT Error: TNLP or ProblemInfo object is NULL in FDEval trampoline.\n");
      return 1; /*  Guard return */
   }

   /*  --- Adjust ROWNO based on Base --- */
   /*  CONOPT always uses 0-based indexing (C_STYLE) regardless of Ipopt's index style.
    *  The bridge converts FORTRAN indices (1-based) from Ipopt to C-style (0-based) for CONOPT.
    *  The objective function is now the last row in the constraint matrix.
    */
   bool is_objective = (ROWNO == problem_info->objective_row_index); /*  CONOPT uses 0-based */
   Ipopt::Index conopt_constraint_idx = -1;
   if (!is_objective) {
      conopt_constraint_idx = ROWNO; /*  CONOPT always uses 0-based indexing */
      if (conopt_constraint_idx < 0 || conopt_constraint_idx >= problem_info->m_split) {
         if (jnlst)
            jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Bridge Error: Invalid ROWNO %d received from CONOPT.\n", ROWNO);
         return 1; /*  Critical error */
      }
      /*  Map CONOPT constraint index to Ipopt constraint index using the mapping */
      /*  Note: original_constraint_map[conopt_constraint_idx] maps to the original Ipopt constraint
       */
   }

   try {
      /*  --- Evaluate Function Value (MODE 1 or 3) --- */
      if (MODE & 1) {
         if (is_objective) {
            EvaluateObjectiveValue(context, tnlp, problem_info, jnlst, X, ERRCNT, G);
         }
         else {
            EvaluateConstraintValue(context, problem_info, jnlst, conopt_constraint_idx, ROWNO, G);
         }
      }

      /*  --- Evaluate Derivatives (MODE 2 or 3) --- */
      if (MODE & 2) {
         if (is_objective) {
               EvaluateObjectiveGradient(context, tnlp, problem_info, jnlst, X, ERRCNT, JAC);
         }
         else {
            result = EvaluateConstraintJacobianRow(
                  context, problem_info, jnlst, conopt_constraint_idx, JAC);
         }
      }
   }
   catch (const std::exception& e) {
      /*  Catch exceptions from user code */
      if (jnlst)
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Exception caught in user NLP evaluation: %s\n", e.what());
      return 1;
   }
   catch (...) {
      if (jnlst)
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Unknown exception caught in user NLP evaluation.\n");
      return 1;
   }

   /*
    * CONOPT Docs: Return non-zero for serious/permanent error.
    * Evaluation errors are handled via ERRCNT.
    */
   return result;
}

/**
 * @brief CONOPT FDEvalIni callback - Function and Derivative Evaluator Initialization
 *
 * This optional callback is called each time the point of interest has changed, and it
 * defines the coming point and tells which rows CONOPT will need during the following
 * calls to FDEval. This implementation caches constraint values and the constraint
 * Jacobian - both of which require evaluating every row at once - if ROWLIST contains at
 * least one constraint row. The objective value/gradient are cheap single calls, so they
 * are simply evaluated on demand in FDEval instead of being cached here.
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
int COI_CALLCONV Conopt_FDEvalIni(const double X[], const int ROWLIST[], int MODE, int LISTSIZE,
      int NUMTHREAD, int IGNERR, int* ERRCNT, int NUMVAR, void* USRMEM) {
   /*  if there are no rows in the ROWLIST, then we don't evaluate any constraints or objective. */
   if (LISTSIZE == 0)
      return 0;

   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!tnlp || !problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: TNLP or ProblemInfo object is NULL in FDEvalIni.\n");
      }
      return 1; /*  Guard return */
   }
   /*  Process based on MODE: 1=function, 2=derivatives, 3=both */
   bool need_function_eval = (MODE & 1);
   bool need_derivative_eval = (MODE & 2);

   if (!need_function_eval && !need_derivative_eval) {
      /*  No evaluation needed */
      return 0; /*  Guard return */
   }

   /*  Get the cache (should already be initialized in OptimizeTNLP) */
   if (!context->fdeval_cache_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: FDEvalCache not initialized in FDEvalIni.\n");
      }
      return 1; /*  Critical error */
   }

   /*  The objective is evaluated live in FDEval (see EvaluateObjectiveValue), so FDEvalIni
    *  only needs to decide whether ROWLIST contains a genuine constraint row at all - eval_g
    *  and eval_jac_g always compute every original constraint at once, so there is no point
    *  calling them when ROWLIST only contains the objective row. */
   bool constraint_requested = false;
   for (int i = 0; i < LISTSIZE; ++i) {
      if (ROWLIST[i] != problem_info->objective_row_index) {
         constraint_requested = true;
         break;
      }
   }

   try {
      if (need_function_eval && constraint_requested) {
         EvaluateAndCacheConstraints(context, tnlp, problem_info, jnlst, X, ERRCNT);
      }

      if (need_derivative_eval && constraint_requested) {
         EvaluateAndCacheJacobian(context, tnlp, problem_info, jnlst, X, ERRCNT);
      }

      if (jnlst) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
               "CONOPT Bridge: FDEvalIni processed ROWLIST of size %d (constraints=%s).\n",
               LISTSIZE, constraint_requested ? "yes" : "no");
      }
   }
   catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Exception in FDEvalIni: %s\n", e.what());
      }
      (*ERRCNT)++;
   }
   catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Unknown exception in FDEvalIni.\n");
      }
      (*ERRCNT)++;
   }

   return 0;
}

/*  === Hessian of Lagrangian: Structure === */
int COI_CALLCONV Conopt_2DLagrStr(
      int HSRW[], int HSCL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM) {
   int result = 0; /*  Default to success */

   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = context->problem_info_;
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Problem info is NULL in 2DLagrStr.\n");
      }
      return 1; /*  Guard return */
   }

   /*  Basic consistency check */
   if (NHESS != problem_info->nnz_h_lag) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: NHESS mismatch in 2DLagrStr. Expected %d, got %d.\n",
               problem_info->nnz_h_lag, NHESS);
      }
      return 1;
   }

   /*  Fill row/col indices in sorted order using precomputed permutation */
   const Ipopt::Index nnz = problem_info->nnz_h_lag;
   const auto& perm = problem_info->hess_perm_sorted_to_orig;
   printf("Hessian Structure:\n");
   for (Ipopt::Index k = 0; k < nnz; ++k) {
      Ipopt::Index orig = perm[k];
      HSRW[k] = static_cast<int>(problem_info->hess_iRow[orig]);
      HSCL[k] = static_cast<int>(problem_info->hess_jCol[orig]);
      printf("%d %d\n", HSRW[k], HSCL[k]);
   }

   if (NODRV) {
      /*  0 indicates that we can compute all second derivatives */
      *NODRV = 0;
   }

   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
            "CONOPT Bridge: 2DLagrStr populated (NHESS=%d).\n", NHESS);
   }

   return result;
}

/*  === Hessian of Lagrangian: Values === */
int COI_CALLCONV Conopt_2DLagrVal(const double X[], const double U[], const int HSRW[],
      const int HSCL[], double HSVL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS,
      void* USRMEM) {
   int result = 0; /*  Default to success */

   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = context->tnlp_;
   Ipopt::IpoptProblemInfo* problem_info = context->problem_info_;
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!tnlp || !problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: TNLP or ProblemInfo is NULL in 2DLagrVal.\n");
      }
      return 1; /*  Guard return */
   }

   /*  Sanity checks */
   if (NUMVAR != problem_info->n || NUMCON != problem_info->m_split ||
         NHESS != problem_info->nnz_h_lag) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Dimension mismatch in 2DLagrVal (n=%d/%d, m_split=%d/%d, "
               "nnz_h=%d/%d).\n",
               NUMVAR, problem_info->n, NUMCON, problem_info->m_split, NHESS,
               problem_info->nnz_h_lag);
      }
      return 1;
   }

   /*  Build obj_factor from the multiplier of the objective pseudo-row */
   Ipopt::Number obj_factor = 0.0;
   if (problem_info->objective_row_index >= 0 && problem_info->objective_row_index < NUMCON) {
      obj_factor = U[problem_info->objective_row_index];
   }

   /*  Scale the objective's contribution to the Hessian by the same magnitude applied to
    *  the objective value/gradient in FDEval, so the second-order model CONOPT builds
    *  stays consistent with the (possibly scaled) first-order information it already has.
    */
   obj_factor *= GetObjectiveScaleMagnitude(problem_info);

   /*  Aggregate split-row multipliers into original constraint multipliers for TNLP */
   std::vector<Ipopt::Number> lambda(problem_info->m, 0.0);
   for (Ipopt::Index split_row = 0; split_row < problem_info->m_split; ++split_row) {
      if (split_row == problem_info->objective_row_index)
         continue; /*  skip objective row */
      int orig_idx = problem_info->original_constraint_map[split_row];
      if (orig_idx >= 0 && orig_idx < problem_info->m) {
         /*  Both parts of a RANGE add to the same original constraint Hessian contribution */
         lambda[orig_idx] += U[split_row];
      }
   }

   /*
    * Request Hessian values from TNLP at X with obj_factor and aggregated lambda
    * Values are expected in the same order as the structure we stored in problem_info_
    * Note: According to Ipopt documentation, structure arrays (iRow, jCol) should only
    * be provided on the first call. Subsequent calls for values should pass nullptr.
    * The structure was already retrieved in RetrieveProblemInfo.
    */
   std::vector<Ipopt::Number> values(problem_info->nnz_h_lag);

   bool ok = false;
   try {
      ok = tnlp->eval_h(problem_info->n, X, true, obj_factor, problem_info->m, lambda.data(), true,
            problem_info->nnz_h_lag,
            nullptr, /*  Structure arrays not needed for value evaluation */
            nullptr, /*  Structure arrays not needed for value evaluation */
            values.data());
      if (ok && GetContext(USRMEM)->stats_) {
         GetContext(USRMEM)->stats_->IncrementHessianEvaluations();
      }
   }
   catch (const std::exception& e) {
      if (jnlst)
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Exception in TNLP::eval_h: %s\n", e.what());
      return 1;
   }
   catch (...) {
      if (jnlst)
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Unknown exception in TNLP::eval_h.\n");
      return 1;
   }

   if (!ok) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP,
               "CONOPT Bridge: TNLP::eval_h returned false in 2DLagrVal.\n");
      }
      return 1; /*  treat as critical for CONOPT */
   }

   /*  Copy values into HSVL according to the sorted order recorded in 2DLagrStr */
   if (problem_info->hess_perm_sorted_to_orig.size() ==
         static_cast<size_t>(problem_info->nnz_h_lag)) {
      for (Ipopt::Index k = 0; k < problem_info->nnz_h_lag; ++k) {
         Ipopt::Index orig = problem_info->hess_perm_sorted_to_orig[k];
         HSVL[k] = values[orig];
      }
   }
   else {
      /*  Fallback: no permutation recorded, copy in original order */
      for (Ipopt::Index k = 0; k < problem_info->nnz_h_lag; ++k) {
         HSVL[k] = values[k];
      }
   }

   if (NODRV) {
      *NODRV = 0;
   }

   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP,
            "CONOPT Bridge: 2DLagrVal computed Hessian values (NHESS=%d).\n", NHESS);
   }

   return result;
}

int COI_CALLCONV Conopt_Status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM) {
   int result = 0; /*  Default to success */

   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!context || !context->status_solution_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Context or StatusSolution is NULL in Status.\n");
      }
      return 1; /*  Critical error */
   }

   /*  OBJVAL is reported by CONOPT in the same (possibly scaled) units that FDEval handed
    *  it - see GetObjectiveScaleMagnitude. Divide back out so conopt_objval_ matches what
    *  the TNLP's own eval_f would return, since that is what finalize_solution's obj_value
    *  parameter is documented to be. The raw, still-scaled value is kept separately for
    *  SolveStatistics::FinalScaledObjective. */
   const double scale = GetObjectiveScaleMagnitude(problem_info);

   /*  Cache the status information for later use in finalize_solution */
   context->status_solution_->conopt_modsta_ = MODSTA;
   context->status_solution_->conopt_solsta_ = SOLSTA;
   context->status_solution_->conopt_iter_ = ITER;
   context->status_solution_->conopt_objval_ = OBJVAL / scale;
   context->status_solution_->conopt_objval_scaled_ = OBJVAL;
   context->status_solution_->status_cached_ = true;

   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
            "CONOPT Bridge: Status cached (MODSTA=%d, SOLSTA=%d, ITER=%d, OBJVAL=%g).\n", MODSTA,
            SOLSTA, ITER, OBJVAL);
   }

   return result;
}

int COI_CALLCONV Conopt_Solution(const double XVAL[], const double XMAR[], const int XBAS[],
      const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[],
      const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM) {
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);

   if (!context || !problem_info || !context->status_solution_) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Context, ProblemInfo, or StatusSolution is NULL in Solution.\n");
      }
      return 1; /*  Critical error */
   }

   /*  Verify dimensions match */
   if (NUMVAR != problem_info->n || NUMCON != problem_info->m_split) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: Dimension mismatch in Solution. "
               "Expected: n=%d, m_split=%d. Got: n=%d, m=%d\n",
               problem_info->n, problem_info->m_split, NUMVAR, NUMCON);
      }
      return 1; /*  Critical error */
   }

   ConoptStatusSolution* status_sol = context->status_solution_;

   /*  Cache the solution data for later use in finalize_solution */
   status_sol->x_solution_.assign(XVAL, XVAL + NUMVAR);
   status_sol->x_marginals_.assign(XMAR, XMAR + NUMVAR);
   status_sol->x_basis_.assign(XBAS, XBAS + NUMVAR);
   status_sol->x_status_.assign(XSTA, XSTA + NUMVAR);

   /*  Map CONOPT constraint data to original constraint indices */
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
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
            "CONOPT Bridge: Solution cached (n=%d, m=%d).\n", NUMVAR, NUMCON);
   }

   return 0; /*  Success */
}

int COI_CALLCONV Conopt_Message(int SMSG, int DMSG, int NMSG, char* MSGV[], void* USRMEM) {
   int result = 0; /*  Default to success */

   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   if (!jnlst)
      return 0;

   /*
    * Process messages for each stream according to CONOPT documentation:
    * SMSG lines -> Screen file (immediate display, progress updates)
    * NMSG lines -> Status file (summary for solution status)
    * DMSG lines -> Documentation file (detailed iteration log and debugging)
    */

   /*  Process Screen file messages (SMSG lines) - immediate display */
   for (int i = 0; i < SMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      /*  Strip trailing newlines */
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
         msg.pop_back();
      }
      /*  Screen messages are for immediate display - use SUMMARY level */
      jnlst->Printf(Ipopt::J_ITERSUMMARY, Ipopt::J_MAIN, "%s\n", msg.c_str());
   }

   /*  Process Status file messages (NMSG lines) - solution summary */
   for (int i = 0; i < NMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
         msg.pop_back();
      }
      /*  Status messages are summaries - use SUMMARY level with SOLUTION category */
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_SOLUTION, "CONOPT Status: %s\n", msg.c_str());
   }

   /*  Process Documentation file messages (DMSG lines) - detailed debugging */
   for (int i = 0; i < DMSG && MSGV[i] != nullptr; ++i) {
      std::string msg = MSGV[i];
      while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
         msg.pop_back();
      }
      /*  Documentation messages are detailed - use DETAILED level with NLP category */
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_NLP, "CONOPT Debug: %s\n", msg.c_str());
   }

   /*  Flush buffer to ensure messages are displayed immediately */
   jnlst->FlushBuffer();

   return result;
}

int COI_CALLCONV Conopt_ErrMsg(int ROWNO, int COLNO, int POSNO, const char* MSG, void* USRMEM) {
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);
   if (!jnlst)
      return 0;

   /*  Clean up the message */
   std::string msg = MSG;
   while (!msg.empty() && (msg.back() == '\n' || msg.back() == '\r')) {
      msg.pop_back();
   }

   /*  Determine error category and level based on the type of error */
   Ipopt::EJournalLevel level = Ipopt::J_ERROR;
   Ipopt::EJournalCategory category =
         Ipopt::J_NLP; /*  Default to NLP category for model-specific errors */

   /*
    * CONOPT always uses C-style 0-based indexing, but if the user's TNLP uses FORTRAN indexing,
    * we need to convert the indices for display. ROWNO and COLNO from CONOPT are C-style (0-based).
    */
   int display_row = ROWNO;
   int display_col = COLNO;
   if (problem_info && problem_info->index_style == Ipopt::FORTRAN_STYLE) {
      /*  Convert C-style indices to FORTRAN-style (1-based) for display */
      if (ROWNO >= 0) {
         display_row = ROWNO + 1;
      }
      if (COLNO >= 0) {
         display_col = COLNO + 1;
      }
   }

   /*
    * Handle special cases for CONOPT (always uses C-style 0-based indexing):
    * COLNO = -1: message about a row (ROWNO between 0 and NumCon-1)
    * ROWNO = -1: message about a column (COLNO between 0 and NumVar-1)
    * ROWNO and COLNO both non-negative: message about Jacobian element or (row,column)-pair
    */

   if (COLNO == -1) {
      /*  Message about a row */
      if (ROWNO >= 0) {
         jnlst->Printf(
               level, category, "CONOPT Row Error (Row %d): %s\n", display_row, msg.c_str());
      }
      else {
         jnlst->Printf(level, category, "CONOPT Row Error: %s\n", msg.c_str());
      }
   }
   else if (ROWNO == -1) {
      /*  Message about a column */
      if (COLNO >= 0) {
         jnlst->Printf(
               level, category, "CONOPT Column Error (Col %d): %s\n", display_col, msg.c_str());
      }
      else {
         jnlst->Printf(level, category, "CONOPT Column Error: %s\n", msg.c_str());
      }
   }
   else if (ROWNO >= 0 && COLNO >= 0) {
      /*  Message about Jacobian element or (row,column)-pair */
      if (POSNO >= 0) {
         /*  Specific Jacobian element */
         jnlst->Printf(level, category, "CONOPT Jacobian Error (Row %d, Col %d, Pos %d): %s\n",
               display_row, display_col, POSNO, msg.c_str());
      }
      else {
         /*  (row,column)-pair */
         jnlst->Printf(level, category, "CONOPT Matrix Error (Row %d, Col %d): %s\n", display_row,
               display_col, msg.c_str());
      }
   }
   else {
      /*  General error message */
      jnlst->Printf(level, category, "CONOPT Error: %s\n", msg.c_str());
   }

   /*  Flush buffer to ensure error messages are displayed immediately */
   jnlst->FlushBuffer();

   return 0;
}

/**
 * @brief CONOPT Progress callback - Algorithmic Progress
 *
 * This optional callback is called at the end of each iteration to report
 * algorithmic progress. It bridges to Ipopt's intermediate_callback method,
 * allowing user code to monitor progress and request early termination.
 *
 * @param LEN_INT Number of elements in INT array (currently 5)
 * @param INT Integer array containing: [ITER, PHASE, NUMINF, NUMNOP, NSUPER]
 * @param LEN_RL Number of elements in RL array (currently 4)
 * @param RL Real array containing: [SUMINF, OBJVAL, RGMAX, STEP]
 * @param X Current point (read-only)
 * @param USRMEM User memory pointer containing IpoptConoptContext
 * @return 0 to continue, non-zero to request termination
 */
int COI_CALLCONV Conopt_Progress(
      int LEN_INT, const int INT[], int LEN_RL, const double RL[], const double X[], void* USRMEM) {
   IpoptConoptContext* context = GetContext(USRMEM);
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   Ipopt::Journalist* jnlst = GetJournalist(USRMEM);
   Ipopt::IpoptProblemInfo* problem_info = GetProblemInfo(USRMEM);

   if (!tnlp || !problem_info) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
               "CONOPT Bridge Error: TNLP or ProblemInfo is NULL in Progress.\n");
      }
      return 0; /*  Continue optimization despite error */
   }

   /*  Validate array lengths */
   if (LEN_INT < 5 || LEN_RL < 4) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_MAIN,
               "CONOPT Bridge Warning: Progress callback received insufficient data "
               "(LEN_INT=%d, LEN_RL=%d).\n",
               LEN_INT, LEN_RL);
      }
      return 0; /*  Continue optimization */
   }

   /*  Extract CONOPT progress information */
   const int iter = INT[0];  /*  ITER: Number of the iteration */
   const int phase = INT[1]; /*  PHASE: Phase of optimization (0-5) */
   /*  INT[2] = NUMINF: Number of infeasible constraints (not used) */
   /*  INT[3] = NUMNOP: Number of non-optimal variables (not used) */
   /*  INT[4] = NSUPER: Number of super-basic variables (not used) */

   const double suminf = RL[0]; /*  SUMINF: Sum of infeasibilities */
   const double objval = RL[1]; /*  OBJVAL: Value of true objective function */
   const double rgmax = RL[2];  /*  RGMAX: Numerically largest reduced gradient */
   const double step = RL[3];   /*  STEP: Optimal steplength */

   /*  Map CONOPT phase to Ipopt AlgorithmMode */
   /*  CONOPT phases: 0=infeasible Newton, 1-2=infeasible GRG, 3-4=feasible GRG, 5=SQP */
   /*  Ipopt modes: RegularMode=0, RestorationPhaseMode=1 */
   Ipopt::AlgorithmMode mode = Ipopt::RegularMode;
   if (phase == 0 || phase == 1 || phase == 2) {
      /*  Phases 0-2 are feasibility phases, similar to Ipopt's restoration phase */
      mode = Ipopt::RestorationPhaseMode;
   }

   /*  Map CONOPT data to Ipopt intermediate_callback parameters */
   /*  Many Ipopt parameters don't have direct CONOPT equivalents, so we use defaults */
   Ipopt::Index ipopt_iter = static_cast<Ipopt::Index>(iter);
   Ipopt::Number ipopt_obj_value = objval;
   Ipopt::Number ipopt_inf_pr = suminf; /*  Primal infeasibility */
   Ipopt::Number ipopt_inf_du = rgmax;  /*  Dual infeasibility (reduced gradient) */
   Ipopt::Number ipopt_mu = 0.0;        /*  Barrier parameter - not available from CONOPT */
   Ipopt::Number ipopt_d_norm = step;   /*  Step length approximates search direction norm */
   Ipopt::Number ipopt_regularization_size = 0.0; /*  Not available from CONOPT */
   Ipopt::Number ipopt_alpha_du = step;           /*  Use step length as approximation */
   Ipopt::Number ipopt_alpha_pr = step;           /*  Use step length as approximation */
   Ipopt::Index ipopt_ls_trials = 0; /*  Line search trials - not available from CONOPT */

   /*  Populate IpoptData from CONOPT progress data */
   PopulateIpoptDataFromConoptProgress(context, problem_info, iter, objval, step);

   bool should_continue = true;
   try {
      should_continue = tnlp->intermediate_callback(mode, ipopt_iter, ipopt_obj_value, ipopt_inf_pr,
            ipopt_inf_du, ipopt_mu, ipopt_d_norm, ipopt_regularization_size, ipopt_alpha_du,
            ipopt_alpha_pr, ipopt_ls_trials, context->ip_data_, nullptr);
   }
   catch (const std::exception& e) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Exception in intermediate_callback: %s\n", e.what());
      }
      /*  On exception, continue optimization */
      should_continue = true;
   }
   catch (...) {
      if (jnlst) {
         jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_USER_APPLICATION,
               "CONOPT Bridge: Unknown exception in intermediate_callback.\n");
      }
      /*  On exception, continue optimization */
      should_continue = true;
   }

   /*  CONOPT: return non-zero to request termination */
   /*  Ipopt: return false to request termination */
   /*  So we invert the boolean */
   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
            "CONOPT Bridge: Progress callback (iter=%d, phase=%d, obj=%g, suminf=%g, "
            "rgmax=%g, continue=%s)\n",
            iter, phase, objval, suminf, rgmax, should_continue ? "yes" : "no");
   }

   return should_continue ? 0 : 1; /*  0 = continue, 1 = stop */
}

/**
 * @brief CONOPT Option callback - Define Options
 *
 * This optional callback is called repeatedly by CONOPT to get non-default option values.
 * It returns one option per call until all options are provided (indicated by blank NAME).
 * The callback reads options from the OptionsList object stored in the context.
 *
 * @param NCALL Input: The sequence number of the call (0-based for C API)
 * @param RVAL Output: The value if the option is a real value
 * @param IVAL Output: The value if the option is an integer value
 * @param LVAL Output: The value if the option is a logical value (0=false, non-zero=true)
 * @param NAME Output: The CONOPT CR-cell name (8 characters, padded with blanks). Blank indicates
 * no more options.
 * @param USRMEM User memory pointer containing IpoptConoptContext
 * @return 0 on success, non-zero on error
 */
int COI_CALLCONV Conopt_Option(
      int NCALL, double* RVAL, int* IVAL, int* LVAL, char* NAME, void* USRMEM) {
   int result = 0; /* Default to success */

   /* Cast USRMEM from double* to void* (CONOPT uses double* for Fortran compatibility) */
   void* usrmem_void = static_cast<void*>(USRMEM);
   IpoptConoptContext* context = GetContext(usrmem_void);
   Ipopt::Journalist* jnlst = GetJournalist(usrmem_void);

   if (!context || !context->options_list_) {
      /* If no context or options list, return blank name (no more options) */
      if (NAME) {
         std::memset(NAME, ' ', 8);
         NAME[8] = '\0'; /* Null terminator to prevent buffer over-reads */
      }
      return 0;
   }

   /* Get the call number (0-based for C API) */
   int call_index = NCALL;

   /* Get the option for this call index */
   std::string opt_name;
   double rval = 0.0;
   int ival = 0;
   int lval = 0;

   /* GetConoptOption sets the correct value based on stored type
    * Since ConoptOption::Type is private, we pass nullptr for type
    */
   bool has_option = context->options_list_->GetConoptOption(
         call_index, opt_name, &rval, &ival, &lval, nullptr);

   if (!has_option) {
      /* No more options - return blank name */
      if (NAME) {
         std::memset(NAME, ' ', 8);
         NAME[8] = '\0'; /* Null terminator to prevent buffer over-reads */
      }
      if (jnlst) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Bridge: Option callback - no more options (NCALL=%d).\n", call_index);
      }
      return 0;
   }

   /* Copy the option name to NAME buffer (must be exactly 8 characters) */
   if (NAME) {
      /* Initialize entire buffer with blanks - CONOPT expects exactly 8 characters, Fortran-style
       */
      std::memset(NAME, ' ', 8);

      if (opt_name.length() > 0) {
         /* Ensure we have exactly 8 characters to copy */
         /* The opt_name should already be padded to 8 by pad_to_8(), but be defensive */
         size_t copy_len = std::min(opt_name.length(), static_cast<size_t>(8));

         /* Copy the characters - use data() to avoid null terminator issues */
         if (copy_len > 0) {
            std::memcpy(NAME, opt_name.data(), copy_len);
         }

         /* If the string was shorter than 8, the remaining bytes are already blanks from memset
          * above */
         /* If the string was exactly 8, we've copied all 8 bytes */
      }
      /* If opt_name was empty, NAME is already all blanks from memset */

      /* CRITICAL: Add null terminator to prevent CONOPT from reading beyond 8 bytes.
       * CONOPT's error reporting code appears to treat NAME as a C string, so we need
       * to null-terminate it even though the API docs say it's Fortran-style (8 bytes).
       * We assume CONOPT allocated at least 9 bytes for the buffer. */
      NAME[8] = '\0';
   }

   /* Set the appropriate value - GetConoptOption already set the correct one based on type
    * CONOPT will use the value matching the CR-cell type, so we set all three
    * (CONOPT ignores values that don't match the CR-cell type)
    */
   if (RVAL)
      *RVAL = rval;
   if (IVAL)
      *IVAL = ival;
   if (LVAL)
      *LVAL = lval;

   if (jnlst) {
      jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
            "CONOPT Bridge: Option callback - providing option %d: %s\n", call_index,
            opt_name.c_str());
   }

   return result;
}

/*
 * ... Implementations for other trampolines (SDDir, etc.) ...
 * These will follow the same pattern: cast USRMEM, call corresponding TNLP method
 * (if one exists), translate arguments/results. Many might be simple stubs
 * if Ipopt doesn't have a direct equivalent.
 */
