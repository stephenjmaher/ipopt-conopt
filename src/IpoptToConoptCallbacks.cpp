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

// --- Implementation of the Trampolines ---

// Note: Conopt_ReadMatrix is tricky because Ipopt doesn't have a direct equivalent.
// Ipopt gets this info via get_nlp_info, get_bounds_info, and get_starting_point.
// This trampoline might need to cache info from previous calls or might not be needed
// if CONOPT can get the matrix structure differently. For now, let's assume it
// needs to be implemented based on cached data or default values.
int COI_CALLCONV Conopt_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPEX[],
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

      // Populate constraint RHS values (using split constraints)
      for (Ipopt::Index i = 0; i < NUMCON; ++i) {
         RHS[i] = problem_info->g_rhs[i];
         ESTA[i] = 0; // Constraint status - CONOPT will set this
      }

      // Populate Jacobian structure (using split Jacobian)
      for (Ipopt::Index k = 0; k < NUMNZ; ++k) {
         ROWNO[k] = problem_info->jac_g_iRow_split[k];
         COLSTA[k] = problem_info->jac_g_jCol_split[k];
         VALUE[k] = 0.0; // Values will be computed by FDEval
         NLFLAG[k] = (problem_info->has_constraint_linearity &&
                      problem_info->const_linearity_split[problem_info->jac_g_iRow_split[k]] == Ipopt::NONLINEAR) ? 1 : 0;
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

   if (!tnlp) {
      // Log error even if journalist is null? Maybe stderr.
      fprintf(stderr, "CONOPT Error: TNLP object is NULL in FDEval trampoline.\n");
      if (ERRCNT) (*ERRCNT)++;
      return 1; // Indicate critical error
   }

   // --- Adjust ROWNO based on Base ---
   // CONOPT doc mentions Base determines if ROWNO is 0/1 based. Ipopt is always C-style (0-based).
   // Assuming CONOPT uses 0-based if C_STYLE was requested.
   // If CONOPT Obj Row is special (e.g., 0) and constraints start at 1, adjust here.
   // Let's assume ROWNO_in = 0 for objective, ROWNO_in = 1..M for constraints based on doc.
   // We will map this to Ipopt's 0..M-1 constraint indexing.
   // **CRITICAL**: Verify CONOPT's objective row index convention! Assuming 0 here.
   bool is_objective = (ROWNO_in == 0);
   Ipopt::Index constraint_idx = -1;
   if (!is_objective) {
      constraint_idx = ROWNO_in - 1; // Assuming constraints are 1-based from CONOPT
      if (constraint_idx < 0 || constraint_idx >= g_m) {
         if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "CONOPT Shim Error: Invalid ROWNO %d received from CONOPT.\n", ROWNO_in);
         if (ERRCNT) (*ERRCNT)++;
         return 1; // Critical error
      }
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
            Ipopt::Number obj_value;
            if (!tnlp->eval_f(g_n, X, new_x, obj_value)) {
               evaluation_errors++;
               if (jnlst) jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Shim: User's eval_f failed.\n");
            } else {
               *G = obj_value;
            }
         } else { // It's a constraint
            // Need temporary storage for *all* constraints
            std::vector<Ipopt::Number> all_g(g_m);
            if (!tnlp->eval_g(g_n, X, new_x, g_m, all_g.data())) {
               evaluation_errors++;
               if (jnlst) jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Shim: User's eval_g failed for row %d.\n", ROWNO_in);
            } else {
               *G = all_g[constraint_idx]; // Store the requested constraint value
            }
         }
      }

      // --- Evaluate Derivatives (MODE 2 or 3) ---
      if (MODE & 2) {
          // Initialize Jacobian row/gradient to zero, as Ipopt only provides non-zeros
          for(Ipopt::Index j=0; j<g_n; ++j) {
              JAC[j] = 0.0;
          }

         if (is_objective) {
            // Get the full gradient into the dense JAC vector provided by CONOPT
            if (!tnlp->eval_grad_f(g_n, X, new_x, JAC)) {
               evaluation_errors++;
               if (jnlst) jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Shim: User's eval_grad_f failed.\n");
            }
         } else { // It's a constraint Jacobian row
            // 1. Get *all* Jacobian values from Ipopt into our cached buffer
            if (!tnlp->eval_jac_g(g_n, X, new_x, g_m, g_nnz_jac_g,
                                 nullptr, nullptr, g_jac_g_values.data())) {
               evaluation_errors++;
               if (jnlst) jnlst->Printf(Ipopt::J_WARNING, Ipopt::J_NLP, "CONOPT Shim: User's eval_jac_g failed fetching values.\n");
            } else {
               // 2. Populate the dense CONOPT JAC vector for the requested row
               for (Ipopt::Index k = 0; k < g_nnz_jac_g; ++k) {
                  // Check if this non-zero belongs to the requested constraint row
                  // Adjusting for potential 1-based indexing from Ipopt structure if needed
                  Ipopt::Index row_idx = g_jac_g_iRow[k] - (g_index_style == Ipopt::TNLP::FORTRAN_STYLE ? 1 : 0);
                  if (row_idx == constraint_idx) {
                     Ipopt::Index col_idx = g_jac_g_jCol[k] - (g_index_style == Ipopt::TNLP::FORTRAN_STYLE ? 1 : 0);
                     if (col_idx >= 0 && col_idx < g_n) {
                         JAC[col_idx] = g_jac_g_values[k];
                     } else {
                        // Should not happen if structure is correct
                         if (jnlst) jnlst->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "CONOPT Shim Error: Invalid column index %d from Ipopt structure.\n", col_idx);
                         evaluation_errors++;
                     }
                  }
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


int COI_CALLCONV Conopt_FDEvalIni(const double X[], const int ROWLIST[], int MODE, int LISTSIZE, int NUMTHREAD,
   int IGNERR, int *ERRCNT, int NUMVAR, void *USRMEM)
{
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   // Ipopt's TNLP doesn't have direct equivalents for FDEvalIni/End.
   // These might be for CONOPT's internal setup/teardown for parallel evaluations.
   // We might not need to call anything on the TNLP object here.
   // Consult CONOPT docs for what this callback is expected to do.
   return 0; // Indicate success
}


int COI_CALLCONV Conopt_FDEvalEnd(int IGNERR, int *ERRCNT, void *USRMEM)
{
   Ipopt::TNLP* tnlp = GetTNLP(USRMEM);
   // See comment in Conopt_FDEvalIni.
   return 0; // Indicate success
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
