#pragma once

/*
 * This is your drop-in replacement for Ipopt/IpoptApplication.hpp
 */

// 1. INCLUDE ORIGINAL IPOPT HEADERS (for utilities)
#include "IpSmartPtr.hpp"
#include "IpReferenced.hpp"
#include "IpReturnCodes.hpp"
#include "IpoptSolveStatistics.hpp" // <-- Include the original stats header
// ... etc ...

// 2. INCLUDE THE CONOPT C-API
#include "conopt.h"

// 3. INCLUDE OUR NEW TRAMPOLINE DECLARATIONS
#include "IpoptToConoptCallbacks.hpp"
#include "IpoptProblemInfo.hpp"

// 4. FORWARD DECLARE THE OTHER SHIM CLASSES
namespace Ipopt {
   class TNLP;
   class Journal;
   class Journalist;
   class RegisteredOptions;
   class SolveStatistics;
}

namespace Ipopt {

   // Forward declare our shim class
   class SolveStatistics;

   /**
    * @brief This is your shim class for IpoptApplication.
    * It holds the CONOPT C-API handle and wires up the callbacks.
    */
   class IpoptApplication : public ReferencedObject {

   private:
      /**
       * @brief The handle to the CONOPT C-API object.
       * (Using your _ suffix convention for member variables)
       */
      coiHandle_t cntvect_;
      SmartPtr<Journalist> jnlst_;
      IpoptConoptContext context_;

      // We'll also need to store the user's TNLP
      SmartPtr<TNLP> tnlp_;

      SmartPtr<SolveStatistics> stats_;

      /**
       * @brief Problem information retrieved from TNLP
       */
      IpoptProblemInfo problem_info_;

   public:
      /**
       * @brief Constructor: Create the CONOPT handle.
       */
      IpoptApplication() :
         cntvect_(nullptr),
         jnlst_(new Journalist())
      {
         COI_Create(&cntvect_);
         jnlst_->SetAllPrintLevels(J_SUMMARY);
         stats_ = new SolveStatistics(); // <-- Create the shim stats object
      }

      // Constructor accepting journalist
      IpoptApplication(SmartPtr<Journalist> jnlst) :
         cntvect_(nullptr),
         jnlst_(jnlst)
      {
           COI_Create(&cntvect_);
           if (!jnlst_) {
               jnlst_ = new Journalist();
               jnlst_->SetAllPrintLevels(J_SUMMARY);
           }
           stats_ = new SolveStatistics(); // <-- Create the shim stats object
      }

      /**
       * @brief Destructor: Free the CONOPT handle.
       */
      virtual ~IpoptApplication() {
         if (cntvect_) {
            COI_Free(&cntvect_);
         }
      }

      /**
       * @brief Initialize: Stubbed for now.
       * We would pass options to CONOPT here.
       */
      ApplicationReturnStatus Initialize() {
         // TODO: Process options, e.g. from a RegisteredOptions shim
         return Solve_Succeeded;
      }

      /**
       * @brief Retrieve problem information from the TNLP object
       * This method calls all the necessary TNLP methods to gather problem data
       */
      ApplicationReturnStatus RetrieveProblemInfo(SmartPtr<TNLP> tnlp) {
         if (!tnlp) {
            return Invalid_Problem_Definition;
         }

         // Clear any existing data
         problem_info_.clear();

         // 1. Get basic problem dimensions
         if (!tnlp->get_nlp_info(problem_info_.n, problem_info_.m,
                                problem_info_.nnz_jac_g, problem_info_.nnz_h_lag,
                                problem_info_.index_style)) {
            if (jnlst_) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Failed to get NLP info from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         // Resize vectors based on dimensions
         problem_info_.resize_vectors();

         // 2. Get bounds information
         if (!tnlp->get_bounds_info(problem_info_.n, problem_info_.x_l.data(),
                                   problem_info_.x_u.data(),
                                   problem_info_.m, problem_info_.g_l.data(),
                                   problem_info_.g_u.data())) {
            if (jnlst_) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Failed to get bounds info from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         // 3. Get starting point
         if (!tnlp->get_starting_point(problem_info_.n, problem_info_.has_initial_x,
                                      problem_info_.x_init.data(),
                                      problem_info_.has_initial_z,
                                      problem_info_.z_L_init.data(),
                                      problem_info_.z_U_init.data(),
                                      problem_info_.m, problem_info_.has_initial_lambda,
                                      problem_info_.lambda_init.data())) {
            if (jnlst_) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Failed to get starting point from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         // 4. Get Jacobian structure (if needed)
         if (problem_info_.nnz_jac_g > 0) {
            if (!tnlp->eval_jac_g(problem_info_.n, problem_info_.x_init.data(), true,
                                 problem_info_.m, problem_info_.nnz_jac_g,
                                 problem_info_.jac_g_iRow.data(),
                                 problem_info_.jac_g_jCol.data(), nullptr)) {
               if (jnlst_) {
                  jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                "CONOPT Shim: Failed to get Jacobian structure from TNLP.\n");
               }
               return Invalid_Problem_Definition;
            }
         }

         // 5. Get Hessian structure (if needed)
         if (problem_info_.nnz_h_lag > 0) {
            if (!tnlp->eval_h(problem_info_.n, problem_info_.x_init.data(), true, 1.0,
                             problem_info_.m, problem_info_.lambda_init.data(), true,
                             problem_info_.nnz_h_lag, problem_info_.hess_iRow.data(),
                             problem_info_.hess_jCol.data(), nullptr)) {
               if (jnlst_) {
                  jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                "CONOPT Shim: Failed to get Hessian structure from TNLP.\n");
               }
               return Invalid_Problem_Definition;
            }
         }

         // 6. Get optional information (these may not be implemented by all TNLPs)

         // Variable linearity
         if (tnlp->get_variables_linearity(problem_info_.n, problem_info_.var_linearity.data())) {
            problem_info_.has_variable_linearity = true;
         }

         // Constraint linearity
         if (tnlp->get_constraints_linearity(problem_info_.m, problem_info_.const_linearity.data())) {
            problem_info_.has_constraint_linearity = true;
         }

         // Number of nonlinear variables
         if (tnlp->get_number_of_nonlinear_variables(problem_info_.num_nonlin_vars)) {
            problem_info_.has_nonlinear_vars = true;
         }

         // 7. Split constraints for CONOPT (constraints with both bounds need to be duplicated)
         problem_info_.split_constraints();

         // Verify we have complete information
         if (!problem_info_.is_complete()) {
            if (jnlst_) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Incomplete problem information retrieved.\n");
            }
            return Invalid_Problem_Definition;
         }

         if (jnlst_) {
            jnlst_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
                          "CONOPT Shim: Problem info retrieved and constraints split successfully.\n%s",
                          problem_info_.to_string().c_str());
         }

         return Solve_Succeeded;
      }

      /**
       * @brief The main "Solve" method.
       * This is where we wire everything together.
       */
      ApplicationReturnStatus OptimizeTNLP(SmartPtr<TNLP> tnlp) {
         if (!cntvect_ || !tnlp || !jnlst_ || !stats_) {
            return Invalid_Problem_Definition;
         }

         // --- Retrieve problem information from TNLP ---
         ApplicationReturnStatus info_status = RetrieveProblemInfo(tnlp);
         if (info_status != Solve_Succeeded) {
            return info_status;
         }

         // --- Prepare the context cookie ---
         context_.tnlp_ = tnlp.get();
         context_.journalist_ = jnlst_.get();
         context_.stats_ = stats_.get();
         context_.problem_info_ = &problem_info_; // Add problem info to context

         COIDEF_UsrMem(cntvect_, &context_);

         // --- Register all our trampolines ---
         // These C functions are defined in IpoptToConoptCallbacks.cpp
         // and will call the virtual methods on the tnlp* cookie.

         // Mandatory Callbacks
         COIDEF_ReadMatrix(cntvect_, Conopt_ReadMatrix);
         COIDEF_FDEval(cntvect_, Conopt_FDEval);
         COIDEF_Status(cntvect_, Conopt_Status);
         COIDEF_Solution(cntvect_, Conopt_Solution);
         COIDEF_Message(cntvect_, Conopt_Message);
         COIDEF_ErrMsg(cntvect_, Conopt_ErrMsg);

         // Optional Callbacks
         COIDEF_FDEvalIni(cntvect_, Conopt_FDEvalIni);
         COIDEF_FDEvalEnd(cntvect_, Conopt_FDEvalEnd);

         // ... (Register all others: Progress, SDDir, etc.) ...

         // --- Pass problem dimensions to CONOPT (use split dimensions) ---
         COIDEF_NumVar(cntvect_, problem_info_.n);
         COIDEF_NumCon(cntvect_, problem_info_.m_split);
         COIDEF_NumNz(cntvect_, problem_info_.nnz_jac_g_split);
         COIDEF_NumNlNz(cntvect_, problem_info_.num_nonlin_vars);
         COIDEF_NumHess(cntvect_, problem_info_.nnz_h_lag);

         // --- Solve the problem ---
         COI_Solve(cntvect_);

         // ... (Get status from Conopt_Solution or Conopt_Status callbacks) ...
         // ... (Populate our stats_ object) ...

         // (Translate CONOPT status to Ipopt status)
         return Solve_Succeeded;
      }

      /**
       * @brief Getter for statistics (stubbed)
       */
      SmartPtr<SolveStatistics> Statistics() {
         return stats_;
      }

      // ... (All other public IpoptApplication methods) ...
   };

} // namespace Ipopt
