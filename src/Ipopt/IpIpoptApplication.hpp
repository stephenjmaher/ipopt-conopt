#pragma once

/*
 * This is your drop-in replacement for Ipopt/IpoptApplication.hpp
 */

// 1. INCLUDE ORIGINAL IPOPT HEADERS (for utilities)
#include "IpSmartPtr.hpp"
#include "IpReferenced.hpp"
#include "IpReturnCodes.hpp"
#include "IpSolveStatistics.hpp"
#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
// ... etc ...

// 2. INCLUDE THE CONOPT C-API
#include "conopt.h"

// 3. INCLUDE OUR NEW TRAMPOLINE DECLARATIONS
#include "IpoptToConoptCallbacks.hpp"
#include "IpoptProblemInfo.hpp"
#include <cassert>
#include <vector>
#include <string>

// 4. FORWARD DECLARE THE OTHER SHIM CLASSES
namespace Ipopt {
   class TNLP;
   class Journal;
   class Journalist;
   class RegisteredOptions;
   class SolveStatistics;

   // Dummy OptionsList class
   class OptionsList : public ReferencedObject {
   public:
      void SetNumericValue(const std::string& name, Number value) {
         // Stubbed - do nothing
      }

      void SetStringValue(const std::string& name, const std::string& value) {
         // Stubbed - do nothing
      }
   };

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
         stats_ = new SolveStatistics(); // <-- Create the SolveStatistics object

         // Set verbose output for debugging
         if (!IsNull(jnlst_)) {
            // Add a console journal with detailed output
            jnlst_->AddFileJournal("console", "stdout", Ipopt::J_DETAILED);
         }
      }

      // Constructor accepting journalist
      IpoptApplication(SmartPtr<Journalist> jnlst) :
         cntvect_(nullptr),
         jnlst_(jnlst)
      {
           COI_Create(&cntvect_);
           if (IsNull(jnlst_)) {
               jnlst_ = new Journalist();
           }
           stats_ = new SolveStatistics(); // <-- Create the SolveStatistics object

           // Set verbose output for debugging
           if (!IsNull(jnlst_)) {
               // Add a console journal with detailed output
               jnlst_->AddFileJournal("console", "stdout", Ipopt::J_DETAILED);
           }
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
         if (IsNull(tnlp)) {
            return Invalid_Problem_Definition;
         }

         // Clear any existing data
         problem_info_.clear();

         // 1. Get basic problem dimensions
         TNLP::IndexStyleEnum index_style;
         if (!tnlp->get_nlp_info(problem_info_.n, problem_info_.m,
                                problem_info_.nnz_jac_g, problem_info_.nnz_h_lag,
                                index_style)) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Failed to get NLP info from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         // Store the index style and check for FORTRAN style (not supported)
         if (index_style == TNLP::FORTRAN_STYLE) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim Error: FORTRAN-style indexing is not supported for C/C++ projects. "
                             "Please use C-style (0-based) indexing in your TNLP implementation.\n");
            }
            return Invalid_Problem_Definition;
         }
         problem_info_.index_style = C_STYLE;

         // Resize vectors based on dimensions
         problem_info_.resize_vectors();

         // 2. Get bounds information
         if (!tnlp->get_bounds_info(problem_info_.n, problem_info_.x_l.data(),
                                   problem_info_.x_u.data(),
                                   problem_info_.m, problem_info_.g_l.data(),
                                   problem_info_.g_u.data())) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Failed to get bounds info from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         // 3. Get starting point
         if (!tnlp->get_starting_point(problem_info_.n, problem_info_.init_x_req,
                                      problem_info_.x_init.data(),
                                      problem_info_.init_z_req,
                                      problem_info_.z_L_init.data(),
                                      problem_info_.z_U_init.data(),
                                      problem_info_.m, problem_info_.init_lambda_req,
                                      problem_info_.lambda_init.data())) {
            if (!IsNull(jnlst_)) {
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
               if (!IsNull(jnlst_)) {
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
               if (!IsNull(jnlst_)) {
                  jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                                "CONOPT Shim: Failed to get Hessian structure from TNLP.\n");
               }
               return Invalid_Problem_Definition;
            }
         }

         // 6. Get optional information (these may not be implemented by all TNLPs)

         // Variable linearity - populate with dummy data since we don't use it
         if (problem_info_.n > 0) {
            std::vector<TNLP::LinearityType> var_linearity(problem_info_.n);
            tnlp->get_variables_linearity(problem_info_.n, var_linearity.data());
            // Fill with dummy data
            for (Index i = 0; i < problem_info_.n; ++i) {
               problem_info_.var_linearity[i] = NONLINEAR; // Default to non-linear
            }
         }

         // Number of nonlinear variables - assume all variables are nonlinear
         problem_info_.num_nonlin_vars = problem_info_.n;
         problem_info_.has_nonlinear_vars = true;

         // 7. Split constraints for CONOPT (constraints with both bounds need to be duplicated)
         problem_info_.split_constraints();

         // Verify we have complete information
         if (!problem_info_.is_complete()) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                             "CONOPT Shim: Incomplete problem information retrieved.\n");
            }
            return Invalid_Problem_Definition;
         }

         if (!IsNull(jnlst_)) {
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
         if (!cntvect_ || IsNull(tnlp) || IsNull(jnlst_) || IsNull(stats_)) {
            return Invalid_Problem_Definition;
         }

         // --- Retrieve problem information from TNLP ---
         ApplicationReturnStatus info_status = RetrieveProblemInfo(tnlp);
         if (info_status != Solve_Succeeded) {
            return info_status;
         }

         // Start timing
         if (!IsNull(stats_)) {
            stats_->StartTiming();
         }

         // --- Prepare the context cookie ---
         context_.tnlp_ = GetRawPtr(tnlp);
         context_.journalist_ = GetRawPtr(jnlst_);
         context_.stats_ = GetRawPtr(stats_);
         context_.problem_info_ = &problem_info_; // Add problem info to context

        // --- Initialize FDEval cache ---
        context_.fdeval_cache_ = new FDEvalCache(problem_info_.m_split, problem_info_.nnz_jac_g, problem_info_.n);

         // creating the status solution structure
         context_.status_solution_ = new ConoptStatusSolution();

         COIDEF_UsrMem(cntvect_, &context_);

         // --- Register all our trampolines ---
         // These C functions are defined in IpoptToConoptCallbacks.cpp
         // and will call the virtual methods on the tnlp* cookie.

         int COI_ERROR = 0;

         // Mandatory Callbacks
         COI_ERROR += COIDEF_ReadMatrix(cntvect_, Conopt_ReadMatrix);
         COI_ERROR += COIDEF_FDEval(cntvect_, Conopt_FDEval);
         COI_ERROR += COIDEF_Status(cntvect_, Conopt_Status);
         COI_ERROR += COIDEF_Solution(cntvect_, Conopt_Solution);
         COI_ERROR += COIDEF_Message(cntvect_, Conopt_Message);
         COI_ERROR += COIDEF_ErrMsg(cntvect_, Conopt_ErrMsg);

         // Optional Callbacks
         COI_ERROR += COIDEF_FDEvalIni(cntvect_, Conopt_FDEvalIni);
         COI_ERROR += COIDEF_FDEvalEnd(cntvect_, Conopt_FDEvalEnd);

         // Check for callback registration errors
         if (COI_ERROR != 0) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                            "CONOPT Shim Error: Callback registration failed with error %d\n", COI_ERROR);
            }
            return Internal_Error;
         }

         // ... (Register all others: Progress, SDDir, etc.) ...

         // --- Set up CONOPT problem information ---
         if (!SetupConoptProblem()) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                            "CONOPT Shim Error: Failed to set up CONOPT problem information.\n");
            }
            return Internal_Error;
         }

         // --- Solve the problem ---
         COI_ERROR = COI_Solve(cntvect_);

         // Check for CONOPT solve errors
         if (COI_ERROR != 0) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                            "CONOPT Shim Error: COI_Solve failed with error %d\n", COI_ERROR);
            }
            return Internal_Error;
         }

         // Stop timing
         if (!IsNull(stats_)) {
            stats_->StopTiming();
         }

         // reporting the final solution to the user
         CallFinalizeSolutionWithCachedData(&context_);

         // --- Populate SolveStatistics with CONOPT data ---
         PopulateSolveStatistics(&context_);

         // --- Cleanup StatusSolution object ---
         if (context_.status_solution_) {
            delete context_.status_solution_;
            context_.status_solution_ = nullptr;
         }

         // --- Cleanup FDEval cache ---
         if (context_.fdeval_cache_) {
            delete context_.fdeval_cache_;
            context_.fdeval_cache_ = nullptr;
         }

         // Return the solve status from the statistics
         if (!IsNull(stats_)) {
            return stats_->SolveStatus();
         }

         return Solve_Succeeded;
      }

      /**
       * @brief Getter for statistics (stubbed)
       */
      SmartPtr<SolveStatistics> Statistics() {
         return stats_;
      }


      /**
       * @brief Get options object (stubbed for now)
       */
      SmartPtr<OptionsList> Options() {
         // Return a dummy options object for now
         static SmartPtr<OptionsList> dummy_options = new OptionsList();
         return dummy_options;
      }

      /**
       * @brief Set up CONOPT problem information
       * This method configures CONOPT with all the necessary problem parameters
       * @return true if successful, false otherwise
       */
      bool SetupConoptProblem() {
         if (!cntvect_) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                            "CONOPT Shim Error: CONOPT handle is NULL in SetupConoptProblem.\n");
            }
            return false;
         }

         // Set problem dimensions
         int COI_ERROR = 0;

         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                         "CONOPT Shim: Setting problem dimensions: n=%d, m_split=%d, nnz_jac_g_split=%d, nnz_h_lag=%d\n",
                         problem_info_.n, problem_info_.m_split, problem_info_.nnz_jac_g_split,
                         problem_info_.nnz_h_lag);
         }

        COI_ERROR += COIDEF_NumVar(cntvect_, problem_info_.n);
        COI_ERROR += COIDEF_NumCon(cntvect_, problem_info_.m_split);
        COI_ERROR += COIDEF_NumNz(cntvect_, problem_info_.nnz_jac_g_split);
        COI_ERROR += COIDEF_NumNlNz(cntvect_, problem_info_.nnz_jac_g_split);
        COI_ERROR += COIDEF_NumHess(cntvect_, problem_info_.nnz_h_lag);

        assert(problem_info_.objective_row_index >= 0);
        COI_ERROR += COIDEF_ObjCon(cntvect_, problem_info_.objective_row_index);

        // Set output control
        COI_ERROR += COIDEF_StdOut(cntvect_, 1);

        COI_ERROR += COIDEF_FVforAll(cntvect_, 1);

         // Check for setup errors
         if (COI_ERROR != 0) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                            "CONOPT Shim Error: Problem setup failed with total error %d\n", COI_ERROR);
            }
            return false;
         }

         return true;
      }

      // ... (All other public IpoptApplication methods) ...
   };

// Factory function
SmartPtr<IpoptApplication> IpoptApplicationFactory();


// Factory function implementation
inline SmartPtr<IpoptApplication> IpoptApplicationFactory() {
   return new IpoptApplication();
}

} // namespace Ipopt
