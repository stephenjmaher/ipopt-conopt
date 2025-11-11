#pragma once

/*
 * This is your drop-in replacement for Ipopt/IpoptApplication.hpp
 */

/*  1. INCLUDE ORIGINAL IPOPT HEADERS (for utilities) */
#include "IpSmartPtr.hpp"
#include "IpReferenced.hpp"
#include "IpReturnCodes.hpp"
#include "IpSolveStatistics.hpp"
#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
#include "IpOptionsList.hpp"
#include "Ipopt/IpIpoptData.hpp"
/*  ... etc ... */

/*  2. INCLUDE THE CONOPT C-API */
#include "conopt.h"

/*  3. INCLUDE OUR NEW TRAMPOLINE DECLARATIONS */
#include "IpoptToConoptCallbacks.hpp"
#include "IpoptProblemInfo.hpp"
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <string>
#include <map>

/*  5. FORWARD DECLARE THE OTHER SHIM CLASSES */
namespace Ipopt {
class TNLP;
class Journal;
class Journalist;
class RegisteredOptions;
class SolveStatistics;
}  // namespace Ipopt

namespace Ipopt {

/*  Forward declare our shim class */
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

   /*  We'll also need to store the user's TNLP */
   SmartPtr<TNLP> tnlp_;

   SmartPtr<SolveStatistics> stats_;

   /**
    * @brief IpoptData instance for callbacks
    */
   SmartPtr<IpoptData> ip_data_;

   /**
    * @brief Problem information retrieved from TNLP
    */
   IpoptProblemInfo problem_info_;

   /**
    * @brief OptionsList that maps Ipopt options to CONOPT
    */
   SmartPtr<OptionsList> options_;

   void InitializeIpoptApplication(bool create_console_out) {
      COI_Create(&cntvect_);
      stats_ = new SolveStatistics();
      if (IsNull(options_))
         options_ = new OptionsList();

      /*  Set infinity value from OptionsList - must be set */
      Number infinity_value;
      if (IsNull(options_) ||
            !options_->GetNumericValue("nlp_upper_bound_inf", infinity_value, "")) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: nlp_upper_bound_inf option is required but not set.\n");
         }
         /*  Cannot return error from constructor, but will fail when OptimizeTNLP is called */
         /*  Leave upper_bound_inf uninitialized (0.0) - will cause error in OptimizeTNLP */
      }
      else {
         problem_info_.upper_bound_inf = infinity_value;
      }

      /*  Set verbose output for debugging */
      if (IsNull(jnlst_)) {
         jnlst_ = new Journalist();
      }
   }

 public:
   /**
    * @brief Constructor: Create the CONOPT handle.
    */
   IpoptApplication(bool create_console_out = true, bool create_empty = false)
       : cntvect_(nullptr), jnlst_(nullptr), stats_(nullptr), options_(nullptr) {
      InitializeIpoptApplication(create_console_out);
   }

   /*  Constructor accepting journalist */
   IpoptApplication(SmartPtr<Journalist> jnlst)
       : cntvect_(nullptr), jnlst_(jnlst), stats_(nullptr), options_(nullptr) {
      InitializeIpoptApplication(true);
   }

   IpoptApplication(SmartPtr<RegisteredOptions> reg_options, SmartPtr<OptionsList> options,
         SmartPtr<Journalist> jnlst)
       : cntvect_(nullptr), jnlst_(jnlst), stats_(nullptr), options_(options) {
      InitializeIpoptApplication(true);
   }

   /**
    * @brief Destructor: Free the CONOPT handle.
    */
   virtual ~IpoptApplication() {
      if (cntvect_) {
         COI_Free(&cntvect_);
      }
   }

   void InitializeDefaultJournal() {
      if (!IsNull(jnlst_)) {
         /*  Use summary-level output by default to reduce verbosity */
         jnlst_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
      }
   }

   ApplicationReturnStatus Initialize(bool allow_clobber = false) {
      InitializeDefaultJournal();

      // Try to open the default options file
      std::ifstream is("ipopt.opt");
      if (!is.is_open()) {
         // This is not an error, just no file to read.
         // (Ipopt behavior: only returns true if file *was* read and parsed OK)
         return Solve_Succeeded;  // Or maybe just return, no file found
      }
      return ParseOptionsStream(is, allow_clobber, GetRawPtr(jnlst_));
   }

   /**
    * @brief Reads options from a user-specified file.
    */
   ApplicationReturnStatus Initialize(const std::string& params_file, bool allow_clobber = false) {
      InitializeDefaultJournal();

      std::ifstream is(params_file);
      if (!is.is_open()) {
         if (!IsNull(jnlst_))
            jnlst_->Printf(
                  J_ERROR, J_MAIN, "Error: Could not open options file: %s\n", params_file.c_str());
         // This is the correct Ipopt return code for a file not found
         return Invalid_Option;
      }
      return ParseOptionsStream(is, allow_clobber, GetRawPtr(jnlst_));
   }

   /**
    * @brief Reads options from any std::istream.
    */
   ApplicationReturnStatus Initialize(std::istream& is, bool allow_clobber = false) {
      InitializeDefaultJournal();

      // This is the base implementation that the file-based methods call
      return ParseOptionsStream(is, allow_clobber, GetRawPtr(jnlst_));
   }

   /** Method to register all Ipopt options. */
   static void RegisterAllIpoptOptions(const SmartPtr<RegisteredOptions>& roptions) {}

   /**
    * @brief Retrieve problem information from the TNLP object
    * This method calls all the necessary TNLP methods to gather problem data
    * @param tnlp The TNLP object to retrieve information from
    * @param infinity The infinity value to use for determining problem type
    */
   ApplicationReturnStatus RetrieveProblemInfo(SmartPtr<TNLP> tnlp, Number infinity) {
      if (IsNull(tnlp)) {
         return Invalid_Problem_Definition;
      }

      /*  Clear any existing data, passing infinity value from options */
      problem_info_.clear(infinity);

      /*  1. Get basic problem dimensions */
      TNLP::IndexStyleEnum index_style;
      if (!tnlp->get_nlp_info(problem_info_.n, problem_info_.m, problem_info_.nnz_jac_g,
                problem_info_.nnz_h_lag, index_style)) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim: Failed to get NLP info from TNLP.\n");
         }
         return Invalid_Problem_Definition;
      }

      /*  Store the index style - we'll handle FORTRAN by converting indices in the shim */
      problem_info_.index_style = (index_style == TNLP::FORTRAN_STYLE) ? FORTRAN_STYLE : C_STYLE;

      if (problem_info_.index_style == FORTRAN_STYLE && !IsNull(jnlst_)) {
         jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Shim: FORTRAN-style indexing detected. "
               "Converting row/column indices to C-style for CONOPT.\n");
      }

      /*  Resize vectors based on dimensions */
      problem_info_.resize_vectors();

      /*  2. Get bounds information */
      if (!tnlp->get_bounds_info(problem_info_.n, problem_info_.x_l.data(),
                problem_info_.x_u.data(), problem_info_.m, problem_info_.g_l.data(),
                problem_info_.g_u.data())) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim: Failed to get bounds info from TNLP.\n");
         }
         return Invalid_Problem_Definition;
      }

      /*  Clamp variable bounds to [-upper_bound_inf, upper_bound_inf] */
      for (Index i = 0; i < problem_info_.n; ++i) {
         if (IsFiniteNumber(problem_info_.x_l[i]) && problem_info_.x_l[i] < -infinity) {
            problem_info_.x_l[i] = -infinity;
         }
         if (IsFiniteNumber(problem_info_.x_u[i]) && problem_info_.x_u[i] > infinity) {
            problem_info_.x_u[i] = infinity;
         }
      }

      /*  3. Get starting point */
      if (!tnlp->get_starting_point(problem_info_.n, problem_info_.init_x_req,
                problem_info_.x_init.data(), problem_info_.init_z_req,
                problem_info_.z_L_init.data(), problem_info_.z_U_init.data(), problem_info_.m,
                problem_info_.init_lambda_req, problem_info_.lambda_init.data())) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim: Failed to get starting point from TNLP.\n");
         }
         return Invalid_Problem_Definition;
      }

      /*  3.5. Get scaling parameters
       *  Prepare arrays for scaling factors
       *
       *  TODO: apply the scaling to the problem
       */
      std::vector<Number> x_scaling_temp(problem_info_.n);
      std::vector<Number> g_scaling_temp(problem_info_.m);
      Number obj_scaling_temp = 1.0;
      bool use_x_scaling_temp = false;
      bool use_g_scaling_temp = false;

      /*  Call get_scaling_parameters - it may return false if not implemented, which is OK */
      if (tnlp->get_scaling_parameters(obj_scaling_temp, use_x_scaling_temp, problem_info_.n,
                x_scaling_temp.data(), use_g_scaling_temp, problem_info_.m,
                g_scaling_temp.data())) {
         /*  Scaling parameters were provided - store them */
         problem_info_.obj_scaling = obj_scaling_temp;
         problem_info_.use_x_scaling = use_x_scaling_temp;
         problem_info_.use_g_scaling = use_g_scaling_temp;

         if (use_x_scaling_temp) {
            problem_info_.x_scaling = x_scaling_temp;
         }
         if (use_g_scaling_temp) {
            problem_info_.g_scaling = g_scaling_temp;
         }

         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                  "CONOPT Shim: Scaling parameters retrieved (obj_scaling=%g, use_x_scaling=%s, \n"
                  "use_g_scaling=%s).\n",
                  obj_scaling_temp, use_x_scaling_temp ? "true" : "false",
                  use_g_scaling_temp ? "true" : "false");
         }
      }
      else {
         /*  Scaling parameters not provided - use defaults */
         problem_info_.obj_scaling = 1.0;
         problem_info_.use_x_scaling = false;
         problem_info_.use_g_scaling = false;
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                  "CONOPT Shim: No scaling parameters provided, using defaults.\n");
         }
      }

      /*  4. Get Jacobian structure (if needed) */
      if (problem_info_.nnz_jac_g > 0) {
         if (!tnlp->eval_jac_g(problem_info_.n, problem_info_.x_init.data(), true, problem_info_.m,
                   problem_info_.nnz_jac_g, problem_info_.jac_g_iRow.data(),
                   problem_info_.jac_g_jCol.data(), nullptr)) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                     "CONOPT Shim: Failed to get Jacobian structure from TNLP.\n");
            }
            return Invalid_Problem_Definition;
         }

         /*  Convert FORTRAN indices (1-based) to C-style (0-based) if needed */
         if (problem_info_.index_style == FORTRAN_STYLE) {
            for (Index k = 0; k < problem_info_.nnz_jac_g; ++k) {
               problem_info_.jac_g_iRow[k] = problem_info_.jac_g_iRow[k] - 1;
               problem_info_.jac_g_jCol[k] = problem_info_.jac_g_jCol[k] - 1;
            }
         }

         /*  4.5. Get nonlinear terms in Jacobian (if available) */
         if (problem_info_.nnz_jac_g > 0) {
            /*  Allocate arrays for nonlinear term indices */
            std::vector<Index> nl_iRow(problem_info_.nnz_jac_g);
            std::vector<Index> nl_jCol(problem_info_.nnz_jac_g);
            Index n_nl_terms = 0;

            /*  Call get_nonlinear_terms - it may return false if not implemented */
            if (tnlp->get_nonlinear_terms(problem_info_.n, problem_info_.m, n_nl_terms,
                      nl_iRow.data(), nl_jCol.data())) {
               /*  Nonlinear terms were provided - store them */
               problem_info_.nonlinear_terms_collected = true;
               problem_info_.n_nl_terms = n_nl_terms;

               /*  Convert FORTRAN indices (1-based) to C-style (0-based) if needed */
               if (problem_info_.index_style == FORTRAN_STYLE) {
                  for (Index k = 0; k < n_nl_terms; ++k) {
                     nl_iRow[k] = nl_iRow[k] - 1;
                     nl_jCol[k] = nl_jCol[k] - 1;
                  }
               }

               /*  Create a mapping from (row, col) pairs to Jacobian entry indices */
               /*  First, build a map from (row, col) to Jacobian index */
               std::map<std::pair<Index, Index>, Index> jacobian_map;
               for (Index k = 0; k < problem_info_.nnz_jac_g; ++k) {
                  jacobian_map[std::make_pair(
                        problem_info_.jac_g_iRow[k], problem_info_.jac_g_jCol[k])] = k;
               }

               /*  Mark nonlinear entries in jac_g_is_nonlinear */
               for (Index k = 0; k < n_nl_terms; ++k) {
                  auto it = jacobian_map.find(std::make_pair(nl_iRow[k], nl_jCol[k]));
                  if (it != jacobian_map.end()) {
                     Index jac_idx = it->second;
                     if (jac_idx >= 0 && jac_idx < problem_info_.nnz_jac_g) {
                        problem_info_.jac_g_is_nonlinear[jac_idx] = true;
                     }
                  }
               }

               if (!IsNull(jnlst_)) {
                  jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                        "CONOPT Shim: Collected %d nonlinear terms from Jacobian.\n",
                        static_cast<int>(n_nl_terms));
               }
            }
            else {
               /*  Nonlinear terms not provided - mark as not collected */
               problem_info_.nonlinear_terms_collected = false;
               problem_info_.n_nl_terms = 0;
               if (!IsNull(jnlst_)) {
                  jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                        "CONOPT Shim: Nonlinear terms not provided by TNLP.\n");
               }
            }
         }
      }

      /*  5. Get Hessian structure (only if nonlinear terms were collected) */
      if (problem_info_.nnz_h_lag > 0) {
         if (problem_info_.nonlinear_terms_collected) {
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

            /*  Convert FORTRAN indices (1-based) to C-style (0-based) if needed */
            if (problem_info_.index_style == FORTRAN_STYLE) {
               for (Index k = 0; k < problem_info_.nnz_h_lag; ++k) {
                  problem_info_.hess_iRow[k] = problem_info_.hess_iRow[k] - 1;
                  problem_info_.hess_jCol[k] = problem_info_.hess_jCol[k] - 1;
               }
            }

            /*  Compute CONOPT ordering permutation once here */
            problem_info_.compute_hessian_permutation();

            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                     "CONOPT Shim: Collected Hessian structure (%d non-zeros).\n",
                     static_cast<int>(problem_info_.nnz_h_lag));
            }
         }
         else {
            /*  Nonlinear terms not collected - skip Hessian collection */
            /*  Set nnz_h_lag to 0 to indicate no Hessian will be used */
            problem_info_.nnz_h_lag = 0;
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
                     "CONOPT Shim: Skipping Hessian collection (nonlinear terms not collected).\n");
            }
         }
      }

      /*  6. Get optional information (these may not be implemented by all TNLPs) */

      /*  Variable linearity - populate with dummy data since we don't use it */
      if (problem_info_.n > 0) {
         std::vector<TNLP::LinearityType> var_linearity(problem_info_.n);
         tnlp->get_variables_linearity(problem_info_.n, var_linearity.data());
         /*  Fill with dummy data */
         for (Index i = 0; i < problem_info_.n; ++i) {
            problem_info_.var_linearity[i] = NONLINEAR; /*  Default to non-linear */
         }
      }

      /*  Number of nonlinear variables - assume all variables are nonlinear */
      problem_info_.num_nonlin_vars = problem_info_.n;
      problem_info_.has_nonlinear_vars = true;

      /*  7. Split constraints for CONOPT (constraints with both bounds need to be duplicated) */
      problem_info_.split_constraints();

      /*  Verify we have complete information */
      if (!problem_info_.is_complete()) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim: Incomplete problem information retrieved.\n");
         }
         return Invalid_Problem_Definition;
      }

      if (!IsNull(jnlst_)) {
         jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Shim: Problem info retrieved and constraints split successfully.\n%s",
               problem_info_.to_string().c_str());
      }

      return Solve_Succeeded;
   }

 private:
   /**
    * @brief Setup solve instances (stats, ip_data, fdeval_cache, status_solution) and context.
    * @param tnlp The TNLP object to use
    * @return true on success, false on error
    */
   bool SetupSolveInstances(SmartPtr<TNLP> tnlp) {
      /*  Create new statistics instance if needed */
      if (IsNull(stats_)) {
         stats_ = new SolveStatistics();
      }

      /*  Start timing */
      if (!IsNull(stats_)) {
         stats_->StartTiming();
      }

      /*  --- Prepare the context cookie --- */
      context_.tnlp_ = GetRawPtr(tnlp);
      context_.journalist_ = GetRawPtr(jnlst_);
      context_.stats_ = GetRawPtr(stats_);
      context_.problem_info_ = &problem_info_;
      context_.options_list_ = GetRawPtr(options_);

      /*  --- Create IpoptData instance --- */
      ip_data_ = new IpoptData();
      context_.ip_data_ = GetRawPtr(ip_data_);

      /*  --- Initialize FDEval cache --- */
      context_.fdeval_cache_ =
            new FDEvalCache(problem_info_.m_split, problem_info_.nnz_jac_g, problem_info_.n);

      /*  --- Create status solution structure --- */
      context_.status_solution_ = new ConoptStatusSolution();

      /*  --- Update USRMEM pointer with context --- */
      COIDEF_UsrMem(cntvect_, &context_);

      return true;
   }

   /**
    * @brief Clean up old solve instances (fdeval_cache, ip_data, status_solution).
    * Statistics are not cleaned up here as they may be accessed via Statistics().
    */
   void CleanupOldInstances() {
      /*  Free old fdeval_cache */
      if (context_.fdeval_cache_) {
         delete context_.fdeval_cache_;
         context_.fdeval_cache_ = nullptr;
      }

      /*  Free old ip_data */
      if (!IsNull(ip_data_)) {
         context_.ip_data_ = nullptr;
         ip_data_ = nullptr;
      }

      /*  Free old status_solution */
      if (context_.status_solution_) {
         delete context_.status_solution_;
         context_.status_solution_ = nullptr;
      }
   }

   /**
    * @brief Register CONOPT callbacks.
    * @return ApplicationReturnStatus indicating success or failure
    */
   ApplicationReturnStatus RegisterConoptCallbacks() {
      int COI_ERROR = 0;

      /*  Mandatory Callbacks */
      COI_ERROR += COIDEF_ReadMatrix(cntvect_, Conopt_ReadMatrix);
      COI_ERROR += COIDEF_FDEval(cntvect_, Conopt_FDEval);
      COI_ERROR += COIDEF_Status(cntvect_, Conopt_Status);
      COI_ERROR += COIDEF_Solution(cntvect_, Conopt_Solution);
      COI_ERROR += COIDEF_Message(cntvect_, Conopt_Message);
      COI_ERROR += COIDEF_ErrMsg(cntvect_, Conopt_ErrMsg);

      /*  Optional Callbacks */
      COI_ERROR += COIDEF_FDEvalIni(cntvect_, Conopt_FDEvalIni);
      COI_ERROR += COIDEF_Option(cntvect_, Conopt_Option);

      /*
       * Hessian of Lagrangian callbacks
       * these are only used if the hessian structure is provided AND nonlinear terms were
       * collected.
       */
      if (problem_info_.nnz_h_lag > 0 && problem_info_.nonlinear_terms_collected) {
         COI_ERROR += COIDEF_2DLagrStr(cntvect_, Conopt_2DLagrStr);
         COI_ERROR += COIDEF_2DLagrVal(cntvect_, Conopt_2DLagrVal);
      }

      /*  Check for callback registration errors */
      if (COI_ERROR != 0) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: Callback registration failed with error %d\n", COI_ERROR);
         }
         return Internal_Error;
      }

      return Solve_Succeeded;
   }

   /**
    * @brief Apply options to CONOPT and solve the problem.
    * @param error_context Context string for error messages
    * @return ApplicationReturnStatus indicating success or failure
    */
   ApplicationReturnStatus ApplyOptionsAndSolve(const char* error_context) {
      /*  --- Apply user options to CONOPT --- */
      if (!IsNull(options_)) {
         if (!options_->ApplyToConopt(cntvect_, GetRawPtr(jnlst_))) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                     "CONOPT Shim Error: Failed to apply options to CONOPT in %s.\n",
                     error_context);
            }
            return Internal_Error;
         }
      }

      /*  --- Solve the problem --- */
      int COI_ERROR = COI_Solve(cntvect_);

      /*  Check for CONOPT solve errors */
      if (COI_ERROR != 0) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: COI_Solve failed in %s with error %d\n", error_context,
                  COI_ERROR);
         }
         return Internal_Error;
      }

      return Solve_Succeeded;
   }

   /**
    * @brief Finalize the solve: stop timing, call finalize_solution, populate statistics, and
    * cleanup.
    * @return ApplicationReturnStatus indicating the solve result
    */
   ApplicationReturnStatus FinalizeAndCleanup() {
      /*  Stop timing */
      if (!IsNull(stats_)) {
         stats_->StopTiming();
      }

      /*  --- Report the final solution to the user --- */
      CallFinalizeSolutionWithCachedData(&context_);

      /*  --- Populate SolveStatistics with CONOPT data --- */
      PopulateSolveStatistics(&context_);

      /*  --- Cleanup temporary instances --- */
      CleanupOldInstances();

      /*  Return the solve status from the statistics */
      if (!IsNull(stats_)) {
         return stats_->SolveStatus();
      }

      return Solve_Succeeded;
   }

 public:
   /**
    * @brief The main "Solve" method.
    * This is where we wire everything together.
    */
   ApplicationReturnStatus OptimizeTNLP(SmartPtr<TNLP> tnlp) {
      if (!cntvect_ || IsNull(tnlp) || IsNull(jnlst_) || IsNull(stats_)) {
         return Invalid_Problem_Definition;
      }

      /*  --- Store the TNLP for later verification in ReOptimizeTNLP --- */
      tnlp_ = tnlp;

      /*  --- Retrieve infinity value from options - must be set --- */
      Number infinity_value;
      if (IsNull(options_) ||
            !options_->GetNumericValue("nlp_upper_bound_inf", infinity_value, "")) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: nlp_upper_bound_inf option is required but not set.\n");
         }
         return Invalid_Option;
      }

      /*  --- Retrieve problem information from TNLP --- */
      ApplicationReturnStatus info_status = RetrieveProblemInfo(tnlp, infinity_value);
      if (info_status != Solve_Succeeded) {
         return info_status;
      }

      /*  --- Setup solve instances and context --- */
      if (!SetupSolveInstances(tnlp)) {
         return Internal_Error;
      }

      /*  --- Register CONOPT callbacks --- */
      ApplicationReturnStatus callback_status = RegisterConoptCallbacks();
      if (callback_status != Solve_Succeeded) {
         return callback_status;
      }

      /*  --- Set up CONOPT problem information --- */
      if (!SetupConoptProblem()) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: Failed to set up CONOPT problem information.\n");
         }
         return Internal_Error;
      }

      /*  --- Apply options and solve --- */
      ApplicationReturnStatus solve_status = ApplyOptionsAndSolve("OptimizeTNLP");
      if (solve_status != Solve_Succeeded) {
         return solve_status;
      }

      /*  --- Finalize and cleanup --- */
      return FinalizeAndCleanup();
   }

   /**
    * @brief Reoptimize the problem without redoing all setup.
    * This method calls COI_Solve without the setup performed in OptimizeTNLP.
    * Options are reapplied to CONOPT before solving in case they were changed.
    * An OptimizeTNLP call must be performed before calling this method.
    * The SolutionStatistics, fdeval_cache, and ip_data are freed and recreated.
    * @return ApplicationReturnStatus indicating the solve result
    */
   ApplicationReturnStatus ReOptimizeTNLP(SmartPtr<TNLP> tnlp) {
      /*  Verify that OptimizeTNLP was called first */
      if (!problem_info_.is_complete()) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: ReOptimizeTNLP called before OptimizeTNLP. "
                  "Problem information is not available.\n");
         }
         return Invalid_Problem_Definition;
      }

      if (!cntvect_ || IsNull(tnlp) || IsNull(jnlst_)) {
         return Invalid_Problem_Definition;
      }

      /*  Verify that the TNLP matches the one used in OptimizeTNLP */
      if (tnlp != tnlp_) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: ReOptimizeTNLP called with different TNLP than "
                  "OptimizeTNLP.\n");
         }
         return Invalid_Problem_Definition;
      }

      /*  --- Update starting point from TNLP --- */
      if (!tnlp->get_starting_point(problem_info_.n, problem_info_.init_x_req,
                problem_info_.x_init.data(), problem_info_.init_z_req,
                problem_info_.z_L_init.data(), problem_info_.z_U_init.data(), problem_info_.m,
                problem_info_.init_lambda_req, problem_info_.lambda_init.data())) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim: Failed to get starting point from TNLP in ReOptimizeTNLP.\n");
         }
         return Invalid_Problem_Definition;
      }

      /*  --- Free old instances --- */
      if (!IsNull(stats_)) {
         stats_ = nullptr;
      }
      CleanupOldInstances();

      /*  --- Setup solve instances and context --- */
      if (!SetupSolveInstances(tnlp)) {
         return Internal_Error;
      }

      /*  --- Apply options and solve (no callback registration or problem setup needed) --- */
      ApplicationReturnStatus solve_status = ApplyOptionsAndSolve("ReOptimizeTNLP");
      if (solve_status != Solve_Succeeded) {
         return solve_status;
      }

      /*  --- Finalize and cleanup --- */
      return FinalizeAndCleanup();
   }

   /** Get the Journalist for printing output */
   virtual SmartPtr<Journalist> Jnlst() {
      return jnlst_;
   }

   /**
    * @brief Getter for statistics (stubbed)
    */
   SmartPtr<SolveStatistics> Statistics() {
      return stats_;
   }

   /**
    * @brief Get options object
    */
   SmartPtr<OptionsList> Options() {
      return options_;
   }

   /** Get the IpoptData Object */
   SmartPtr<IpoptData> IpoptDataObject() {
      return ip_data_;
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

      /*  Set problem dimensions */
      int COI_ERROR = 0;

      if (!IsNull(jnlst_)) {
         jnlst_->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Shim: Setting problem dimensions: n=%d, m_split=%d, nnz_jac_g_split=%d, "
               "nnz_h_lag=%d\n",
               problem_info_.n, problem_info_.m_split, problem_info_.nnz_jac_g_split,
               problem_info_.nnz_h_lag);
      }

      COI_ERROR += COIDEF_NumVar(cntvect_, problem_info_.n);
      COI_ERROR += COIDEF_NumCon(cntvect_, problem_info_.m_split);
      COI_ERROR += COIDEF_NumNz(cntvect_, problem_info_.nnz_jac_g_split);
      /*  Use split nonlinear count if nonlinear terms were collected, otherwise use total split
       * count */
      int num_nl_nz = problem_info_.nonlinear_terms_collected ? problem_info_.nnz_jac_g_split_nl
                                                              : problem_info_.nnz_jac_g_split;
      COI_ERROR += COIDEF_NumNlNz(cntvect_, num_nl_nz);
      COI_ERROR += COIDEF_NumHess(cntvect_, problem_info_.nnz_h_lag);

      assert(problem_info_.objective_row_index >= 0);
      COI_ERROR += COIDEF_ObjCon(cntvect_, problem_info_.objective_row_index);
      COI_ERROR += COIDEF_OptDir(cntvect_, -1);

      /*  Set output control */
      COI_ERROR += COIDEF_StdOut(cntvect_, 1);

      COI_ERROR += COIDEF_FVforAll(cntvect_, 1);

      /*  Check for setup errors */
      if (COI_ERROR != 0) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: Problem setup failed with total error %d\n", COI_ERROR);
         }
         return false;
      }

      return true;
   }

 private:
   // Helper to trim whitespace from start and end
   std::string trim(const std::string& str) {
      size_t first = str.find_first_not_of(" \t\n\r");
      if (std::string::npos == first)
         return "";
      size_t last = str.find_last_not_of(" \t\n\r");
      return str.substr(first, (last - first + 1));
   }

   // Helper to convert to lower-case for case-insensitive matching
   std::string to_lower(const std::string& s) {
      std::string out = s;
      std::transform(
            out.begin(), out.end(), out.begin(), [](unsigned char c) { return std::tolower(c); });
      return out;
   }

   /**
    * @brief Core logic for parsing an options stream (like ipopt.opt)
    */
   ApplicationReturnStatus ParseOptionsStream(
         std::istream& is, bool allow_clobber, Journalist* jnlst) {
      std::string line;
      int line_num = 0;

      while (std::getline(is, line)) {
         line_num++;

         // 1. Remove comments
         size_t comment_pos = line.find('#');
         if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
         }
         line = trim(line);
         if (line.empty())
            continue;

         // 2. Split line at first space or '='
         size_t eq_pos = line.find_first_of("= \t");
         if (eq_pos == std::string::npos || eq_pos == 0) {
            if (jnlst)
               jnlst->Printf(J_WARNING, J_MAIN, "Skipping malformed option line %d: %s\n", line_num,
                     line.c_str());
            continue;
         }

         std::string name = to_lower(trim(line.substr(0, eq_pos)));
         std::string value_str = trim(line.substr(eq_pos + 1));

         // 3. Check allow_clobber logic
         Index int_val;
         Number num_val;
         std::string str_val;
         bool is_set = options_->GetIntegerValue(name, int_val, "") ||
               options_->GetNumericValue(name, num_val, "") ||
               options_->GetStringValue(name, str_val, "");

         if (is_set && !allow_clobber) {
            if (jnlst)
               jnlst->Printf(J_ERROR, J_MAIN,
                     "Option %s on line %d conflicts with a pre-set value (allow_clobber=false).\n",
                     name.c_str(), line_num);
            return Invalid_Option;  // Fail loudly
         }

         // 4. Set the value by trying to parse its type
         // This replicates Ipopt's behavior of auto-detecting type
         std::stringstream ss(value_str);
         Index int_val_check;
         ss >> int_val_check;
         if (!ss.fail() && ss.eof()) {
            // Successfully parsed as an integer
            options_->SetIntegerValue(name, int_val_check);
         }
         else {
            ss.clear();
            ss.str(value_str);
            Number num_val_check;
            ss >> num_val_check;
            if (!ss.fail() && ss.eof()) {
               // Successfully parsed as a double
               options_->SetNumericValue(name, num_val_check);
            }
            else {
               // Treat it as a string
               options_->SetStringValue(name, value_str);
            }
         }
      }  // end while loop

      return Solve_Succeeded;
   }

   /*  ... (All other public IpoptApplication methods) ... */
};

/*  Factory function */
SmartPtr<IpoptApplication> IpoptApplicationFactory();

/*  Factory function implementation */
inline SmartPtr<IpoptApplication> IpoptApplicationFactory() {
   return new IpoptApplication();
}

} /*  namespace Ipopt */
