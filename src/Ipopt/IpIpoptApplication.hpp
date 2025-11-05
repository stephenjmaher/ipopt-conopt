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

      /*  Set verbose output for debugging */
      if (IsNull(jnlst_)) {
         jnlst_ = new Journalist();

         /*  Use summary-level output by default to reduce verbosity */
         jnlst_->AddFileJournal("console", "stdout", Ipopt::J_SUMMARY);
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

   ApplicationReturnStatus Initialize(bool allow_clobber = false) {
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
      // This is the base implementation that the file-based methods call
      return ParseOptionsStream(is, allow_clobber, GetRawPtr(jnlst_));
   }

   /** Method to register all Ipopt options. */
   static void RegisterAllIpoptOptions(const SmartPtr<RegisteredOptions>& roptions) {}

   /**
    * @brief Retrieve problem information from the TNLP object
    * This method calls all the necessary TNLP methods to gather problem data
    */
   ApplicationReturnStatus RetrieveProblemInfo(SmartPtr<TNLP> tnlp) {
      if (IsNull(tnlp)) {
         return Invalid_Problem_Definition;
      }

      /*  Clear any existing data */
      problem_info_.clear();

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
         jnlst_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN,
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
      }

      /*  5. Get Hessian structure (if needed) */
      if (problem_info_.nnz_h_lag > 0) {
         if (!tnlp->eval_h(problem_info_.n, problem_info_.x_init.data(), true, 1.0, problem_info_.m,
                   problem_info_.lambda_init.data(), true, problem_info_.nnz_h_lag,
                   problem_info_.hess_iRow.data(), problem_info_.hess_jCol.data(), nullptr)) {
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

      /*  --- Retrieve problem information from TNLP --- */
      ApplicationReturnStatus info_status = RetrieveProblemInfo(tnlp);
      if (info_status != Solve_Succeeded) {
         return info_status;
      }

      /*  Start timing */
      if (!IsNull(stats_)) {
         stats_->StartTiming();
      }

      /*  --- Prepare the context cookie --- */
      context_.tnlp_ = GetRawPtr(tnlp);
      context_.journalist_ = GetRawPtr(jnlst_);
      context_.stats_ = GetRawPtr(stats_);
      context_.problem_info_ = &problem_info_;      /*  Add problem info to context */
      context_.options_list_ = GetRawPtr(options_); /*  Add options list to context */

      /*  --- Create IpoptData instance --- */
      ip_data_ = new IpoptData();
      context_.ip_data_ = GetRawPtr(ip_data_);

      /*  --- Initialize FDEval cache --- */
      context_.fdeval_cache_ =
            new FDEvalCache(problem_info_.m_split, problem_info_.nnz_jac_g, problem_info_.n);

      /*  creating the status solution structure */
      context_.status_solution_ = new ConoptStatusSolution();

      COIDEF_UsrMem(cntvect_, &context_);

      /*
       * --- Register all our trampolines ---
       * These C functions are defined in IpoptToConoptCallbacks.cpp
       * and will call the virtual methods on the tnlp* cookie.
       */

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
       * these are only used if the hessian structure is provided.
       */
      if (problem_info_.nnz_h_lag > 0) {
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

      /*  ... (Register all others: SDDir, etc.) ... */

      /*  --- Apply user options to CONOPT --- */
      if (!IsNull(options_)) {
         if (!options_->ApplyToConopt(cntvect_, GetRawPtr(jnlst_))) {
            if (!IsNull(jnlst_)) {
               jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                     "CONOPT Shim Error: Failed to apply options to CONOPT.\n");
            }
            return Internal_Error;
         }
      }

      /*  --- Set up CONOPT problem information --- */
      if (!SetupConoptProblem()) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: Failed to set up CONOPT problem information.\n");
         }
         return Internal_Error;
      }

      /*  --- Solve the problem --- */
      COI_ERROR = COI_Solve(cntvect_);

      /*  Check for CONOPT solve errors */
      if (COI_ERROR != 0) {
         if (!IsNull(jnlst_)) {
            jnlst_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN,
                  "CONOPT Shim Error: COI_Solve failed with error %d\n", COI_ERROR);
         }
         return Internal_Error;
      }

      /*  Stop timing */
      if (!IsNull(stats_)) {
         stats_->StopTiming();
      }

      /*  reporting the final solution to the user */
      CallFinalizeSolutionWithCachedData(&context_);

      /*  --- Populate SolveStatistics with CONOPT data --- */
      PopulateSolveStatistics(&context_);

      /*  --- Cleanup StatusSolution object --- */
      if (context_.status_solution_) {
         delete context_.status_solution_;
         context_.status_solution_ = nullptr;
      }

      /*  --- Cleanup FDEval cache --- */
      if (context_.fdeval_cache_) {
         delete context_.fdeval_cache_;
         context_.fdeval_cache_ = nullptr;
      }

      /*  --- Cleanup IpoptData --- */
      context_.ip_data_ = nullptr;
      ip_data_ = nullptr;

      /*  Return the solve status from the statistics */
      if (!IsNull(stats_)) {
         return stats_->SolveStatus();
      }

      return Solve_Succeeded;
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
      COI_ERROR += COIDEF_NumNlNz(cntvect_, problem_info_.nnz_jac_g_split);
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
