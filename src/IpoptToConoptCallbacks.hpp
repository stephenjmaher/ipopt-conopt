/**
 * @file IpoptToConoptCallbacks.hpp
 * @brief Declares the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */
#ifndef IPOPT_TO_CONOPT_CALLBACKS_HPP
#define IPOPT_TO_CONOPT_CALLBACKS_HPP

#include "conopt.h" // For COI_CALLCONV and C-API types
#include <vector>

// Forward declare Ipopt classes needed by the struct
namespace Ipopt {
   class TNLP;
   class Journalist;
   class SolveStatistics;
   struct IpoptProblemInfo;
}


/**
 * @brief Struct to hold cached status and solution data from CONOPT callbacks.
 * This stores the status information and solution values for later use in finalize_solution.
 */
struct ConoptStatusSolution {
   // Status information
   bool status_cached_;                   // Whether status has been cached
   int conopt_modsta_;                    // CONOPT model status
   int conopt_solsta_;                    // CONOPT solver status
   int conopt_iter_;                      // CONOPT iteration count
   double conopt_objval_;                 // CONOPT objective value

   // Solution data
   bool solution_cached_;                 // Whether solution has been cached
   std::vector<double> x_solution_;       // Final variable values
   std::vector<double> x_marginals_;      // Variable marginals (z_L/z_U)
   std::vector<int> x_basis_;             // Variable basis indicators
   std::vector<int> x_status_;            // Variable status indicators
   std::vector<double> y_solution_;       // Final constraint values
   std::vector<double> y_marginals_;      // Constraint marginals (lambda)
   std::vector<int> y_basis_;             // Constraint basis indicators
   std::vector<int> y_status_;            // Constraint status indicators

   /**
    * @brief Constructor for ConoptStatusSolution
    */
   ConoptStatusSolution()
      : status_cached_(false),
        conopt_modsta_(0),
        conopt_solsta_(0),
        conopt_iter_(0),
        conopt_objval_(0.0),
        solution_cached_(false) {
   }

   /**
    * @brief Reset all cached data
    */
   void reset() {
      status_cached_ = false;
      solution_cached_ = false;
      conopt_modsta_ = 0;
      conopt_solsta_ = 0;
      conopt_iter_ = 0;
      conopt_objval_ = 0.0;
      x_solution_.clear();
      x_marginals_.clear();
      x_basis_.clear();
      x_status_.clear();
      y_solution_.clear();
      y_marginals_.clear();
      y_basis_.clear();
      y_status_.clear();
   }
};

/**
 * @brief Struct to hold cached constraint values and jacobian for FDEvalIni optimization.
 * This stores the results of constraint evaluations and jacobian with constant lookup time.
 */
struct FDEvalCache {
   std::vector<double> constraint_values_;  // Cached constraint values (size = numcons)
   std::vector<bool> constraint_valid_;     // Validity flags for each constraint (size = numcons)
   double objective_value_;                // Cached objective value
   bool objective_valid_;                  // Whether objective value is valid
   std::vector<double> objective_gradient_; // Cached objective gradient (size = n)
   bool objective_gradient_valid_;         // Whether objective gradient is valid
   std::vector<double> jacobian_values_;   // Cached jacobian values (size = nnz_jac_g)
   std::vector<bool> jacobian_valid_;      // Validity flags for jacobian entries
   bool jacobian_cached_;                  // Whether jacobian has been cached
   int num_constraints_;                   // Number of constraints (for bounds checking)
   int nnz_jacobian_;                      // Number of non-zero jacobian entries
   int num_variables_;                     // Number of variables (for objective gradient)

   /**
    * @brief Constructor for FDEvalCache
    * @param num_constraints Number of constraints to allocate space for
    * @param nnz_jacobian Number of non-zero jacobian entries to allocate space for
    * @param num_variables Number of variables (for objective gradient)
    */
   FDEvalCache(int num_constraints, int nnz_jacobian, int num_variables)
      : num_constraints_(num_constraints),
        nnz_jacobian_(nnz_jacobian),
        num_variables_(num_variables),
        objective_value_(0.0),
        objective_valid_(false),
        objective_gradient_valid_(false),
        jacobian_cached_(false) {
      constraint_values_.resize(num_constraints, 0.0);
      constraint_valid_.resize(num_constraints, false);
      objective_gradient_.resize(num_variables, 0.0);
      jacobian_values_.resize(nnz_jacobian, 0.0);
      jacobian_valid_.resize(nnz_jacobian, false);
   }

   /**
    * @brief Reset all validity flags to false
    */
   void invalidateAll() {
      std::fill(constraint_valid_.begin(), constraint_valid_.end(), false);
      std::fill(jacobian_valid_.begin(), jacobian_valid_.end(), false);
      objective_valid_ = false;
      objective_gradient_valid_ = false;
      jacobian_cached_ = false;
   }

   /**
    * @brief Cache jacobian values
    * @param jacobian_values Vector containing jacobian values
    */
   void cacheJacobian(const std::vector<double>& jacobian_values) {
      if (jacobian_values.size() == static_cast<size_t>(nnz_jacobian_)) {
         jacobian_values_ = jacobian_values;
         std::fill(jacobian_valid_.begin(), jacobian_valid_.end(), true);
         jacobian_cached_ = true;
      }
   }

   /**
    * @brief Get cached jacobian value for a given index
    * @param jacobian_idx The jacobian index to get the value for
    * @param value Output parameter for the cached value
    * @return true if the value is valid and cached, false otherwise
    */
   bool getCachedJacobianValue(int jacobian_idx, double& value) const {
      if (jacobian_idx < 0 || jacobian_idx >= nnz_jacobian_) {
         return false;
      }
      if (jacobian_valid_[jacobian_idx]) {
         value = jacobian_values_[jacobian_idx];
         return true;
      }
      return false;
   }

   /**
    * @brief Check if jacobian has been cached
    * @return true if jacobian is cached, false otherwise
    */
   bool isJacobianCached() const {
      return jacobian_cached_;
   }

   /**
    * @brief Cache objective gradient values
    * @param gradient_values Vector containing objective gradient values
    */
   void cacheObjectiveGradient(const std::vector<double>& gradient_values) {
      if (gradient_values.size() == static_cast<size_t>(num_variables_)) {
         objective_gradient_ = gradient_values;
         objective_gradient_valid_ = true;
      }
   }

   /**
    * @brief Get cached objective gradient value for a given variable index
    * @param var_idx The variable index to get the gradient for
    * @param value Output parameter for the cached value
    * @return true if the value is valid and cached, false otherwise
    */
   bool getCachedObjectiveGradientValue(int var_idx, double& value) const {
      if (var_idx < 0 || var_idx >= num_variables_) {
         return false;
      }
      if (objective_gradient_valid_) {
         value = objective_gradient_[var_idx];
         return true;
      }
      return false;
   }

   /**
    * @brief Check if objective gradient has been cached
    * @return true if objective gradient is cached, false otherwise
    */
   bool isObjectiveGradientCached() const {
      return objective_gradient_valid_;
   }
};

/**
 * @brief Struct to hold the context needed by the C trampolines.
 * This will be passed as the void* USRMEM cookie.
 */
typedef struct {
   Ipopt::TNLP* tnlp_;
   Ipopt::Journalist* journalist_;
   Ipopt::SolveStatistics* stats_;
   Ipopt::IpoptProblemInfo* problem_info_;
   FDEvalCache* fdeval_cache_;            // Cache for FDEvalIni optimization
   ConoptStatusSolution* status_solution_; // Cache for status and solution data
} IpoptConoptContext;

/**
 * @brief Cleanup function for IpoptConoptContext.
 * Frees any allocated memory in the context.
 */
void CleanupIpoptConoptContext(IpoptConoptContext* context);

/**
 * @brief Get cached constraint value for a given row.
 * @param context The context containing the cache
 * @param row_idx The row index to get the value for
 * @param value Output parameter for the cached value
 * @return true if the value is valid and cached, false otherwise
 */
bool GetCachedConstraintValue(IpoptConoptContext* context, int row_idx, double& value);

/**
 * @brief Get cached objective value.
 * @param context The context containing the cache
 * @param value Output parameter for the cached objective value
 * @return true if the objective value is valid and cached, false otherwise
 */
bool GetCachedObjectiveValue(IpoptConoptContext* context, double& value);

/**
 * @brief Get cached jacobian value for a given index.
 * @param context The context containing the cache
 * @param jacobian_idx The jacobian index to get the value for
 * @param value Output parameter for the cached value
 * @return true if the value is valid and cached, false otherwise
 */
bool GetCachedJacobianValue(IpoptConoptContext* context, int jacobian_idx, double& value);

/**
 * @brief Check if jacobian has been cached.
 * @param context The context containing the cache
 * @return true if jacobian is cached, false otherwise
 */
bool IsJacobianCached(IpoptConoptContext* context);

/**
 * @brief Get cached objective gradient value for a given variable index.
 * @param context The context containing the cache
 * @param var_idx The variable index to get the gradient for
 * @param value Output parameter for the cached value
 * @return true if the value is valid and cached, false otherwise
 */
bool GetCachedObjectiveGradientValue(IpoptConoptContext* context, int var_idx, double& value);

/**
 * @brief Check if objective gradient has been cached.
 * @param context The context containing the cache
 * @return true if objective gradient is cached, false otherwise
 */
bool IsObjectiveGradientCached(IpoptConoptContext* context);

/**
 * @brief Call finalize_solution with cached status and solution data.
 * This should be called at the end of OptimizeTNLP after both Status and Solution
 * callbacks have been invoked by CONOPT.
 * @param context The context containing the cached data
 * @return true on success, false on error
 */
bool CallFinalizeSolutionWithCachedData(IpoptConoptContext* context);

/**
 * @brief Populate SolveStatistics object with data from CONOPT.
 * This extracts information from the cached status and solution data
 * and populates the Ipopt SolveStatistics object.
 * @param context The context containing the cached data and stats object
 * @return true on success, false on error
 */
bool PopulateSolveStatistics(IpoptConoptContext* context);

/*
 * These are the C-style trampoline functions.
 * Their implementations will cast the void* USRMEM to an Ipopt::TNLP*
 * and call the corresponding Ipopt virtual method.
 */

#ifdef __cplusplus
extern "C" {
#endif

int COI_CALLCONV Conopt_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[], int TYPE[],
   double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[], int NUMVAR, int NUMCON,
   int NUMNZ, void *USRMEM);

int COI_CALLCONV Conopt_FDEval(const double X[], double *G, double JAC[], int ROWNO, const int JACNUM[], int MODE,
   int IGNERR, int *ERRCNT, int NUMVAR, int NUMJAC, int THREAD, void *USRMEM);

int COI_CALLCONV Conopt_FDEvalIni(const double X[], const int ROWLIST[], int MODE, int LISTSIZE, int NUMTHREAD,
   int IGNERR, int *ERRCNT, int NUMVAR, void *USRMEM);

int COI_CALLCONV Conopt_FDEvalEnd(int IGNERR, int *ERRCNT, void *USRMEM);

int COI_CALLCONV Conopt_Status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void *USRMEM);

int COI_CALLCONV Conopt_Solution(const double XVAL[], const double XMAR[], const int XBAS[], const int XSTA[],
   const double YVAL[], const double YMAR[], const int YBAS[], const int YSTA[], int NUMVAR, int NUMCON, void *USRMEM);

int COI_CALLCONV Conopt_Message(int SMSG, int DMSG, int NMSG, char *MSGV[], void *USRMEM);

int COI_CALLCONV Conopt_ErrMsg(int ROWNO, int COLNO, int POSNO, const char *MSG, void *USRMEM);

// ... and so on for all other callbacks (Progress, Option, TriOrd, SDDir, etc.) ...

#ifdef __cplusplus
}
#endif

#endif // IPOPT_TO_CONOPT_CALLBACKS_HPP
