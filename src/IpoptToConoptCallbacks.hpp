/**
 * @file IpoptToConoptCallbacks.hpp
 * @brief Declares the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */
#ifndef IPOPT_TO_CONOPT_CALLBACKS_HPP
#define IPOPT_TO_CONOPT_CALLBACKS_HPP

#include "conopt.h" /*  For COI_CALLCONV and C-API types */
#include <cassert>
#include <vector>
#include <cstddef>
#include <utility>

/*  Forward declare Ipopt classes needed by the struct */
namespace Ipopt {
class TNLP;
class Journalist;
class SolveStatistics;
class OptionsList;
class IpoptData;
struct IpoptProblemInfo;
} /*  namespace Ipopt */

/**
 * @brief Struct to hold cached status and solution data from CONOPT callbacks.
 * This stores the status information and solution values for later use in finalize_solution.
 */
struct ConoptStatusSolution {
   /*  Status information */
   bool status_cached_;   /*  Whether status has been cached */
   int conopt_modsta_;    /*  CONOPT model status */
   int conopt_solsta_;    /*  CONOPT solver status */
   int conopt_iter_;      /*  CONOPT iteration count */
   double conopt_objval_; /*  Objective value in the TNLP's original (unscaled) units */
   double conopt_objval_scaled_; /*  Raw OBJVAL as reported by CONOPT, i.e. scaled by the
                                    obj_scaling magnitude applied in FDEval */

   /*  Solution data */
   bool solution_cached_;            /*  Whether solution has been cached */
   std::vector<double> x_solution_;  /*  Final variable values */
   std::vector<double> x_marginals_; /*  Variable marginals (z_L/z_U) */
   std::vector<int> x_basis_;        /*  Variable basis indicators */
   std::vector<int> x_status_;       /*  Variable status indicators */
   std::vector<double> y_solution_;  /*  Final constraint values */
   std::vector<double> y_marginals_; /*  Constraint marginals (lambda) */
   std::vector<int> y_basis_;        /*  Constraint basis indicators */
   std::vector<int> y_status_;       /*  Constraint status indicators */

   /**
    * @brief Constructor for ConoptStatusSolution
    */
   ConoptStatusSolution()
       : status_cached_(false), conopt_modsta_(0), conopt_solsta_(0), conopt_iter_(0),
         conopt_objval_(0.0), conopt_objval_scaled_(0.0), solution_cached_(false) {}

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
      conopt_objval_scaled_ = 0.0;
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
   std::vector<double> constraint_values_; /*  Cached constraint values, indexed by ORIGINAL
                                              constraint index (size = m, not m_split) - the
                                              split-row translation happens on demand in
                                              FDEval via original_constraint_map */
   std::vector<double> jacobian_values_;   /*  Cached jacobian values (size = nnz_jac_g) */
   bool jacobian_cached_;                  /*  Whether jacobian has been cached */
   int num_constraints_;                   /*  Number of original constraints (for bounds
                                               checking; equal to problem_info->m) */
   int nnz_jacobian_;                       /*  Number of non-zero jacobian entries */

   /**
    * @brief Constructor for FDEvalCache
    * @param num_constraints Number of original constraints (problem_info->m) to allocate
    * space for
    * @param nnz_jacobian Number of non-zero jacobian entries to allocate space for
    */
   FDEvalCache(int num_constraints, int nnz_jacobian)
       : jacobian_cached_(false), num_constraints_(num_constraints), nnz_jacobian_(nnz_jacobian) {
      constraint_values_.resize(num_constraints, 0.0);
      jacobian_values_.resize(nnz_jacobian, 0.0);
   }

   /**
    * @brief Cache jacobian values
    * @param jacobian_values Vector containing jacobian values
    */
   void cacheJacobian(const std::vector<double>& jacobian_values) {
      if (jacobian_values.size() == static_cast<size_t>(nnz_jacobian_)) {
         jacobian_values_ = jacobian_values;
         jacobian_cached_ = true;
      }
   }

   /**
    * @brief Get cached jacobian value for a given index. Trusts that the whole jacobian was
    * populated correctly by cacheJacobian() and that jacobian_idx is in range - callers
    * should check isJacobianCached() first. Bounds are only checked via assert(), which
    * compiles away entirely in release (NDEBUG) builds, so this is just an array read.
    * @param jacobian_idx The jacobian index to get the value for
    * @return the cached value
    */
   double getCachedJacobianValue(int jacobian_idx) const {
      assert(jacobian_idx >= 0 && jacobian_idx < nnz_jacobian_);
      return jacobian_values_[jacobian_idx];
   }

   /**
    * @brief Check if jacobian has been cached
    * @return true if jacobian is cached, false otherwise
    */
   bool isJacobianCached() const {
      return jacobian_cached_;
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
   FDEvalCache* fdeval_cache_;             /*  Cache for FDEvalIni optimization */
   ConoptStatusSolution* status_solution_; /*  Cache for status and solution data */
   Ipopt::OptionsList* options_list_;      /*  Pointer to OptionsList to retrieve optfile name */
   Ipopt::IpoptData* ip_data_;             /*  IpoptData instance for callbacks */
   std::pair<bool, bool> solution_stored_; /*  Whether the stored solution is related to function or
                                              derivative evaluations */
   double* solution_;                      /*  Last point FDEvalIni was called with (length n) */
} IpoptConoptContext;

/**
 * @brief Cleanup function for IpoptConoptContext.
 * Frees any allocated memory in the context.
 */
void CleanupIpoptConoptContext(IpoptConoptContext* context);

/**
 * @brief Get cached constraint value for a given ORIGINAL constraint index (not a CONOPT
 * split-row index - see FDEvalCache::constraint_values_).
 *
 * Trusts that context/row_idx are valid (checked only via assert(), which compiles away in
 * release/NDEBUG builds) rather than branching on them at runtime, so this is just an array
 * read. Defined here (not in the .cpp) and marked inline so the body is visible at every
 * call site for the compiler to actually fold away, rather than merely hinting at it.
 * @param context The context containing the cache
 * @param row_idx The original constraint index to get the value for
 * @return the cached value
 */
inline double GetCachedConstraintValue(IpoptConoptContext* context, int row_idx) {
   assert(context && context->fdeval_cache_);
   FDEvalCache* cache = context->fdeval_cache_;
   assert(row_idx >= 0 && row_idx < cache->num_constraints_);
   return cache->constraint_values_[row_idx];
}

/**
 * @brief Get cached jacobian value for a given index.
 *
 * Trusts that context/jacobian_idx are valid (checked only via assert(), which compiles
 * away in release/NDEBUG builds) rather than branching on them at runtime, so this is just
 * an array read. See GetCachedConstraintValue above for why this is defined inline here.
 * @param context The context containing the cache
 * @param jacobian_idx The jacobian index to get the value for
 * @return the cached value
 */
inline double GetCachedJacobianValue(IpoptConoptContext* context, int jacobian_idx) {
   assert(context && context->fdeval_cache_);
   return context->fdeval_cache_->getCachedJacobianValue(jacobian_idx);
}

/**
 * @brief Check if jacobian has been cached.
 * @param context The context containing the cache
 * @return true if jacobian is cached, false otherwise
 */
bool IsJacobianCached(IpoptConoptContext* context);

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

int COI_CALLCONV Conopt_ReadMatrix(double LOWER[], double CURR[], double UPPER[], int VSTA[],
      int TYPE[], double RHS[], int ESTA[], int COLSTA[], int ROWNO[], double VALUE[], int NLFLAG[],
      int NUMVAR, int NUMCON, int NUMNZ, void* USRMEM);

int COI_CALLCONV Conopt_FDEval(const double X[], double* G, double JAC[], int ROWNO,
      const int JACNUM[], int MODE, int IGNERR, int* ERRCNT, int NUMVAR, int NUMJAC, int THREAD,
      void* USRMEM);

int COI_CALLCONV Conopt_FDEvalIni(const double X[], const int ROWLIST[], int MODE, int LISTSIZE,
      int NUMTHREAD, int IGNERR, int* ERRCNT, int NUMVAR, void* USRMEM);

int COI_CALLCONV Conopt_FDEvalEnd(int IGNERR, int* ERRCNT, void* USRMEM);

int COI_CALLCONV Conopt_Status(int MODSTA, int SOLSTA, int ITER, double OBJVAL, void* USRMEM);

int COI_CALLCONV Conopt_Solution(const double XVAL[], const double XMAR[], const int XBAS[],
      const int XSTA[], const double YVAL[], const double YMAR[], const int YBAS[],
      const int YSTA[], int NUMVAR, int NUMCON, void* USRMEM);

int COI_CALLCONV Conopt_Message(int SMSG, int DMSG, int NMSG, char* MSGV[], void* USRMEM);

int COI_CALLCONV Conopt_ErrMsg(int ROWNO, int COLNO, int POSNO, const char* MSG, void* USRMEM);

/*  Hessian of the Lagrangian structure and values */
int COI_CALLCONV Conopt_2DLagrStr(
      int HSRW[], int HSCL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM);

int COI_CALLCONV Conopt_2DLagrVal(const double X[], const double U[], const int HSRW[],
      const int HSCL[], double HSVL[], int* NODRV, int NUMVAR, int NUMCON, int NHESS, void* USRMEM);

/*  Progress callback for intermediate iteration reporting */
int COI_CALLCONV Conopt_Progress(
      int LEN_INT, const int INT[], int LEN_RL, const double RL[], const double X[], void* USRMEM);

/*  Option callback for setting CONOPT options */
int COI_CALLCONV Conopt_Option(
      int NCALL, double* RVAL, int* IVAL, int* LVAL, char* NAME, void* USRMEM);

/*  ... and so on for all other callbacks (TriOrd, SDDir, etc.) ... */

#ifdef __cplusplus
}
#endif

#endif /*  IPOPT_TO_CONOPT_CALLBACKS_HPP */
