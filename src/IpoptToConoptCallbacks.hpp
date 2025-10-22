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
 * @brief Struct to hold cached constraint values for FDEvalIni optimization.
 * This stores the results of constraint evaluations with constant lookup time.
 */
struct FDEvalCache {
   std::vector<double> constraint_values_;  // Cached constraint values (size = numcons)
   std::vector<bool> constraint_valid_;     // Validity flags for each constraint (size = numcons)
   double objective_value_;                // Cached objective value
   bool objective_valid_;                  // Whether objective value is valid
   int num_constraints_;                   // Number of constraints (for bounds checking)

   /**
    * @brief Constructor for FDEvalCache
    * @param num_constraints Number of constraints to allocate space for
    */
   FDEvalCache(int num_constraints)
      : num_constraints_(num_constraints),
        objective_value_(0.0),
        objective_valid_(false) {
      constraint_values_.resize(num_constraints, 0.0);
      constraint_valid_.resize(num_constraints, false);
   }

   /**
    * @brief Reset all validity flags to false
    */
   void invalidateAll() {
      std::fill(constraint_valid_.begin(), constraint_valid_.end(), false);
      objective_valid_ = false;
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
