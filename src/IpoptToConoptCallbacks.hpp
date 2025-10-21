/**
 * @file IpoptToConoptCallbacks.hpp
 * @brief Declares the C-style trampolines to bridge Ipopt::TNLP to the CONOPT C-API.
 */
#ifndef IPOPT_TO_CONOPT_CALLBACKS_HPP
#define IPOPT_TO_CONOPT_CALLBACKS_HPP

#include "conopt.h" // For COI_CALLCONV and C-API types

// Forward declare Ipopt classes needed by the struct
namespace Ipopt {
   class TNLP;
   class Journalist;
   class SolveStatistics;
   struct IpoptProblemInfo;
}

/**
 * @brief Struct to hold the context needed by the C trampolines.
 * This will be passed as the void* USRMEM cookie.
 */
typedef struct {
   Ipopt::TNLP* tnlp_;
   Ipopt::Journalist* journalist_;
   Ipopt::SolveStatistics* stats_;
   Ipopt::IpoptProblemInfo* problem_info_;
} IpoptConoptContext;

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
