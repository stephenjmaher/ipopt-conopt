// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter          IBM    2005-08-15
// Modified for CONOPT integration

#ifndef __IPSOLVESTATISTICS_HPP__
#define __IPSOLVESTATISTICS_HPP__

#include "IpReferenced.hpp"
#include "IpSmartPtr.hpp"
#include "IpUtils.hpp"
#include <chrono>
#include <ctime>

namespace Ipopt
{

class IPOPTLIB_EXPORT SolveStatistics : public ReferencedObject
{
public:

   SolveStatistics() : num_iters_(0),
       total_cpu_time_(0.0),
       total_sys_time_(0.0),
       total_wallclock_time_(0.0),
       num_obj_evals_(0),
       num_constr_evals_(0),
       num_obj_grad_evals_(0),
       num_constr_jac_evals_(0),
       num_hess_evals_(0),
       scaled_obj_val_(0.0),
       obj_val_(0.0),
       scaled_dual_inf_(0.0),
       dual_inf_(0.0),
       scaled_constr_viol_(0.0),
       constr_viol_(0.0),
       scaled_bound_viol_(0.0),
       bound_viol_(0.0),
       scaled_compl_(0.0),
       compl_(0.0),
       scaled_kkt_error_(0.0),
       kkt_error_(0.0),
       solve_status_(Solve_Succeeded)
   {
      // Constructor implementation - for CONOPT integration, we populate data manually
   }

   virtual ~SolveStatistics()
   { }

   virtual Index IterationCount() const
   {
      return num_iters_;
   }

   virtual Number TotalCpuTime() const
   {
      return total_cpu_time_;
   }

   virtual Number TotalSysTime() const
   {
      return total_sys_time_;
   }

   virtual Number TotalWallclockTime() const
   {
      return total_wallclock_time_;
   }

   virtual void NumberOfEvaluations(
      Index& num_obj_evals,
      Index& num_constr_evals,
      Index& num_obj_grad_evals,
      Index& num_constr_jac_evals,
      Index& num_hess_evals
   ) const
   {
      num_obj_evals = num_obj_evals_;
      num_constr_evals = num_constr_evals_;
      num_obj_grad_evals = num_obj_grad_evals_;
      num_constr_jac_evals = num_constr_jac_evals_;
      num_hess_evals = num_hess_evals_;
   }

   virtual void Infeasibilities(
      Number& dual_inf,
      Number& constr_viol,
      Number& varbounds_viol,
      Number& complementarity,
      Number& kkt_error
   ) const
   {
      dual_inf = dual_inf_;
      constr_viol = constr_viol_;
      varbounds_viol = bound_viol_;
      complementarity = compl_;
      kkt_error = kkt_error_;
   }

   virtual void ScaledInfeasibilities(
      Number& scaled_dual_inf,
      Number& scaled_constr_viol,
      Number& scaled_varbounds_viol,
      Number& scaled_complementarity,
      Number& scaled_kkt_error
   ) const
   {
      scaled_dual_inf = scaled_dual_inf_;
      scaled_constr_viol = scaled_constr_viol_;
      scaled_varbounds_viol = scaled_bound_viol_;
      scaled_complementarity = scaled_compl_;
      scaled_kkt_error = scaled_kkt_error_;
   }

   virtual Number FinalObjective() const
   {
      return obj_val_;
   }

   virtual Number FinalScaledObjective() const
   {
      return scaled_obj_val_;
   }

   // Additional method to get solve status (not in original but useful for CONOPT)
   virtual SolverReturn SolveStatus() const
   {
      return solve_status_;
   }

   // Timing methods for CONOPT integration
   void StartTiming();
   void StopTiming();

   // Setter methods for populating the statistics (not in original but needed for CONOPT)
   void SetIterationCount(Index iters)
   {
      num_iters_ = iters;
   }


   void SetInfeasibilities(
      Number dual_inf,
      Number constr_viol,
      Number bound_viol,
      Number compl,
      Number kkt_error
   )
   {
      dual_inf_ = dual_inf;
      constr_viol_ = constr_viol;
      bound_viol_ = bound_viol;
      compl_ = compl;
      kkt_error_ = kkt_error;
   }

   void SetScaledInfeasibilities(
      Number scaled_dual_inf,
      Number scaled_constr_viol,
      Number scaled_bound_viol,
      Number scaled_compl,
      Number scaled_kkt_error
   )
   {
      scaled_dual_inf_ = scaled_dual_inf;
      scaled_constr_viol_ = scaled_constr_viol;
      scaled_bound_viol_ = scaled_bound_viol;
      scaled_compl_ = scaled_compl;
      scaled_kkt_error_ = scaled_kkt_error;
   }

   void SetFinalObjective(Number obj_val)
   {
      obj_val_ = obj_val;
   }

   void SetFinalScaledObjective(Number scaled_obj_val)
   {
      scaled_obj_val_ = scaled_obj_val;
   }

   void SetSolveStatus(SolverReturn status)
   {
      solve_status_ = status;
   }

   // Evaluation counting methods for CONOPT integration
   void IncrementObjectiveEvaluations(Index count = 1)
   {
      num_obj_evals_ += count;
   }

   void IncrementConstraintEvaluations(Index count = 1)
   {
      num_constr_evals_ += count;
   }

   void IncrementObjectiveGradientEvaluations(Index count = 1)
   {
      num_obj_grad_evals_ += count;
   }

   void IncrementConstraintJacobianEvaluations(Index count = 1)
   {
      num_constr_jac_evals_ += count;
   }

   void IncrementHessianEvaluations(Index count = 1)
   {
      num_hess_evals_ += count;
   }

   // Timing methods for CONOPT integration
   void StartTiming()
   {
      start_cpu_time_ = GetCurrentCpuTime();
      start_wall_time_ = GetCurrentWallTime();
      start_sys_time_ = GetCurrentSysTime();
   }

   void StopTiming()
   {
      total_cpu_time_ = GetCurrentCpuTime() - start_cpu_time_;
      total_wallclock_time_ = GetCurrentWallTime() - start_wall_time_;
      total_sys_time_ = GetCurrentSysTime() - start_sys_time_;
   }

private:
   // Helper methods for timing
   Number GetCurrentCpuTime() const
   {
      // Use clock() for CPU time
      return static_cast<Number>(clock()) / static_cast<Number>(CLOCKS_PER_SEC);
   }

   Number GetCurrentWallTime() const
   {
      // Use std::chrono for wall clock time
      auto now = std::chrono::high_resolution_clock::now();
      auto duration = now.time_since_epoch();
      return std::chrono::duration_cast<std::chrono::microseconds>(duration).count() / 1000000.0;
   }

   Number GetCurrentSysTime() const
   {
      // For system time, we'll use a simplified approach
      // In a real implementation, you might want to use getrusage() or similar
      return GetCurrentCpuTime(); // Simplified - same as CPU time
   }

private:

   SolveStatistics();

   SolveStatistics(
      const SolveStatistics&
   );

   void operator=(
      const SolveStatistics&
   );

   Index num_iters_;
   /* Total CPU time */
   Number total_cpu_time_;
   /* Total system time */
   Number total_sys_time_;
   /* Total wall clock time */
   Number total_wallclock_time_;
   Index num_obj_evals_;
   Index num_constr_evals_;
   Index num_obj_grad_evals_;
   Index num_constr_jac_evals_;
   Index num_hess_evals_;

   Number scaled_obj_val_;
   Number obj_val_;
   Number scaled_dual_inf_;
   Number dual_inf_;
   Number scaled_constr_viol_;
   Number constr_viol_;
   Number scaled_bound_viol_;
   Number bound_viol_;
   Number scaled_compl_;
   Number compl_;
   Number scaled_kkt_error_;
   Number kkt_error_;

   // Additional field for solve status (not in original but useful for CONOPT)
   SolverReturn solve_status_;

   // Timing fields for CONOPT integration
   Number start_cpu_time_;
   Number start_wall_time_;
   Number start_sys_time_;
};

} // namespace Ipopt

#endif
