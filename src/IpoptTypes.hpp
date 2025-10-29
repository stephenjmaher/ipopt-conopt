#ifndef IPOPT_TYPES_HPP
#define IPOPT_TYPES_HPP

#include <limits>
#include <cmath>

namespace Ipopt {

// Basic types used by IPOPT - keeping these as requested
typedef int Index;
typedef double Number;

// Constants
const Number IPOPT_INFINITY = std::numeric_limits<Number>::infinity();
const Number IPOPT_NEG_INFINITY = -std::numeric_limits<Number>::infinity();

// Utility function to check if a number is finite
inline bool IsFiniteNumber(Number val) {
   return std::isfinite(val);
}

// Note: SolverReturn and other enums are now provided by the real IPOPT headers
// This file only provides the basic type definitions and utility functions

}  // namespace Ipopt

#endif  // IPOPT_TYPES_HPP