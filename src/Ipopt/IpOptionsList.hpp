#pragma once

/*
 * This is your drop-in replacement for Ipopt/IpOptionsList.hpp
 */

/* 1. INCLUDE ORIGINAL IPOPT HEADERS (for utilities) */
#include "IpReferenced.hpp"
#include "IpTypes.hpp"
#include "IpJournalist.hpp" /* For logging errors */

/* 2. INCLUDE THE CONOPT C-API (for COIDEF_... functions) */
#include "conopt.h"

/* 3. INCLUDE STANDARD C++ HEADERS */
#include <string>
#include <map>
#include <vector>
#include <algorithm> /* For to-lower */
#include <cctype>    /* For to-lower */

namespace Ipopt {

/**
 * @brief This is your shim class for OptionsList.
 * It stores option values set by the user and provides a method
 * to apply them to the CONOPT C-API.
 */
class OptionsList : public ReferencedObject {

 private:
   /* Store the options as the user sets them */
   std::map<std::string, Index> integer_options_;
   std::map<std::string, Number> numeric_options_;
   std::map<std::string, std::string> string_options_;

   /* Store CONOPT options that need to be set via Option callback */
   struct ConoptOption {
      std::string name; /* CONOPT CR-cell name (8 chars, padded with blanks) */
      enum Type { INTEGER, REAL, LOGICAL } type;
      union {
         int ival;
         double rval;
         int lval; /* 0 = false, non-zero = true */
      } value;
   };
   std::vector<ConoptOption> conopt_options_;

   /* Helper to make option matching case-insensitive */
   std::string to_lower(const std::string& s) const {
      std::string out = s;
      std::transform(
            out.begin(), out.end(), out.begin(), [](unsigned char c) { return std::tolower(c); });
      return out;
   }

   /* Helper to pad string to 8 characters with blanks */
   std::string pad_to_8(const std::string& s) const {
      std::string padded = s;
      padded.resize(8, ' ');
      return padded;
   }

   /* Helper to combine prefix and tag for option lookup */
   std::string make_key(const std::string& tag, const std::string& prefix) const {
      if (prefix.empty()) {
         return to_lower(tag);
      }
      return to_lower(prefix + tag);
   }

 public:
   /* --- Public API (Replicating Ipopt) --- */

   OptionsList() {}
   virtual ~OptionsList() {}

   bool SetIntegerValue(
         const std::string& tag, Index value, bool allow_clobber = true, bool dont_print = false) {
      std::string key = to_lower(tag);
      if (!allow_clobber && integer_options_.find(key) != integer_options_.end()) {
         return false;
      }
      integer_options_[key] = value;
      return true;
   }

   bool SetNumericValue(
         const std::string& tag, Number value, bool allow_clobber = true, bool dont_print = false) {
      std::string key = to_lower(tag);
      if (!allow_clobber && numeric_options_.find(key) != numeric_options_.end()) {
         return false;
      }
      numeric_options_[key] = value;
      return true;
   }

   bool SetStringValue(const std::string& tag, const std::string& value, bool allow_clobber = true,
         bool dont_print = false) {
      std::string key = to_lower(tag);
      if (!allow_clobber && string_options_.find(key) != string_options_.end()) {
         return false;
      }
      string_options_[key] = to_lower(value);
      return true;
   }

   bool SetStringValueIfUnset(const std::string& tag, const std::string& value,
         bool allow_clobber = true, bool dont_print = false) {
      std::string key = to_lower(tag);
      if (string_options_.find(key) == string_options_.end()) {
         string_options_[key] = to_lower(value);
      }
      return true;
   }

   bool SetNumericValueIfUnset(
         const std::string& tag, Number value, bool allow_clobber = true, bool dont_print = false) {
      std::string key = to_lower(tag);
      if (numeric_options_.find(key) == numeric_options_.end()) {
         numeric_options_[key] = value;
      }
      return true;
   }

   bool SetIntegerValueIfUnset(
         const std::string& tag, Index value, bool allow_clobber = true, bool dont_print = false) {
      std::string key = to_lower(tag);
      if (integer_options_.find(key) == integer_options_.end()) {
         integer_options_[key] = value;
      }
      return true;
   }

   bool GetIntegerValue(const std::string& tag, Index& value, const std::string& prefix) const {
      std::string key = make_key(tag, prefix);
      auto it = integer_options_.find(key);
      if (it != integer_options_.end()) {
         value = it->second;
         return true;
      }
      return false;
   }

   bool GetNumericValue(const std::string& tag, Number& value, const std::string& prefix) const {
      std::string key = make_key(tag, prefix);
      auto it = numeric_options_.find(key);
      if (it != numeric_options_.end()) {
         value = it->second;
         return true;
      }
      return false;
   }

   bool GetStringValue(
         const std::string& tag, std::string& value, const std::string& prefix) const {
      std::string key = make_key(tag, prefix);
      auto it = string_options_.find(key);
      if (it != string_options_.end()) {
         value = it->second;
         return true;
      }
      return false;
   }

   bool GetEnumValue(const std::string& tag, Index& value, const std::string& prefix) const {
      std::string key = make_key(tag, prefix);
      auto it = integer_options_.find(key);
      if (it != integer_options_.end()) {
         value = it->second;
         return true;
      }
      return false;
   }

   bool GetBoolValue(const std::string& tag, bool& value, const std::string& prefix) const {
      std::string key = make_key(tag, prefix);
      auto it = integer_options_.find(key);
      if (it != integer_options_.end()) {
         value = (it->second != 0);
         return true;
      }
      return false;
   }

   /*
    * @brief Resets the value of a parameter to its default.
    *
    * In CONOPT, the default values are not known, and can not be
    * retrieved. So this method just returns false to indicate to the
    * user that the option is not unset.
    *
    * TODO: consider whether this should just return true, and ignore
    * the action
    */
   bool UnsetValue(const std::string& tag) {
      return false;
   }

   /**
    * @brief Get a CONOPT option by index (for Option callback)
    * @param index The index of the option (0-based, matches NCALL)
    * @param name Output: CONOPT CR-cell name (8 chars, padded)
    * @param rval Output: Real value (if type is REAL, otherwise unchanged)
    * @param ival Output: Integer value (if type is INTEGER, otherwise unchanged)
    * @param lval Output: Logical value (if type is LOGICAL, otherwise unchanged)
    * @param type Output: The type of the option (INTEGER, REAL, or LOGICAL)
    * @return true if option exists, false if index is out of range
    */
   bool GetConoptOption(int index, std::string& name, double* rval, int* ival, int* lval,
         ConoptOption::Type* type) const {
      if (index < 0 || index >= static_cast<int>(conopt_options_.size())) {
         return false;
      }
      const ConoptOption& opt = conopt_options_[index];
      name = opt.name;
      if (type)
         *type = opt.type;
      switch (opt.type) {
      case ConoptOption::INTEGER:
         if (ival)
            *ival = opt.value.ival;
         break;
      case ConoptOption::REAL:
         if (rval)
            *rval = opt.value.rval;
         break;
      case ConoptOption::LOGICAL:
         if (lval)
            *lval = opt.value.lval;
         break;
      }
      return true;
   }

   /**
    * @brief Get the number of CONOPT options stored
    * @return Number of options that will be provided via Option callback
    */
   int GetNumConoptOptions() const {
      return static_cast<int>(conopt_options_.size());
   }

   /* --- Internal Shim Method --- */

   /**
    * @brief Translates and applies all stored options to the CONOPT handle.
    * This is called by the IpoptApplication shim before COI_Solve().
    * @param cntvect CONOPT handle
    * @param jnlst Journalist for logging
    */
   bool ApplyToConopt(coiHandle_t cntvect, Journalist* jnlst) {
      if (!cntvect)
         return false;

      /* Clear any previously stored options */
      conopt_options_.clear();

      /* --- Apply Integer Options --- */
      for (const auto& pair : integer_options_) {
         const std::string& name = pair.first;
         Index value = pair.second;

         if (name == "max_iter") {
            COIDEF_ItLim(cntvect, value);
         }
         else if (name == "threads") {
            COIDEF_ThreadS(cntvect, value);
            COIDEF_ThreadF(cntvect, value);
            COIDEF_Thread2D(cntvect, value);
         }
         else if (name == "print_level") {
            /* Map Ipopt's 0-12 scale to CONOPT's 0-4 (or 0/1) scale */
            int conopt_level = (value == 0) ? 0 : 1; /* Simple 0=off, >0=on */
            COIDEF_StdOut(cntvect, conopt_level);
         }
         else if (name == "conopt_errlim") { /* Custom option for errlim */
            COIDEF_ErrLim(cntvect, value);
         }
         else if (name == "conopt_debug2d") { /* Custom option for debug2d */
            COIDEF_Debug2D(cntvect, value);
         }
         else {
            /* Store options that need to be set via Option callback */
            if (name == "print_frequency_iter") {
               ConoptOption opt;
               opt.name = pad_to_8("logfreq");
               opt.type = ConoptOption::INTEGER;
               opt.value.ival = static_cast<int>(value);
               conopt_options_.push_back(opt);
            }
         }
      }

      /* --- Apply Numeric Options --- */
      for (const auto& pair : numeric_options_) {
         const std::string& name = pair.first;
         Number value = pair.second;

         if (name == "max_cpu_time") {
            COIDEF_ResLim(cntvect, value);
         }
         else if (name == "conopt_maxheap") {  // Custom option
            COIDEF_MaxHeap(cntvect, value);
         }
         else if (name == "tol") {
            ConoptOption opt;
            opt.name = pad_to_8("RTREDG");
            opt.type = ConoptOption::REAL;
            opt.value.rval = value;
            conopt_options_.push_back(opt);
         }
         else if (name == "constr_viol_tol") {
            ConoptOption opt;
            opt.name = pad_to_8("RTNWMA");
            opt.type = ConoptOption::REAL;
            opt.value.rval = value;
            conopt_options_.push_back(opt);
         }
      }

      /* --- Apply String Options --- */
      for (const auto& pair : string_options_) {
         const std::string& name = pair.first;
         const std::string& value = pair.second;

         if (name == "derivative_test") {
            int conopt_val = 0;
            if (value == "first-order" || value == "first-point")
               conopt_val = -1;
            if (value == "only-some-derivatives")
               conopt_val = 1; /* Or some other value */
            COIDEF_DebugFV(cntvect, conopt_val);
         }
      }

      /* --- Always set internal compatibility options --- */
      COIDEF_FVincLin(cntvect, 1);
      COIDEF_FVforAll(cntvect, 1);

      if (jnlst && !conopt_options_.empty()) {
         jnlst->Printf(Ipopt::J_DETAILED, Ipopt::J_MAIN,
               "CONOPT Shim: Stored %d options to be set via Option callback.\n",
               static_cast<int>(conopt_options_.size()));
      }

      return true;
   }

}; /* class OptionsList */

} /* namespace Ipopt */
