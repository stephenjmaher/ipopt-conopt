# Ipopt-to-CONOPT Bridge

A bridge library that allows optimization problems defined using the Ipopt `TNLP` interface to be solved using CONOPT instead of Ipopt. This provides a seamless way to use CONOPT as a drop-in replacement for Ipopt in existing codebases.

The bridge is not fully featured and there may be features in Ipopt used in an application that are not implemented. If
this is the case, please create an issue requesting the feature.

## Overview

This project implements a bridge between two optimization solvers:

- **Ipopt**: An open-source nonlinear optimization solver -- the source code is available at https://github.com/coin-or/Ipopt
- **CONOPT**: A commercial nonlinear optimization solver from GAMS -- available at https://conopt.gams.com/

While CONOPT is a commercial solver, it is free to use for academic purposes. The licensing details can be found at https://conopt.gams.com/licensing/.

The bridge translates Ipopt's C++ `TNLP` (Tagged Nonlinear Programming) interface to CONOPT's C API, allowing you to solve your optimization problems with CONOPT while maintaining compatibility with Ipopt-based code.

## Key Features

- **Drop-in Replacement**: Provides a shim `IpoptApplication` class that mimics Ipopt's interface
- **Automatic Problem Transformation**: Handles constraint splitting, objective function conversion, and Hessian reordering
- **Callback Translation**: Implements all required CONOPT C API callbacks (ReadMatrix, FDEval, Solution, Status, etc.)
- **Status Code Mapping**: Converts CONOPT status codes to Ipopt-compatible return codes

## Project Structure

### Source Code (`src/`)

- **`IpoptToConoptCallbacks.hpp/cpp`**: Core bridge implementation
  - C-style trampoline functions that implement CONOPT's callback interface
  - Converts between Ipopt `TNLP` methods and CONOPT callbacks
  - Handles caching for performance optimization

- **`IpoptProblemInfo.hpp`**: Problem information structure
  - Stores problem dimensions, bounds, Jacobian/Hessian structures
  - Handles constraint splitting (range constraints -> two inequalities)
  - Manages mappings between original and split constraint formulations

- **`IpoptTypes.hpp`**: Basic type definitions
  - `Index` and `Number` type aliases
  - Utility functions for infinity handling

- **`Ipopt/`**: Ipopt interface headers
  - `IpIpoptApplication.hpp`: Application shim class
  - `IpSolveStatistics.hpp`: Statistics tracking
  - `IpOptionsList.hpp`: Options management

## Building Your Project with the Bridge

### Prerequisites

To build a project using this Ipopt-to-CONOPT bridge, you will need:

1.  **A C++11 (or newer) Compiler** (e.g., g++, Clang, MSVC).

2.  **CONOPT Solver:** The CONOPT C-API headers (e.g., `conopt.h`) and the compiled binary library (e.g., `libconopt.so` or `conopt.lib`).

3.  **Original Ipopt Installation:** The full Ipopt header files (e.g., `IpTypes.hpp`, `IpSmartPtr.hpp`, etc.) and the compiled binary library (e.g., `libipopt.so` or `libipopt.a`).

4.  **This Bridge's Source Code:** All the header (`.hpp`) files and the single implementation file (`IpoptToConoptCallbacks.cpp`) from this repository.

### Integration

Integrating the bridge into your project involves three main steps:

1.  **Compile the Bridge:** Add the bridge's implementation file, `IpoptToConoptCallbacks.cpp`, to the list of source files your build system compiles.

2.  **Set Include Path Order:** Configure your build system to find the **bridge's headers *before*** the original Ipopt headers.

3.  **Link Libraries:** Link your final executable against the compiled bridge object, the CONOPT library, and the Ipopt library.

**Note:** You **must** link against the original Ipopt library (`libipopt.so` or `libipopt.a`). This bridge re-uses Ipopt's underlying object model (like `ReferencedObject` and `SmartPtr`), which requires linking to the original Ipopt library to resolve base class destructors and other utility functions.

### Modifications to existing source code

There are some minor changes required to your existing Ipopt interface source code. This is mainly related to removing
heading includes that are not necessary when using CONOPT. 

- Remove the include of header files that conflict with the Ipopt-to-CONOPT bridge. Some examples are: `IpTNLPAdapter.hpp`
  and `IpOrigIpoptNLP.hpp`. These examples have come up in my testing, so there are likely more.

- Remove references to `IpIpoptCalculatedQuantities` in callbacks. The `IpIpoptCalculatedQuantities` object is not
  populated by the Ipopt-to-CONOPT bridge. This is mainly because many of the calculated quantities are related to the
  Ipopt algorithm. As such, these are not relevant for CONOPT. It is necessary to remove references to
  `IpIpoptCalculatedQuantities` in the `finalize_solution` and `intermediate_callback` callbacks. In these callbacks
  `IpIpoptCalculatedQuantities` is passed as a null pointer.

- Provide the nonlinear Jacobian entries using the `get_nonlinear_terms` method. This is optional, and only needed if
  CONOPT should use second order information. While the nonlinear terms in the Jacobian do not need to be specified for
  the function and derivative evaluations, they are needed when specifying the structure of the Hessian of the
  Lagrangian. When the user provides the structure of the Hessian in `eval_h`, without specifying the nonlinear terms,
  there will be a conflict in CONOPT w.r.t to the linear terms in the Jacobian. In the bridge, the default is to assume
  that all terms are nonlinear in the Jacobian. This is needed because CONOPT doesn't evaluate linear terms, but Ipopt
  does. However, assuming all terms are nonlinear means that the structure of the Hessian doesn't match with the
  structure of the Jacobian. As such, a method has been added to the `TNLP` class, `get_nonlinear_terms` that can be
  implemented to inform CONOPT of the terms in the Jacobian that are nonlinear. Then there is no conflict between the
  Jacobian and the Hessian, and CONOPT can use second order information.

### CMake Example

Using CMake is the recommended way to build your project. The following `CMakeLists.txt` demonstrates the correct setup.

```cmake
cmake_minimum_required(VERSION 3.10)
project(MyIpoptProject)

# ------------------------------------------------------------------
# 1. DEFINE YOUR PROJECT VARIABLES
# (Set these paths based on your system)
# ------------------------------------------------------------------

# Your application's source files
set(MY_APP_SOURCES
    src/MyProblem.cpp
    src/main.cpp
)

# Path to the root of this bridge repository
# This directory should contain 'IpoptToConoptCallbacks.cpp' and the 'Ipopt/' header folder
set(IP2CO_DIR /path/to/ipopt-conopt-bridge)

# Path to your CONOPT C-API installation
set(CONOPT_DIR /path/to/conopt)

# Path to your original Ipopt installation
set(IPOPT_DIR /path/to/ipopt_install)

# ------------------------------------------------------------------
# 2. CONFIGURE THE EXECUTABLE
# ------------------------------------------------------------------

# Add your executable
add_executable(my_app
    ${MY_APP_SOURCES}
    
    # Add the bridge's implementation file to be compiled
    ${IP2CO_DIR}/IpoptToConoptCallbacks.cpp
)

# ------------------------------------------------------------------
# 3. SET INCLUDE PATH ORDER (CRITICAL)
# ------------------------------------------------------------------
target_include_directories(my_app PRIVATE
    
    # 1. The bridge's headers MUST come FIRST
    ${IP2CO_DIR}
    
    # 2. The original Ipopt headers come SECOND
    ${IPOPT_DIR}/include/coin-or
    
    # 3. The CONOPT C-API headers
    ${CONOPT_DIR}/include
)

# ------------------------------------------------------------------
# 4. LINK LIBRARIES
# ------------------------------------------------------------------

# Add the library search paths
target_link_directories(my_app PRIVATE
    ${CONOPT_DIR}/lib
    ${IPOPT_DIR}/lib
)

# Link your app against CONOPT and IPOPT
target_link_libraries(my_app PRIVATE
    conopt # Assumes libconopt.so or conopt.lib
    ipopt  # Assumes libipopt.so or libipopt.a
    # Add other necessary libraries (e.g., m, dl)
    m
    dl
)
```

### Makefile Example

If you are not using CMake, you can use a traditional `Makefile`. The key is to add the bridge's `.cpp` file to your sources and set the include path order using `-I` flags.

```makefile
# ------------------------------------------------------------------
# 1. DEFINE YOUR PROJECT VARIABLES
# (Set these paths based on your system)
# ------------------------------------------------------------------

# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -Wall -g # -g for debugging

# Your application's executable name
TARGET = my_app

# Path to the root of this bridge repository
IP2CO_DIR = /path/to/ipopt-conopt-bridge

# Path to your CONOPT C-API installation
CONOPT_DIR = /path/to/conopt

# Path to your original Ipopt installation
IPOPT_DIR = /path/to/ipopt_install

# ------------------------------------------------------------------
# 2. DEFINE SOURCES AND OBJECTS
# ------------------------------------------------------------------

# Add the bridge's implementation file to the list of sources
SRCS = \
    src/MyProblem.cpp \
    src/main.cpp \
    ${IP2CO_DIR}/IpoptToConoptCallbacks.cpp

# Create a list of object files
OBJS = $(SRCS:.cpp=.o)

# ------------------------------------------------------------------
# 3. SET INCLUDE PATH ORDER (CRITICAL)
# ------------------------------------------------------------------
# The bridge's headers MUST come FIRST.
INCLUDE_PATHS = \
    -I${IP2CO_DIR} \
    -I${IPOPT_DIR}/include/coin-or \
    -I${CONOPT_DIR}/include

# Add include paths to compiler flags
CXXFLAGS += ${INCLUDE_PATHS}

# ------------------------------------------------------------------
# 4. SET LINKER FLAGS
# ------------------------------------------------------------------
LDFLAGS = \
    -L${CONOPT_DIR}/lib \
    -L${IPOPT_DIR}/lib \
    -Wl,-rpath,${CONOPT_DIR}/lib \
    -Wl,-rpath,${IPOPT_DIR}/lib \
    -lconopt -lipopt -lm -ldl

# ------------------------------------------------------------------
# 5. BUILD RULES
# ------------------------------------------------------------------

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Generic rule for compiling .cpp files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
```


## Examples

Two example programs demonstrating how to use the bridge. These examples are taken directly from the Ipopt repository.
The Makefiles have been modified to demonstrate how to build the examples using the bridge.

### 1. **Cpp_example** - Simple Tutorial Example

A minimal example solving:

$$
\begin{align}
\min \quad & -(x_2-2)^2 \\
\text{s.t.} \quad & x_1^2 + x_2 = 1 \\
& -1 \leq x_1 \leq 1
\end{align}
$$

**Files:**
- `cpp_example.cpp`: Main program
- `MyNLP.hpp/cpp`: Simple TNLP implementation

### 2. **hs071_cpp** - Hock-Schittkowski Problem 71

A standard test problem:

$$
\begin{align}
\min \quad & x_1 x_4 (x_1 + x_2 + x_3) + x_3 \\
\text{s.t.} \quad & x_1 x_2 x_3 x_4 \geq 25 \\
& x_1^2 + x_2^2 + x_3^2 + x_4^2 = 40 \\
& 1 \leq x_1, x_2, x_3, x_4 \leq 5
\end{align}
$$

**Files:**
- `hs071_main.cpp`: Main program with options configuration
- `hs071_nlp.hpp/cpp`: HS071 TNLP implementation

### Building the Examples

To build an example, navigate to the example directory and run `make` with all required path variables:

```bash
cd example/hs071_cpp
make CONOPT_DIR=/path/to/conopt \
     IP2CO_DIR=/path/to/ipopt-conopt/src \
     IPOPT_DIR=/path/to/ipopt \
     OPT=opt
```

**Build Options:**
- `OPT=opt`: Optimized build (default, uses `-O2 -DNDEBUG`)
- `OPT=dbg`: Debug build (uses `-g -O0 -DDEBUG`)

**Example with absolute paths:**
```bash
cd example/hs071_cpp
make CONOPT_DIR=/opt/conopt \
     IP2CO_DIR=/home/user/ipopt-conopt/src \
     IPOPT_DIR=/usr/local \
     OPT=opt
```

The Makefiles use the following include order (important for header resolution):
1. CONOPT include directory
2. Bridge source directory (`src/`)
3. Bridge Ipopt interface directory (`src/Ipopt/`)
4. IPOPT include directory

This order ensures that the bridge's shim headers are found before the original IPOPT headers.


## How It Works

### Problem Transformation

The bridge performs several transformations to convert Ipopt problems to CONOPT format:

1. **Constraint Splitting**: Range constraints (with both lower and upper bounds) are split into two separate inequality constraints
2. **Objective Handling**: The objective function is treated as a special "free" constraint row
3. **Hessian Reordering**: CONOPT requires Hessian entries sorted by column then row
4. **Index Conversion**: Handles both C-style (0-based) and Fortran-style (1-based) indexing

### Callback Implementation: Mapping Ipopt to CONOPT

This bridge works by implementing the CONOPT C-API callback functions (the "trampolines"). These trampolines are registered with CONOPT, and when called, they cast the `void* USRMEM` cookie back into an `IpoptConoptContext*` object. This context allows the trampolines to access the user's `Ipopt::TNLP` and `Ipopt::Journalist` objects, bridging the gap between the two APIs.

The core translation logic is as follows:

* **`Conopt_ReadMatrix`**
  This callback is used by CONOPT to load the entire problem definition at once. It maps to several `Ipopt::TNLP` methods (`get_bounds_info`, `get_starting_point`, and `eval_jac_g` for structure). To handle this, the `IpoptApplication` bridge calls these `TNLP` methods *before* the solve to pre-cache all problem data. This trampoline then simply copies that cached data into the arrays provided by CONOPT.

* **`Conopt_FDEvalIni`**
  These are CONOPT callbacks that are called before and after the function evaluation rounds. Since CONOPT executes a
  GRG algorithm, compared to an interior point algorithm, the function evaluation process is a little different. The
  main difference is that CONOPT can request the evaluation of a subset of rows, as opposed to all rows. To avoid
  calling the evaluation methods too often, `Conopt_FDEvalIni` executes `tnlp->eval_f` and `tnlp->eval_g`, for the
  objective and constraint functions, and `tnlp->eval_grad_f` and `tnlp->eval_jac_g`, for the objective and constraint
  derivatives. The results of these evaluations are cached and then used in `Conopt_FDEval`.

* **`Conopt_FDEval`**
  Extracts the evaluation results from the cache and return this to CONOPT. Since CONOPT may only request evaluations
  from a subset of rows, caching the results in `Conopt_FDEvalIni` is critically important.

* **`Conopt_2DLagrStr` / `Conopt_2DLagrVal`**
  These callbacks map directly to Ipopt's two-stage Hessian evaluation (`tnlp->eval_h`).

  * `Conopt_2DLagrStr` is called first to get the sparsity structure (iRow, jCol) of the Hessian.

  * `Conopt_2DLagrVal` is called later to get the numerical values for that structure, passing in the current `obj_factor` and `lambda` (multipliers).

* **`Conopt_Progress`**
  This provides intermediate iteration updates. `Conopt_Progress` is the primary iteration log, and its trampoline calls `tnlp->intermediate_callback`, allowing the user to stop the solve.

* **`Conopt_Status`**
  This reports the final solver status. The trampoline uses the `MODSTA`/`SOLSTA` codes from this callback to populate the `SolveStatistics` bridge.

* **`ConGpt_Solution`**
  This callback delivers the final solution. The trampoline receives CONOPT's solution arrays (`XVAL`, `YMAR`, etc.) and uses them to call `tnlp->finalize_solution` and to populate the `SolveStatistics` bridge with the final primal and dual values.

* **`Conopt_Message` / `Conopt_ErrMsg`**
  These are the logging callbacks. The trampoline intercepts all messages and routes them to the `Ipopt::Journalist` object (from the `USRMEM` context), which respects the user's `print_level` settings.

* **`Conopt_Option`**
  This is an optional callback for handling solver options. If used, the trampoline looks up the option requested by CONOPT (e.g., "tol_opt") in the `OptionsList` bridge and returns the corresponding value (e.g., from Ipopt's "tol").

### Status Code Mapping

The bridge translates CONOPT's dual status code system (MODSTA and SOLSTA) to Ipopt's `ApplicationReturnStatus` enum. CONOPT uses:
- **MODSTA (Model Status)**: Describes the state of the model/solution (optimal, infeasible, unbounded, etc.)
- **SOLSTA (Solver Status)**: Describes how the solver terminated (normal completion, iteration limit, error, etc.)

The mapping logic prioritizes SOLSTA, then considers MODSTA for detailed interpretation. Below is the complete mapping table:

#### Normal Completion (SOLSTA = 1)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| 1 | Optimal | `Solve_Succeeded` |
| 2 | Locally Optimal | `Solve_Succeeded` |
| 7 | Feasible Solution | `Solve_Succeeded` |
| 4 | Locally Infeasible | `Infeasible_Problem_Detected` |
| 5 | Infeasible | `Infeasible_Problem_Detected` |
| Other | Other completion status | `Solved_To_Acceptable_Level` |

#### Iteration Interrupt (SOLSTA = 2)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | Iteration limit reached | `Maximum_Iterations_Exceeded` |

#### Resource/CPU Time Limit (SOLSTA = 3)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | CPU time limit exceeded | `Maximum_CpuTime_Exceeded` |

#### Terminated by Solver (SOLSTA = 4)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| 3 | Unbounded (termination) | `Diverging_Iterates` |
| 4 | Locally Infeasible (termination) | `Restoration_Failed` |
| 5 | Infeasible (termination) | `Infeasible_Problem_Detected` |
| 6 | Intermediate Infeasible (termination) | `Search_Direction_Becomes_Too_Small` |
| 13 | Error No Solution (termination) | `Error_In_Step_Computation` |
| Other | Other termination reason | `Internal_Error` |

#### Evaluation Error Limit (SOLSTA = 5)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | Evaluation error limit exceeded | `Invalid_Number_Detected` |

#### Error (SOLSTA = 6)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | Solver error | `Internal_Error` |

#### User Interrupt (SOLSTA = 8)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | User requested stop | `User_Requested_Stop` |

#### Out of Memory (SOLSTA = 9)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| Any | Out of memory | `Insufficient_Memory` |

#### System Error / Invalid Setup (SOLSTA = 10)

| MODSTA | CONOPT Meaning | Ipopt `ApplicationReturnStatus` |
|--------|----------------|----------------------------------|
| 13 | Invalid problem definition | `Invalid_Problem_Definition` |
| Other | Other system error | `Internal_Error` |

#### Special Case: User Stop (MODSTA = 10)

If MODSTA = 10 (regardless of SOLSTA), the bridge returns `User_Requested_Stop` to indicate explicit user interruption.

#### ApplicationReturnStatus to SolverReturn Conversion

The bridge also converts `ApplicationReturnStatus` to Ipopt's `SolverReturn` enum (used in `finalize_solution`):

| ApplicationReturnStatus | SolverReturn |
|------------------------|--------------|
| `Solve_Succeeded` | `SUCCESS` |
| `Solved_To_Acceptable_Level` | `STOP_AT_ACCEPTABLE_POINT` |
| `Infeasible_Problem_Detected` | `LOCAL_INFEASIBILITY` |
| `Search_Direction_Becomes_Too_Small` | `STOP_AT_TINY_STEP` |
| `Diverging_Iterates` | `DIVERGING_ITERATES` |
| `User_Requested_Stop` | `USER_REQUESTED_STOP` |
| `Feasible_Point_Found` | `FEASIBLE_POINT_FOUND` |
| `Maximum_Iterations_Exceeded` | `MAXITER_EXCEEDED` |
| `Maximum_CpuTime_Exceeded` | `CPUTIME_EXCEEDED` |
| `Error_In_Step_Computation` | `ERROR_IN_STEP_COMPUTATION` |
| `Invalid_Number_Detected` | `INVALID_NUMBER_DETECTED` |
| `Internal_Error` | `INTERNAL_ERROR` |

## License

See `LICENSE` file for license information.

## Contributing

This bridge is designed to maintain compatibility with Ipopt's interface while leveraging CONOPT's solver capabilities. When contributing:

- Maintain API compatibility with Ipopt's `TNLP` interface
- Ensure proper status code mappings
- Test with the provided examples
- Document any constraint transformations or limitations
