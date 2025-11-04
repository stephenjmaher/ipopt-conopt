# Ipopt-to-CONOPT Bridge

A bridge library that allows optimization problems defined using the Ipopt `TNLP` interface to be solved using CONOPT instead of Ipopt. This provides a seamless way to use CONOPT as a drop-in replacement for Ipopt in existing codebases.

## Overview

This project implements a bridge between two optimization solvers:

- **Ipopt**: An open-source nonlinear optimization solver
- **CONOPT**: A commercial nonlinear optimization solver from GAMS

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
  - `IpTNLP.hpp`: TNLP base class interface
  - `IpIpoptApplication.hpp`: Application shim class
  - `IpSolveStatistics.hpp`: Statistics tracking
  - `IpOptionsList.hpp`: Options management

### Examples (`example/`)

Two example programs demonstrating how to use the bridge. These examples are taken directly from the Ipopt repository.
The Makefiles have been modified to demonstrate how to build the examples using the bridge.

#### 1. `Cpp_example/` - Simple Tutorial Example

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

#### 2. `hs071_cpp/` - Hock-Schittkowski Problem 71

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

## How It Works

### Problem Transformation

The bridge performs several transformations to convert Ipopt problems to CONOPT format:

1. **Constraint Splitting**: Range constraints (with both lower and upper bounds) are split into two separate inequality constraints
2. **Objective Handling**: The objective function is treated as a special "free" constraint row
3. **Hessian Reordering**: CONOPT requires Hessian entries sorted by column then row
4. **Index Conversion**: Handles both C-style (0-based) and Fortran-style (1-based) indexing

### Callback Implementation

The bridge implements all required CONOPT callbacks:

- **`Conopt_ReadMatrix`**: Initial problem setup, extracts problem structure from TNLP
- **`Conopt_FDEval`**: Evaluates constraint values and Jacobian
- **`Conopt_FDEvalIni/End`**: Batch evaluation optimization
- **`Conopt_Status`**: Receives CONOPT status updates
- **`Conopt_Solution`**: Receives final solution
- **`Conopt_2DLagrStr/Val`**: Hessian of Lagrangian structure and values
- **`Conopt_Progress`**: Intermediate iteration reporting
- **`Conopt_Option`**: Option handling

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

## Building

### Prerequisites

- **CONOPT**: CONOPT solver library (C API) - the library files and headers are expected in `external/conopt-linux-x86_64/`
- **IPOPT**: IPOPT library - compiled and installed IPOPT with headers
- **C++ Compiler**: C++11 or later (tested with `g++`)

### Setting Up Library and Include Paths

When building projects that use this bridge, you need to explicitly specify paths to both CONOPT and IPOPT libraries and includes. All paths must be provided - there are no defaults.

#### Required Build Variables

The build system requires the following variables to be set:

1. **`CONOPT_DIR`**: Path to CONOPT installation directory
   - Example: `/path/to/conopt` or `/opt/conopt`
   - Must contain `include/` subdirectory with `conopt.h` (CONOPT C API header)
   - Must contain `lib/` subdirectory with `libconopt.so` and related shared libraries
   - Used for both include paths (`CONOPT_DIR/include`) and library paths (`CONOPT_DIR/lib`)

2. **`IP2CO_INC`**: Path to bridge source/include directory
   - Example: `/path/to/ipopt-conopt/src`
   - Contains bridge headers (`IpoptToConoptCallbacks.hpp`, `IpoptProblemInfo.hpp`, etc.)
   - Also contains bridge source file (`IpoptToConoptCallbacks.cpp`) that must be compiled

3. **`IPOPT_DIR`**: Path to IPOPT installation directory
   - Example: `/path/to/ipopt` or `/usr/local`
   - Must contain `include/coin-or/` subdirectory with IPOPT headers
   - Must contain `lib/` subdirectory with `libipopt.so` and related shared libraries
   - Used for both include paths (`IPOPT_DIR/include/coin-or`) and library paths (`IPOPT_DIR/lib`)

Note that when building your Ipopt project with the Ipopt-CONOPT bridge, it is important to specify the include path for
the bridge before the Ipopt include path. Additionally, it is still necessary to link against the Ipopt library, because
some data types and structures are defined by Ipopt, for example `Index` and `Number`.

#### Building the Examples

To build an example, navigate to the example directory and run `make` with all required path variables:

```bash
cd example/hs071_cpp
make CONOPT_DIR=/path/to/conopt \
     IP2CO_INC=/path/to/ipopt-conopt/src \
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
     IP2CO_INC=/home/user/ipopt-conopt/src \
     IPOPT_DIR=/usr/local \
     OPT=opt
```

The Makefiles use the following include order (important for header resolution):
1. CONOPT include directory
2. Bridge source directory (`src/`)
3. Bridge Ipopt interface directory (`src/Ipopt/`)
4. IPOPT include directory

This order ensures that the bridge's shim headers are found before the original IPOPT headers.

#### Linking

The linker requires both CONOPT and IPOPT libraries. The Makefiles link them in this order:
- `-lconopt` (CONOPT C API library)
- `-lipopt` (IPOPT library)
- Standard C++ libraries (`-lstdc++`, `-lm`, `-ldl`)

Runtime library paths are set using `-Wl,--rpath` flags to ensure the correct libraries are found at runtime.

### Building Your Own Project

To integrate this bridge into your own project, you must compile the bridge's implementation file (`IpoptToConoptCallbacks.cpp`) and link it with your code.

#### Step 1: Set Up Include Paths

Add the following include paths to your compiler flags:

```bash
-I$(CONOPT_DIR)/include \
-I$(IP2CO_INC) \
-I$(IP2CO_INC)/Ipopt \
-I$(IPOPT_DIR)/include/coin-or
```

#### Step 2: Compile the Bridge Implementation

**You must compile `IpoptToConoptCallbacks.cpp`** as part of your build process:

```bash
g++ -c $(IP2CO_INC)/IpoptToConoptCallbacks.cpp \
    -I$(CONOPT_DIR)/include \
    -I$(IP2CO_INC) \
    -I$(IP2CO_INC)/Ipopt \
    -I$(IPOPT_DIR)/include/coin-or \
    -o IpoptToConoptCallbacks.o
```

Alternatively, add it to your object file list:
```makefile
OBJS = your_source.o \
       $(IP2CO_INC)/IpoptToConoptCallbacks.o
```

The Makefiles in the examples demonstrate this pattern. Note that the bridge is **not** a header-only library - the `IpoptToConoptCallbacks.cpp` file contains the actual implementation of the callback functions and must be compiled.

#### Step 3: Compile Your Source Files

Compile your own source files with the same include paths:
```bash
g++ -c your_source.cpp \
    -I$(CONOPT_DIR)/include \
    -I$(IP2CO_INC) \
    -I$(IP2CO_INC)/Ipopt \
    -I$(IPOPT_DIR)/include/coin-or \
    -o your_source.o
```

#### Step 4: Link Everything Together

Link all object files (including `IpoptToConoptCallbacks.o`) against CONOPT and IPOPT libraries:
```bash
g++ -o your_program \
    your_source.o \
    IpoptToConoptCallbacks.o \
    -L$(CONOPT_DIR)/lib -L$(IPOPT_DIR)/lib \
    -lconopt -lipopt -lstdc++ -lm -ldl \
    -Wl,--rpath,$(CONOPT_DIR)/lib -Wl,--rpath,$(IPOPT_DIR)/lib
```

**Important**: The `IpoptToConoptCallbacks.cpp` file is required - without it, the linker will fail with undefined references to the CONOPT callback functions.

## License

See `LICENSE` file for license information.

## Contributing

This bridge is designed to maintain compatibility with Ipopt's interface while leveraging CONOPT's solver capabilities. When contributing:

- Maintain API compatibility with Ipopt's `TNLP` interface
- Ensure proper status code mappings
- Test with the provided examples
- Document any constraint transformations or limitations
