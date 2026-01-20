cmake_minimum_required(VERSION 3.16) # Bumped for better LTO support
project(agl_cpp)

# 1. Force Release Build
# This ensures assertions are disabled and optimization flags are active.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# 2. Interprocedural Optimization (LTO/IPO)
# This is the most critical optimization for C++/Python bindings. 
# It allows the compiler to inline code across different .cpp files.
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()

# --- PACKAGES ---
find_package(pybind11 REQUIRED)
find_package(Eigen3 3.3 REQUIRED) # Specify version if possible
find_package(OpenMP REQUIRED)

# Manual Chemfiles Search (Conda optimized)
find_library(CHEMFILES_LIB_FILE 
    NAMES chemfiles 
    HINTS $ENV{CONDA_PREFIX}/lib 
    NO_DEFAULT_PATH
)
# Fallback
if(NOT CHEMFILES_LIB_FILE)
    find_library(CHEMFILES_LIB_FILE NAMES chemfiles)
endif()
if(NOT CHEMFILES_LIB_FILE)
    message(FATAL_ERROR "Could not find chemfiles library.")
endif()

# --- MODULE ---
# 'THIN_LTO' creates smaller binaries and compiles faster while keeping speed
pybind11_add_module(${PROJECT_NAME} agl_binding.cpp)

# --- OPTIMIZATION FLAGS ---

# 3. Compiler Specific Optimizations
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE
        /O2     # Maximize speed
        /Ob2    # Aggressive inlining
        /Oi     # Use intrinsic functions
        /Ot     # Favor fast code
        /fp:fast # Fast floating point model
        /MP     # Multiprocessor compilation (build speed)
    )
else()
    target_compile_options(${PROJECT_NAME} PRIVATE
        -O3             # Maximize optimization level
        -march=native   # Optimize for the CPU found on THIS machine
        -fPIC           # Position Independent Code
        
        # Aggressive Math Optimizations (Use with caution if high precision is required)
        -ffast-math     
        -funroll-loops  
        
        # Visibility hidden reduces binary size and load time for Python modules
        -fvisibility=hidden
    )
endif()

# 4. Critical Definitions for Eigen
# Without these, Eigen performs costly range checks and debug operations.
target_compile_definitions(${PROJECT_NAME} PRIVATE
    NDEBUG
    EIGEN_NO_DEBUG
)

# --- INCLUDE & LINK ---
target_include_directories(${PROJECT_NAME} PRIVATE 
    $ENV{CONDA_PREFIX}/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${CHEMFILES_LIB_FILE}
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
)