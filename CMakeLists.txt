cmake_minimum_required(VERSION 3.19.3)
project(FPE)

# Use C++14 standard to make it compatible with BioDynaMo later on.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

## ------------------------------ CNPY -----------------------------------------
## Since there is not support for the AAA algorithm, we use python to compute 
## the expansion coefficients c_\inf, c_k, and d_k. We export them in the `npy`
## format and use CNPY to read them from the CPP part or the program.

set(CNPY_PATH "/usr/local/include") # PATH to source directore, e.g. cloned git
set(CNPY_LIB "/usr/local/lib") # PATH to installed directory
if(APPLE)
  set(CNPY_LIB_NAME "libcnpy.dylib")
else()
  set(CNPY_LIB_NAME "libcnpy.so")
endif(APPLE)
include_directories(${CNPY_PATH})
add_library(cnpy SHARED IMPORTED)
set_target_properties(cnpy PROPERTIES
  IMPORTED_LOCATION "${CNPY_LIB}/${CNPY_LIB_NAME}"
  INTERFACE_INCLUDE_DIRECTORIES "${CNPY_PATH}"
)


## -------------------------- USE MFEM -----------------------------------------
## MFEM is used for the finite element discretization, computation of mass and 
## stiffness matrix, and general sparse matrix-vector operations. The cmake 
## code of this section was provided by MFEM / LLNL directly.

# Use MFEM in this project.
# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
  find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
  "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
  find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
message(STATUS "Mfem compiler is: ${MFEM_CXX_COMPILER}")
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

message(STATUS "Mfem include dir is: ${MFEM_INCLUDE_DIRS}")
include_directories(${MFEM_INCLUDE_DIRS})
message(STATUS "Mfem library is: ${MFEM_LIBRARIES}")

## ---------------------------  Compile Main  ----------------------------------

# Consider all files in src/ for diffusion simulation.
include_directories("src")
file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cc)
message(STATUS "Building diffusion simulation. Considering the follwing: \n" 
        "  .cc: ${PROJECT_SOURCES} \n" "  .h: ${PROJECT_HEADERS}")
# define target executable for diffusion simulation
add_executable(
  ${CMAKE_PROJECT_NAME}
  ${PROJECT_HEADERS} ${PROJECT_SOURCES}
)
# Link executable against mfem
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${MFEM_LIBRARIES} cnpy)
