cmake_minimum_required(VERSION 3.25)
project(hydroisomerization_v1 VERSION 0.1.0 LANGUAGES C Fortran)

option(HYDROISO_ENABLE_TESTS "Build V1 tests" ON)
find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecserial sunmatrixdense sunlinsoldense)

set(HYDROISO_CORE_SOURCES
  src/kinds.f90
  src/constants.f90
  src/model_types.f90
  src/chemistry/species.f90
  src/chemistry/stoichiometry.f90
  src/chemistry/arrhenius.f90
  src/chemistry/kinetics.f90
  src/transport/knudsen.f90
  src/numerics/fvm_1d.f90
  src/pellet/pellet_1d.f90
  src/pellet/effectiveness.f90
  src/io/config.f90
  src/io/output.f90
)

add_library(hydroiso_core STATIC ${HYDROISO_CORE_SOURCES})
set_target_properties(hydroiso_core PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)
target_compile_options(hydroiso_core PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-Wall;-Wextra;-Wimplicit-interface;-fcheck=all;-fbacktrace>
)

add_library(hydroiso_cvode STATIC
  src/numerics/solver_api.f90
  src/numerics/cvode_backend.f90
)
set_target_properties(hydroiso_cvode PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)
target_include_directories(hydroiso_cvode PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_link_libraries(hydroiso_cvode PUBLIC
  hydroiso_core
  SUNDIALS::fcvode_mod_shared
  SUNDIALS::fnvecserial_mod_shared
  SUNDIALS::fsunmatrixdense_mod_shared
  SUNDIALS::fsunlinsoldense_mod_shared
)

add_executable(hydroiso_v1 src/main.f90)
target_include_directories(hydroiso_v1 PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_link_libraries(hydroiso_v1 PRIVATE hydroiso_cvode)

if(HYDROISO_ENABLE_TESTS)
  enable_testing()
  add_executable(test_core tests/fortran/test_core.f90)
  target_include_directories(test_core PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")
  target_link_libraries(test_core PRIVATE hydroiso_core)
  add_test(NAME core_formulas COMMAND test_core)
endif()
