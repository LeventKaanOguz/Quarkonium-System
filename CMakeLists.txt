cmake_minimum_required(VERSION 3.30)
project(QuarkoniumGEM VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(QGEM_SANITIZE "Build with ASan + UBSan" OFF)
# CERN ROOT (Minuit2) is optional: AUTO uses it when found, ON requires it.
set(QGEM_WITH_ROOT AUTO CACHE STRING "Use CERN ROOT/Minuit2: AUTO, ON or OFF")
set_property(CACHE QGEM_WITH_ROOT PROPERTY STRINGS AUTO ON OFF)

# No version argument: Eigen 3.4 and 5.x both satisfy it (a "3.4" request rejects 5.x).
find_package(Eigen3 5.0 REQUIRED NO_MODULE)
find_package(GTest REQUIRED)
find_package(SQLite3 REQUIRED)  # target SQLite3::SQLite3 needs CMake >= 4.1

add_compile_options(-Wall -Wextra -Wpedantic)
# GCC reports -Wmaybe-uninitialized inside Eigen's lazy product evaluators
# (X^T A X); every frame of the report is Eigen-internal: a known false
# positive.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-Wno-maybe-uninitialized)
endif()
if(QGEM_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# Library: every src/**/*.cpp except apps. Header-only until the first .cpp exists.
file(GLOB_RECURSE QGEM_SOURCES CONFIGURE_DEPENDS src/*.cpp)
list(FILTER QGEM_SOURCES EXCLUDE REGEX "src/(apps|root)/")
if(QGEM_SOURCES)
  add_library(qgem STATIC ${QGEM_SOURCES})
  target_include_directories(qgem PUBLIC src)
  target_link_libraries(qgem PUBLIC Eigen3::Eigen SQLite3::SQLite3)
else()
  add_library(qgem INTERFACE)
  target_include_directories(qgem INTERFACE src)
  target_link_libraries(qgem INTERFACE Eigen3::Eigen SQLite3::SQLite3)
endif()

add_executable(qgem_cli src/apps/qgem.cpp)
set_target_properties(qgem_cli PROPERTIES OUTPUT_NAME qgem)
target_link_libraries(qgem_cli PRIVATE qgem)

# Optional ROOT bridge: only src/root/ and the CLI see ROOT headers (Minuit2 is
# plain C++; ROOT itself is built as C++20, the project is C++23).
if(NOT QGEM_WITH_ROOT STREQUAL "OFF" AND NOT QGEM_SANITIZE)
  if(QGEM_WITH_ROOT STREQUAL "ON")
    find_package(ROOT 6.30 CONFIG REQUIRED COMPONENTS Minuit2)
  else()
    find_package(ROOT 6.30 CONFIG QUIET COMPONENTS Minuit2)
  endif()
endif()
if(ROOT_FOUND)
  message(STATUS "CERN ROOT ${ROOT_VERSION}: building the Minuit2 bridge")
  add_library(qgem_root STATIC src/root/minuit2_engine.cpp)
  target_link_libraries(qgem_root PUBLIC qgem ROOT::Minuit2)
  target_compile_definitions(qgem_root PUBLIC QGEM_HAVE_ROOT)
  target_link_libraries(qgem_cli PRIVATE qgem_root)
else()
  message(STATUS "CERN ROOT not used (QGEM_WITH_ROOT=${QGEM_WITH_ROOT})")
endif()
target_compile_definitions(qgem_cli PRIVATE QGEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

enable_testing()
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS tests/test_*.cpp)
add_executable(qgem_tests ${TEST_SOURCES})
target_link_libraries(qgem_tests PRIVATE qgem GTest::gtest_main)
target_compile_definitions(qgem_tests PRIVATE QGEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
include(GoogleTest)
gtest_discover_tests(qgem_tests)

if(ROOT_FOUND)
  add_executable(qgem_root_tests tests/root/test_minuit2.cpp)
  target_link_libraries(qgem_root_tests PRIVATE qgem_root GTest::gtest_main)
  target_compile_definitions(qgem_root_tests PRIVATE QGEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  gtest_discover_tests(qgem_root_tests)
endif()

# The analysis statistics (model_selection.py) when Python 3 with numpy and
# scipy is available.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
                  RESULT_VARIABLE QGEM_PY_NUMPY OUTPUT_QUIET ERROR_QUIET)
  if(QGEM_PY_NUMPY EQUAL 0)
    add_test(NAME analysis_model_selection
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_SOURCE_DIR}/analysis/test_model_selection.py)
  endif()
endif()

# Slow physics validation of the statistics (bootstrap, MCMC, exact leave-one-out refits): ctest label "slow".
file(GLOB VALIDATION_SOURCES CONFIGURE_DEPENDS tests/validation/*.cpp)
if(VALIDATION_SOURCES)
  add_executable(qgem_validation ${VALIDATION_SOURCES})
  target_link_libraries(qgem_validation PRIVATE qgem GTest::gtest_main)
  target_compile_definitions(qgem_validation PRIVATE QGEM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  add_test(NAME validation COMMAND qgem_validation)
  set_tests_properties(validation PROPERTIES LABELS slow TIMEOUT 1800)
endif()
